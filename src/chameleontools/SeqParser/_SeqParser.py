"""
TABI - UCSC iGEM 2023

StealthGenome: Prep Genome for Stealth and Codon Frequency analysis

PlasmidParser: Identifies Stealth-mutable regions of an annotated plasmid

Author: Tyler Gaw
"""
# Usage to be written, don't wanna rn

from Bio import SeqIO
from Bio.SeqFeature import SimpleLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Iterator,Literal
from chameleontools.ORFfinder import ORFfinder
from chameleontools.FastAreader import FastAreader
import sys


class StealthGenome:
    """
    Reads genome file in genbank or FASTA format,

    Stores Genome sequence and CDS regions
    """

    def _validArg(self, genome: str):
        "helper func: checks valid args"
        g_valid = False, False, False
        g_valid = (
            genome.endswith(".fasta")
            or genome.endswith(".gb")
            or genome.endswith(".gbk")
        )
        return g_valid

    def _fastaGenome(self, fasta: str):
        reader = FastAreader(fasta)
        CDS_records = []
        sequence = []
        for frag, (_, seq) in enumerate(reader.readFasta()):
            genome = Seq(seq)
            orfs = ORFfinder(
                seq, longestGene=True, minGene=200
            ).get_genes()  # array of estimated genes
            for i, j in enumerate(orfs):
                "Converting ORFfinder output to Biopython SimpleLocation"
                loc = SimpleLocation(j[0] - 1, j[1], strand=1 if j[3] > 0 else -1)

                "Building Biopython SeqRecord for CDS"
                rec = SeqRecord(loc.extract(genome))
                rec.description = "Predicted Protein Coding Sequence"
                rec.name = "Predicted CDS HEADER"
                rec.id = f"FRAG_{frag} CDS_{i} for CODON STATS"

                CDS_records.append(rec)
            sequence.append(genome)
        return CDS_records, sequence

    ############################################################################################
    # ^^^ Private Helper Functions ^^^
    ############################################################################################

    def __init__(self, genome_infile) -> None:
        """
        Constructor: Validates inputs, writes to files in FASTA format, saves file names

        Can be used for Stealth_CDS_Extract to produce a CDS generator for input file

        Args:
        genome_infile (str): genome file in genbank or FastA format

        Class Attr.:

        self.getGenome() -> array of Seq() 

        self.getCDS() -> Iterator of Biopython SeqRecord for each CDS | list of SeqRecord if using FastA file
        
        self.filetype() -> type of input file | "fasta" if FastA file, "genbank" if GenBank file
        """
        if genome_infile != None:
            g_bool = self._validArg(genome_infile)

            if not g_bool:
                print(
                    f"Genome File needs to be FASTA file (.fasta) or GenBank file (.gb || .gbk). Got {genome_infile}",
                    file=sys.stderr,
                )
                exit()

            if genome_infile.endswith(".fasta"):
                self.cds_sequence, self.genome_sequence = self._fastaGenome(
                    genome_infile
                )
                self.input = "fasta"
            else:
                self.cds_sequence = self.Stealth_CDS_Extract(genome_infile)
                self.genome_sequence = self.genome_sequence
                self.input = "genbank"

    def Stealth_CDS_Extract(self, gbfile: str) -> Iterator[SeqRecord]:
        """
        Reads a genbank file and yields CDS/ORF features for Codon Optimization

        A generator that yields CDS features in FASTA format

        Used to produce data for codon usage statistics
        """
        gb = SeqIO.parse(gbfile, "genbank")
        count = 0
        for rec in gb:
            self.genome_sequence = [rec.seq]  # define self.genome_sequence from genbank record
            for feature in rec.features:
                if feature.type == "CDS" or feature.type == "ORF":
                    count += 1
                    out = feature.extract(rec)

                    """Set the description to 'product' qualifier if avaiable. Else use 'label' qualifier or [null] if labels are not present"""
                    out.description = feature.qualifiers.get(
                        "product", feature.qualifiers.get("label", ["null"])
                    )[0]

                    out.name = "TEMP CDS HEADER"
                    out.id = f"CDS_{count} CODON STATS"
                    yield out
            return  # Only accepts single entry genbank records

    def filetype(self) -> Literal["fasta","genbank"]:
        return self.input
    
    def getGenome(self) -> list[Seq]:
        return self.genome_sequence
    
    def getCDS(self) -> Iterator[SeqRecord] | list[SeqRecord]:
        return self.cds_sequence
        
class PlasmidParse:
    """
    Reads in an annotated plasmid file in genbank format,
    parses mutable regions from CDS/ORF annotations,
    takes statistics for total mutable range
    """

    def _validArg(self, plasmid_infile: str) -> bool:
        "helper func: checks for valid input"
        return plasmid_infile.endswith(".gb") or plasmid_infile.endswith(".gbk")

    def _parsePlasmid(
        self, gbfile: str
    ) -> tuple[list[SimpleLocation], list[SimpleLocation]]:
        "helper func: Filters for CDS and noncoding sequence"
        avoid, cds = [], []
        gb = SeqIO.parse(gbfile, "genbank")
        for rec in gb:
            self.record = rec
            for feature in rec.features:
                if feature.type in {"source", "gene"}:
                    continue
                if feature.type in {"CDS", "ORF"}:
                    cds.append((feature.location, feature.location.start % 3))
                else:
                    avoid.append((feature.location, feature.location.start % 3))
            break  # Only handles single entry genbank records
        return avoid, cds

    def _trimStart(
        self, loc_list: list[SimpleLocation], cds: SimpleLocation
    ) -> list[SimpleLocation]:
        "helper func: trims the start codon of CDS, handles fwd/rev strands"
        ref = (
            SimpleLocation(cds.start, cds.start + 3, cds.strand)
            if cds.strand > 0
            else SimpleLocation(cds.end - 3, cds.end, cds.strand)
        )
        for i in range(len(loc_list)):
            loc = loc_list[i]
            if loc.start < ref.start and loc.end > ref.end:
                loc_list[i] = SimpleLocation(loc.start, ref.start, loc.strand)
                loc_list.insert(i + 1, SimpleLocation(ref.end, loc.end, loc.strand))
                break
            if loc.start in ref:
                loc_list[i] = SimpleLocation(ref.end, loc.end, loc.strand)
                break
            elif loc.end in ref or loc.end == ref.end:
                loc_list[i] = SimpleLocation(loc.start, ref.start, loc.strand)
                break
        return

    def _frameAdj(self, loc_list: list[SimpleLocation], frame: int) -> SimpleLocation:
        "helper func: trims mutable regions to lie in frame"
        new = []
        for loc in loc_list:
            st = loc.start
            ed = loc.end
            for _ in range(1, 3):
                if st % 3 != frame:
                    st += 1
                if ed % 3 != frame:
                    ed -= 1
            if (ed - st) < 3:
                continue
            adj = SimpleLocation(st, ed, loc.strand)
            new.append(adj)
        return new

    def _removeOverlap(
        self, loc_list: list[SimpleLocation], cds: SimpleLocation
    ) -> list[SimpleLocation]:
        "helper func: removes overlapping cds regions"
        ret = []
        for loc in loc_list:
            new = [[loc.start, loc.end]]
            if loc.start in cds and (loc.end in cds or loc.end == cds.end):
                "completely removes region"
                continue
            if loc.start < cds.start and loc.end > cds.end:
                "splits range into two"
                dummy = [cds.end, loc.end]
                new[-1][1] = cds.start
                new.append(dummy)
            elif loc.start in cds:
                "truncates start"
                new[-1][0] = cds.end
            elif loc.end in cds or loc.end == cds.end:
                "truncates end"
                new[-1][1] = cds.start + 1
            ret.extend([SimpleLocation(i, j, loc.strand) for i, j in new])
        return ret

    def _windowCorrection(
        self,
        old: list[list[SimpleLocation], int, SimpleLocation],
        st_bound: int,
        ed_bound: int,
    ) -> list[list[SimpleLocation], int, SimpleLocation]:
        "helper func: corrects for out-frame overlapping CDS"
        loc, fm, parent = old
        new = SimpleLocation(st_bound, ed_bound, parent.strand)
        dummy = []

        for i in loc:
            if st_bound > i.end:
                continue
            if i.start in new and (i.end in new or i.end == new.end):
                dummy.append(i)
            elif i.start in new:
                dummy.append(SimpleLocation(i.start, ed_bound, i.strand))
                continue
            elif i.end in new or i.end == new.end:
                dummy.append(SimpleLocation(st_bound, i.end, i.strand))
                continue

        return [dummy, fm, new]

    ############################################################################################
    # Private Helper Functions
    ############################################################################################

    def __init__(self, plasmid_infile: str) -> None:
        """
        Constructor: Validates inputs, parses mutable regions

        Args:
        plasmid_infiled (str): genome file in FASTA or genbank format

        Class Methods:
        
        self.getGenBank() -> Biopython SeqRecord() of input plasmid, contains feature annotations

        self.getSeq() -> Biopython Seq() of input Plasmid

        self.regions() -> array of Biopython SimpleLocation() ranges of mutable codons

        self.mutableCount() -> final mutable regions on plasmid

        self.unmutableCount() -> number of removed basepairs from total CDS coverage

        self.regionCount() -> Total CDS coverage before extracting mutable regions
        """
        p_bool = self._validArg(plasmid_infile)
        if not p_bool:
            print(
                f"Plasmid Input File needs to be GenBank file (.gb || .gbk). Got: {plasmid_infile}",
                file=sys.stderr,
            )
            exit()

        "Did this cause it was ugly when using tuple unpacking"
        temp_parse = self._parsePlasmid(plasmid_infile)

        temp = self.defineMutable(temp_parse)
        
        self.record = self.record 
        
        self.plasmid_sequence = self.record.seq # self.record -> Bio.SeqRecord of plasmid

        self.mutable_regions = temp[0]
        self.total_mutable = temp[1]
        self.removed = temp[2]

        self.total_coverage = self.total_mutable + self.removed

    def defineMutable(
        self, location: tuple[list[SimpleLocation]]
    ) -> list[SimpleLocation]:
        "Finds non-overlapped mutable regions in an annotated plasmid"

        avoid, cdsRegions = location
        temp = []  # [[[SimpLocation List],frame,parent CDS]]

        "Removes mutable cds regions from overlapping non-coding sequence"
        for cds, frame in cdsRegions:
            if len(cds) <= 3:
                continue
            bounds = [[cds.start, cds.end]]
            for loc, _ in avoid:
                if cds.start > loc.start and cds.end < loc.end:
                    "case: interior ORF to non-coding region"
                    bounds = [[-1, -1]]
                    break
                if bounds[-1][0] < loc.start and bounds[-1][1] > loc.end:
                    "case: non-coding region interior to ORF"
                    new_bound = [loc.end, cds.end]
                    bounds[-1][1] = loc.start
                    bounds.append(new_bound)
                elif bounds[-1][0] in loc:
                    "case: clipped start"
                    bounds[-1][0] = loc.end
                elif bounds[-1][1] in loc or bounds[-1][1] == loc.end:
                    "case: clipped end"
                    bounds[-1][1] = loc.start
            if (bounds[-1][1] - bounds[-1][0]) <= 3:
                del bounds[-1]
            if len(bounds) <= 0:
                continue
            mut_range = [
                SimpleLocation(start=rng[0], end=rng[1], strand=cds.strand)
                for rng in bounds
            ]
            temp.append([mut_range, frame, cds])

        "Trimming start codons"
        for i in range(len(temp)):
            tLoc = temp[i][0]
            frm = temp[i][1]
            reg = temp[i][2]
            window = temp[i + 1 :]
            for j in range(len(window)):
                _, parent_frame, del_cds = window[j]
                if (del_cds.start in reg and del_cds.strand > 0) or (
                    del_cds.end in reg and del_cds.strand < 0
                ):
                    self._trimStart(temp[i][0], del_cds)
                    tLoc = temp[i][0]
            if reg.strand > 0:
                trimmed = SimpleLocation(reg.start + 3, reg.end, reg.strand)
                if tLoc[0].start == reg.start:
                    temp[i][0][0] = SimpleLocation(
                        trimmed.start, tLoc[0].end, tLoc[0].strand
                    )
            else:
                trimmed = SimpleLocation(reg.start, reg.end - 3, reg.strand)
                if tLoc[-1].end == reg.end:
                    temp[i][0][-1] = SimpleLocation(
                        tLoc[-1].start, trimmed.end, tLoc[-1].strand
                    )

        "sorts regions by CDS length, compares larger CDS to smaller CDS to remove overlap"
        temp = sorted(temp, key=lambda x: len(x[2]), reverse=True)

        "remove overlaps"
        for i in range(len(temp)):
            frag = temp[i][0]
            frm = temp[i][1]
            region = temp[i][2]
            window = temp[i + 1 :]
            if region == None:
                continue
            for j in range(len(window)):
                _, parent_frame, parent_cds = window[j]
                if parent_cds == None:
                    continue
                if parent_cds.start in region or parent_cds.end in region:
                    if parent_frame != frm:
                        temp[i][0] = self._removeOverlap(frag, parent_cds)
                        frag = temp[i][0]
                if parent_cds.start in region and parent_cds.end in region:
                    temp[i + j + 1][2] = None
                elif parent_cds.start in region:
                    temp[i + j + 1] = self._windowCorrection(
                        window[j], region.end, parent_cds.end
                    )
                elif parent_cds.end in region:
                    temp[i + j + 1] = self._windowCorrection(
                        window[j], parent_cds.start, region.start
                    )

        "frame adjusts mutable fragments to lie in frame of parent CDS"
        for i in range(len(temp)):
            temp[i][0] = self._frameAdj(temp[i][0], temp[i][1])

        mutable_regions = sorted(
            [region for CDS, _, tf in temp for region in CDS if tf is not None],
            key=lambda x: x.start,
        )

        "mutable region statistics - total covered region and mutable regions"
        total_coverage = sum([len(x[-1]) for x in temp if x[-1] is not None])
        removed = total_coverage - sum([len(x) for x in mutable_regions])
        return mutable_regions, total_coverage, removed

    def getSeq(self) -> Seq:
        return self.plasmid_sequence
    
    def regions(self) -> list[SeqRecord]:
        return self.mutable_regions
    
    def mutableCount(self) -> int:
        return self.total_mutable
    
    def regionCount(self) -> int:
        return self.total_coverage
    
    def unmutableCount(self) -> int:
        return self.removed 
    
    def getGenBank(self) -> SeqRecord:
        return self.record
    

__all__ = ["StealthGenome","PlasmidParse"]
from Bio import SeqIO
from Bio.SeqFeature import SimpleLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Iterator, Literal
from chameleontools.ORFfinder import ORFfinder
from chameleontools.FastAreader import FastAreader
import sys


class StealthGenome:
    """
    Reads genome file in genbank or FASTA format,

    Stores Genome sequence and CDS regions
    """

    def _validArg(self, genome: str):
        """"Helper method, checks valid args"

        Args:
            genome (str): input file name

        Returns:
            bool: valid input
        """
        g_valid = False, False, False
        g_valid = (
            genome.endswith(".fasta")
            or genome.endswith(".fa")
            or genome.endswith(".gb")
            or genome.endswith(".gbk")
        )
        return g_valid

    def _fastaGenome(self, fasta: str):
        """Helper method, handles reading in of FastA files

        Args:
            fasta (str): FastA filename

        Returns:
            tuple[list[Bio.SeqRecord.SeqRecord],Bio.Seq.Seq]: List of CDS regions and genome sequence
        """
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

    def __init__(self, genome_infile: str) -> None:
        """Validates inputs and stores genome sequence and CDS regions

        Args:
            genome_infile (str): genome file in GenBank or FastA format
            
        Usage:

            self.getGenome() -> array of Seq()

            self.getCDS() -> List of SeqRecord

            self.filetype() -> type of input file | "fasta" if FastA file, "genbank" if GenBank file
        """
        if genome_infile != None:
            g_bool = self._validArg(genome_infile)

            if not g_bool:
                print(
                    f"Genome File needs to be FASTA file (.fasta || .fa) or GenBank file (.gb || .gbk). Got {genome_infile}",
                    file=sys.stderr,
                )
                exit()

            if genome_infile.endswith(".fasta") or genome_infile.endswith(".fa"):
                self.cds_sequence, self.genome_sequence = self._fastaGenome(
                    genome_infile
                )
                self.input = "fasta"
            else:
                self.cds_sequence,self.genome_sequence = self._Stealth_CDS_Extract(genome_infile)
                self.input = "genbank"

    def _Stealth_CDS_Extract(self, gbfile: str) -> tuple[list[Seq],Iterator[SeqRecord]]:
        """Generator for CDSs in a GenBank file
        
        Args:
            gbfile (str): input GenBank file
        """
        gb = SeqIO.parse(gbfile, "genbank") 
        for rec in gb:
            genome_sequence = [
                rec.seq
            ]  # define self.genome_sequence from genbank record
            def _generator() -> Iterator[SeqRecord]:
                count = 0
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
            return  _generator(), genome_sequence  # Only accepts single entry genbank records

    def filetype(self) -> Literal["fasta", "genbank"]:
        """Returns filetype of input file

        Returns:
            Literal[str]: "fasta" or "genbank"
        """
        return self.input

    def getGenome(self) -> list[Seq]:
        """Returns genome sequence(s)

        Returns:
            list[Seq]: Sequence(s) of genome stored in a list
        """
        return self.genome_sequence

    def getCDS(self) -> list[SeqRecord]:
        """Returns iterable of CDS regions

        Returns:
            list[SeqRecord]: list of CDSs

        Yields:
            Iterator[SeqRecord]: Iterable of CDSs
        """
        return self.cds_sequence


class PlasmidParse:
    """Reads in an annotated plasmid file in genbank format,
    parses mutable regions from CDS/ORF annotations,
    takes statistics for total mutable range
    """

    def _validArg(self, plasmid_infile: str) -> bool:
        """Helper method, checks for valid input
        
        Args:
            plasmid_infile (str): filename for input file
            
        Returns:
            bool: boolean of valid filename """
        return plasmid_infile.endswith(".gb") or plasmid_infile.endswith(".gbk")

    def _parsePlasmid(
        self, gbfile: str
    ) -> tuple[list[SimpleLocation], list[SimpleLocation]]:
        """Helper func, filters for CDS and noncoding sequence

        Args:
            gbfile (str): input plasmid sequence

        Returns:
            tuple[list[SimpleLocation], list[SimpleLocation]]: regions to avoid and CDS regions
        """
        avoid, cds = set(),set() # sets to avoid duplicate annotation regions
        gb = SeqIO.parse(gbfile, "genbank")
        for rec in gb:
            self.record = rec
            for feature in rec.features:
                '''SimpleLocation unhashable, break down into data and recreate object'''
                simpLoc = feature.location
                simpLoc_data = (simpLoc.start,simpLoc.end,simpLoc.strand)
                feat_data = (simpLoc_data, simpLoc.start % 3) # keeps location, frame
                if feature.type in self.ignore_annotations:
                    continue
                if feature.type in self.keep_annotations: 
                    cds.add(feat_data)
                else:
                    avoid.add(feat_data)
            break  # Only handles single entry genbank records
        
        '''recreate SimpleLocation objects'''
        avoid = [(SimpleLocation(x[0],x[1],x[2]),y) for x,y in avoid]
        cds = [(SimpleLocation(x[0],x[1],x[2]),y) for x,y in cds] 
        return avoid, cds

    def _trimStart(
        self, loc_list: list[SimpleLocation], cds: SimpleLocation
    ) -> list[SimpleLocation]:
        """Helper method, trims the start codon of CDS from others that overlap, handles fwd/rev strands

        Args:
            loc_list (list[SimpleLocation]): list of SimpleLocation objects to subtract from
            cds (SimpleLocation): CDS to subtract

        Returns:
            list[SimpleLocation]: list of trimmed SimpleLocation objects
        """
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

    def _frameAdj(self, loc_list: list[SimpleLocation], frame: int) -> list[SimpleLocation]:
        """Helper func, trims mutable regions to lie in frame

        Args:
            loc_list (list[SimpleLocation]): List of SimpleLocations to frame adjust
            frame (int): Frame to adjust to

        Returns:
            List[SimpleLocation]: List of frame-adjusted SimpleLocations
        """
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
        """Helper method, removes overlapping cds regions

        Args:
            loc_list (list[SimpleLocation]): List of locations to subtract from
            cds (SimpleLocation): CDS to subtract 

        Returns:
            list[SimpleLocation]: List of non-overlapping regions of CDS
        """
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
        """Helper method, corrects for out-frame overlapping CDS

        Args:
            old (list[list[SimpleLocation], int, SimpleLocation]): list of current mutable regions
            st_bound (int): beginning of overlapping CDS
            ed_bound (int): end of overlapping CDS

        Returns:
            list[list[SimpleLocation], int, SimpleLocation]: New list of mutable regions
        """
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

    def __init__(self, plasmid_infile: str, keep_annotation = [], ignore_annotation = []) -> None:
        """Validates inputs, parses out mutable regions

        Args:
            plasmid_infiled (str): plasmid file in FASTA or genbank format

        Usage:
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

        self.keep_annotations = {"CDS","ORF","gene"}
        self.ignore_annotations = {"source"}
        
        for annotation in keep_annotation:
            self.keep_annotations.add(annotation)
        for annotation in ignore_annotation:
            self.ignore_annotations.add(annotation)
        
        "Did this cause it was ugly when using tuple unpacking"
        temp_parse = self._parsePlasmid(plasmid_infile)

        temp = self.defineMutable(temp_parse)

        self.record = self.record

        self.plasmid_sequence = (
            self.record.seq
        )  # self.record -> Bio.SeqRecord of plasmid

        self.mutable_regions = temp[0]
        self.total_coverage = temp[1]
        self.removed = temp[2]

        self.total_mutable = self.total_coverage - self.removed

    def defineMutable(
        self, location: tuple[list[SimpleLocation]]
    ) -> tuple[list[SimpleLocation],int,int]:
        """
        Parses out all mutable regions based on plasmid annotations
        """
        
        
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
                    if parent_frame != frm or parent_cds.strand != region.strand:
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
        """Get the sequence of the plasmid

        Returns:
            Seq: Plasmid NT sequence
        """
        return self.plasmid_sequence

    def regions(self) -> list[SimpleLocation]:
        """Get the list of mutable regions

        Returns:
            list[SimpleLocation]: SimpleLocations defining mutable regions
        """
        return self.mutable_regions

    def mutableCount(self) -> int:
        """Get length of mutable regions

        Returns:
            int: NT count of mutable regions
        """
        return self.total_mutable

    def regionCount(self) -> int:
        """Get total length of CDS regions

        Returns:
            int: NT count of CDS regions
        """
        return self.total_coverage

    def unmutableCount(self) -> int:
        """Get length of unremoved regions

        Returns:
            int: NT count of unmutable CDS regions
        """
        return self.removed

    def getGenBank(self) -> SeqRecord:
        """Get GenBank record of plasmid infile

        Returns:
            SeqRecord: SeqRecord of input plasmid
        """
        return self.record


__all__ = ["StealthGenome", "PlasmidParse"]

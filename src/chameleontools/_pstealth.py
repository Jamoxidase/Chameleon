import chameleontools.SeqParser as SeqParser
import chameleontools.StealthParser as StealthParser
import chameleontools.Stealth as Stealth
import chameleontools.CodonAnalyzer as CodonAnalyzer
import chameleontools.ChromatoSeq as ChromatoSeq
import Bio.Seq as Seq, Bio.SeqFeature as SeqFeature
import Bio.SeqIO as SeqIO
import sys
from tqdm import tqdm


class pipeline:
    def _openOutfile(self, outfile: str):
        if outfile == None:
            return sys.stdout
        else:
            return open(outfile, "w")

    def __init__(
        self,
        genome_infile: str,
        plasmid_infile: str,
        outfile: str,
        z_score: float,
        pseudo: int,
        kMax: int,
        kMin: int,
        palindrome: bool = False,
        keep: list[str] = [],
        ignore: list[str] = []
    ) -> None:
        self.loading_bar = tqdm(total=100, desc="Reading Files")  # CLI loading bar

        self.stealth_input = SeqParser.StealthGenome(genome_infile)
        self.plasmid_input = SeqParser.PlasmidParse(plasmid_infile, keep_annotation= keep, ignore_annotation= ignore)

        self.loading_bar.update(15)

        self.stealth_motifs = self.load_stealth(z_score, pseudo, kMax, kMin, palindrome)

        self.unremoved_motifs = 0

        self.total_motifs = 0

        self.run(outfile)

    def load_stealth(
        self,
        z_score: float,
        pseudo: int,
        kMax: int,
        kMin: int,
        palindrome: bool = False,
    ):
        self.loading_bar.set_description("Performing Stealth Analysis")

        stealth_out = Stealth.StealthV0(self.stealth_input, z_score, pseudo, kMax, kMin)
        self.loading_bar.update(10)
        return StealthParser.ParseStealth(stealth_out, palindrome)

    def codon_analysis(self) -> CodonAnalyzer.CDSanalyzer:
        self.loading_bar.set_description("Gathering Codon Statistics")

        codon_usage = CodonAnalyzer.CDSanalyzer(self.stealth_input)
        self.loading_bar.update(20)
        return codon_usage

    def hide_motif(self, region: list[SeqFeature.SimpleLocation]):
        plasmid_seq = self.plasmid_input.getSeq()
        codon_freq = self.codon_analysis().getFrequency()
        codon_optimizer = CodonAnalyzer.CodonOptimizer(codon_freq)
        pattern_remover = ChromatoSeq.patternConstrainer(
            self.stealth_motifs, codon_freq
        )
        motif_check = ChromatoSeq.MotifChecker(self.stealth_motifs)

        def _hide_motif_helper(aaSeq: str) -> Seq.Seq:
            roundsData = []
            for _ in range(5):
                seed = codon_optimizer.assembleSeed(aaSeq)
                optimization = pattern_remover.optimizeSequence(seed)

                mortality = motif_check.checkMotifs(optimization)
                roundsData.append((mortality, optimization))
            lowest_mortality = min(roundsData, key=lambda x: x[0])
            self.unremoved_motifs += lowest_mortality[0]
            best_seq = lowest_mortality[1]
            return Seq.Seq(best_seq)

        self.loading_bar.set_description("Removing Motifs")

        modified_blocks = []  # list of SimpleLocations, mapped to modified Seq()
        
        loading_progress = 48 / (sum([len(x) for x in region])) # more dynamic loading bar for motif removal
        
        for s_loc in region:
            coding_seq = s_loc.extract(plasmid_seq)
            self.total_motifs += motif_check.checkMotifs(coding_seq)
            aa_seq = coding_seq.translate()
            stealth_seq = _hide_motif_helper(aa_seq)
            if s_loc.strand == -1:
                stealth_seq = stealth_seq.reverse_complement()
            modified_blocks.append([s_loc, stealth_seq])

            self.loading_bar.update((loading_progress * len(s_loc)))

        return modified_blocks

    def run(self, outfile: str):
        mutable_regions = self.plasmid_input.regions()
        blocks = self.hide_motif(mutable_regions)
        plasmid_sequence = self.plasmid_input.getSeq()
        new_seq = Seq.Seq("")
        cut_start = 0

        self.loading_bar.set_description("Assembling Stealth Plasmid")

        for s_loc, seq in blocks:
            new_seq += plasmid_sequence[cut_start : s_loc.start]
            new_seq += seq
            cut_start = s_loc.end

            self.loading_bar.update((7 / len(blocks)))

        new_seq += plasmid_sequence[cut_start:]

        stealth_record = self.plasmid_input.getGenBank()
        stealth_record.seq = new_seq  # Overwrite old sequence with Stealth sequence

        with self._openOutfile(outfile) as fd:
            SeqIO.write(stealth_record, fd, "genbank")

        self.loading_bar.close()

    def writeOut(self):
        total_motifs = self.total_motifs
        remaining_motifs = self.unremoved_motifs
        removed_motifs = total_motifs - remaining_motifs

        plasmid_name = self.plasmid_input.getGenBank().name
        
        total_cds = self.plasmid_input.regionCount()
        mutable = self.plasmid_input.mutableCount()

        print(
            f"""{plasmid_name}: {mutable} NT were mutable from {total_cds} NT of CDS coverge\nRemoved {removed_motifs} of {total_motifs} motifs from {plasmid_name}. {remaining_motifs} motifs remain.""",
            file=sys.stderr,
        )


import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Stealth optimizes an input plasmid sequence to a given organism's genome",
        usage=f"""pstealth --genome (-g) <genome infile> --plasmid (-p) <plasmid infile> --outfile -o [outfile | default: stdout] [-zPMmrs]
            Optional Args:
                --zScore (-z) -> [zscore cutoff value | default: -4]
                --pseudo (-P) [pseudo-count value | default: 0]
                --max (-M) [maximum motif size | default: 8]
                --min (-m) [minimum motif size | default: 1]
                --palindrome (-r) [Remove RC palindromes only | default: off]
                --silent (-s) [Hide report message | default: show]
                --keep (-k) [Adds annotations to consider in Mutable Regions | default = ('ORF','gene','CDS')]
                --ignore (-i) [Adds annotations to ignore when defining Mutable Regions | default = ('source')]""",
    )
    parser.add_argument(
        "-g",
        "--genome",
        type=str,
        action="store",
        help="input genome file in GenBank or FastA format",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outfile",
        default=None,
        type=str,
        action="store",
        help="output filename, default = STDOUT",
    )
    parser.add_argument(
        "-p",
        "--plasmid",
        type=str,
        action="store",
        help="input plasmid file in GenBank format",
        required=True,
    )
    parser.add_argument(
        "-z",
        "--zScore",
        type=float,
        default=-4,
        action="store",
        help="set Z-score cutoff for underrepresented motifs, default = -4",
    )
    parser.add_argument(
        "-P",
        "--pseudo",
        type=int,
        default=0,
        action="store",
        help="set pseudo-count value, default = 0",
    )
    parser.add_argument(
        "-M",
        "--max",
        type= int,
        default=8,
        choices=range(2, 10),
        help="Set maximum motif size, default = 8",
    )
    parser.add_argument(
        "-m",
        "--min",
        type=int,
        default=1,
        choices=range(1, 9),
        help="Set minimum motif size, default = 1",
    )
    parser.add_argument(
        "-r",
        "--palindrome",
        default=False,
        action="store_true",
        help="Hide only RC palindromes, default = off",
    )
    parser.add_argument(
        "-s",
        "--silent",
        default=False,
        action="store_true",
        help="Disable final modification report, default = show",
    )
    
    parser.add_argument(
        "-k",
        "--keep",
        type= str,
        default= None,
        action="store",
        nargs = '+',
        help="Add plasmid annotations to what are considered mutable regions- case sensitive, Default = {'ORF','gene','CDS'}",
    )
    
    parser.add_argument(
        "-i",
        "--ignore",
        type= str,
        default= None,
        action="store",
        nargs = '+',
        help="Add plasmid annotations to ignore when defining mutable regions- case sensitive, Default = {'source'}",
    )

    args = parser.parse_args()

    genome_infile = args.genome
    plasmid_infile = args.plasmid
    outfile = args.outfile
    z_score = args.zScore
    pseudo = args.pseudo
    kMax = args.max
    kMin = args.min
    palindrome = args.palindrome
    keep_annotation = args.keep if args.keep is not None else []
    ignore_annotation = args.ignore if args.ignore is not None else []
    
    dummy = pipeline(
        genome_infile, plasmid_infile, outfile, z_score, pseudo, kMax, kMin, palindrome, keep_annotation, ignore_annotation
    )

    if not args.silent:
        dummy.writeOut()

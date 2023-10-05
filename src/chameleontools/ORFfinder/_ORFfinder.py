
class ORFfinder:
    """ORF finder class that scans through a genome and reports ORFs from all frames
    """

    def __init__(
        self,
        seq: str,
        longestGene=False,
        minGene=100,
        starts: set = {"ATG"},
        stops: set = {"TAA", "TAG", "TGA"},
    ):
        """Finds all potential ORFs in a given sequence

        Args:
            seq (str): genomic sequence
            longestGene (bool, optional): Toggles finding largest gene per ORF. Defaults to False.
            minGene (int, optional): Gene length report cutoff. Defaults to 100.
            starts (set, optional): Define start codons. Defaults to {"ATG"}.
            stops (set, optional): Define stop codons. Defaults to {"TAA", "TAG", "TGA"}.
        """
        assert minGene >= 0
        self.geneCandidates = []
        format = seq.strip().upper().replace(" ", "")
        flipTable = format[::-1].maketrans(
            "ATCG", "TAGC"
        )  # Builds reverse compliment strand

        self.fiveThree = format
        self.threeFive = format[::-1].translate(flipTable)
        self.validStart = starts
        self.validStops = stops
        self.minGene = minGene
        self.longestGene = longestGene

        self._analysis()

    def _geneFinder(self, seq, rev, longGene):
        """Helper function, collects starts and identifies stops
        """
        start = set()
        for frame in range(0, 3):
            start.clear()
            for i in range(frame, len(seq), 3):
                codon = seq[i : i + 3]

                if codon in self.validStart:
                    start.add(i)

                if codon in self.validStops and len(start) > 0:
                    _frame = -1 - frame if rev else frame + 1
                    self._addGene(start, i, _frame, longGene)
                    start = set()

            if len(start) > 0:
                _frame = -1 - frame if rev else frame + 1
                self._addGene(start, len(seq) - 3, _frame, longGene)

    def _addGene(self, start, i, frame, longGene):
        """Helper Function, does all gene-candidate calculations 
        """
        for s in sorted(start):
            length = i - s + 3
            if frame < 0:  # Handles reverse strand calcs
                startPos = len(self.fiveThree) - i - 2
                stopPos = len(self.fiveThree) - s
            else:  # Handles top strand
                startPos = s + 1
                stopPos = i + 3
            if length >= self.minGene:
                self.geneCandidates.append([startPos, stopPos, length, frame])
            if longGene:  # LongestGene flag
                return

    def _analysis(self):
        """Runs ORF finder on input sequence
        """
        self._geneFinder(self.fiveThree, False, self.longestGene)
        self._geneFinder(self.threeFive, True, self.longestGene)

        # sorts geneCandidates list, sorting by gene length
        self.geneCandidates = [
            x
            for x in sorted(
                self.geneCandidates, key=lambda entry: entry[2], reverse=True
            )
            if x[2] % 3 == 0 ]

    def get_genes(self) -> list[tuple[int, int, int, int]]:
        """Fetches all gene candidates
        
        Returns:
            list[tuple[int, int, int, int]]: List of gene candidates
        """
        return self.geneCandidates


__all__ = ["ORFfinder"]

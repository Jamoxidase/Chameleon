import chameleontools.FastAreader as FastAreader
import chameleontools.SeqParser as SeqParser
import Bio.SeqRecord
import random


class InvalidAA(KeyError):
    """
    Inherits from KeyError
    """
    def __init__(self, message) -> None:
        super().__init__(message)

class CDSanalyzer:
    """Counts up and calculates frequency of codons given a CDS

    Returns:
        CDSanalyzer: Contains counts and frequencies of codons from given CDS collection
    """
    dnaCodonTable = {
        # DNA codon table
        # T
        "TTT": "F",
        "TCT": "S",
        "TAT": "Y",
        "TGT": "C",  # TxT
        "TTC": "F",
        "TCC": "S",
        "TAC": "Y",
        "TGC": "C",  # TxC
        "TTA": "L",
        "TCA": "S",
        "TAA": "*",
        "TGA": "*",  # TxA
        "TTG": "L",
        "TCG": "S",
        "TAG": "*",
        "TGG": "W",  # TxG
        # C
        "CTT": "L",
        "CCT": "P",
        "CAT": "H",
        "CGT": "R",  # CxT
        "CTC": "L",
        "CCC": "P",
        "CAC": "H",
        "CGC": "R",  # CxC
        "CTA": "L",
        "CCA": "P",
        "CAA": "Q",
        "CGA": "R",  # CxA
        "CTG": "L",
        "CCG": "P",
        "CAG": "Q",
        "CGG": "R",  # CxG
        # A
        "ATT": "I",
        "ACT": "T",
        "AAT": "N",
        "AGT": "S",  # AxT
        "ATC": "I",
        "ACC": "T",
        "AAC": "N",
        "AGC": "S",  # AxC
        "ATA": "I",
        "ACA": "T",
        "AAA": "K",
        "AGA": "R",  # AxA
        "ATG": "M",
        "ACG": "T",
        "AAG": "K",
        "AGG": "R",  # AxG
        # G
        "GTT": "V",
        "GCT": "A",
        "GAT": "D",
        "GGT": "G",  # GxT
        "GTC": "V",
        "GCC": "A",
        "GAC": "D",
        "GGC": "G",  # GxC
        "GTA": "V",
        "GCA": "A",
        "GAA": "E",
        "GGA": "G",  # GxA
        "GTG": "V",
        "GCG": "A",
        "GAG": "E",
        "GGG": "G",  # GxG
    }

    validAA = set(dnaCodonTable.values())

    def __init__(
        self, input: str | SeqParser.StealthGenome | list[Bio.SeqRecord.SeqRecord] = None
    ) -> None:
        """Creates dictionaries to store and calculate codon counts and frequencies

        Args:
            input (str | SeqParser.StealthGenome | list[Bio.SeqRecord.SeqRecord]): FastA file or an iterable sample of CDS 
        """
        self.codonUsage = {
            "F": {"TTT": 0, "TTC": 0},
            "L": {"TTA": 0, "TTG": 0, "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0},
            "S": {"TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0, "AGT": 0, "AGC": 0},
            "Y": {"TAT": 0, "TAC": 0},
            "C": {"TGT": 0, "TGC": 0},
            "W": {"TGG": 0},
            "P": {"CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0},
            "H": {"CAT": 0, "CAC": 0},
            "Q": {"CAA": 0, "CAG": 0},
            "R": {"CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0, "AGA": 0, "AGG": 0},
            "I": {"ATT": 0, "ATC": 0, "ATA": 0},
            "M": {"ATG": 0},
            "T": {"ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0},
            "N": {"AAT": 0, "AAC": 0},
            "K": {"AAA": 0, "AAG": 0},
            "V": {"GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0},
            "A": {"GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0},
            "D": {"GAT": 0, "GAC": 0},
            "E": {"GAA": 0, "GAG": 0},
            "G": {"GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0},
            "*": {"TAA": 0, "TAG": 0, "TGA": 0},
        }
        self.count = 0
        self.codonFreq = {k: {} for k in self.codonUsage.keys()}
        self._read_file(input) if type(input) == str else self._read_StealthGenome(
            input
        )
        self._calculate()

    def _read_file(self, input: str) -> None:
        """Helper method, reads in FastA file

        Args:
            input (str): FastA file
        """
        freader = FastAreader.FastAreader(input)
        for _, cds in freader.readFasta():
            cds = cds.seq
            for i in range(len(cds), 3):
                codon = input[i : i + 3]
                codon_to_aa = self.dnaCodonTable[codon]
                self.codonUsage[codon_to_aa][codon] += 1

    def _read_StealthGenome(
        self, input: SeqParser.StealthGenome | list[Bio.SeqRecord.SeqRecord]
    ) -> None:
        """Helper method, reads CDS from StealthGenome object or iterable of SeqRecord

        Args:
            input (SeqParser.StealthGenome | list[Bio.SeqRecord.SeqRecord]): Iterable SeqRecord of CDS
        """
        input = input.cds_sequence if type(input) == SeqParser.StealthGenome else input
        for cds in input:
            cds = str(cds.seq)
            for i in range(0, len(cds), 3):
                codon = cds[i : i + 3]
                codon_to_aa = self.dnaCodonTable[codon]
                self.codonUsage[codon_to_aa][codon] += 1

    def _calculate(self):
        """Helper method, calculates codon frequency from given codon counts
        """
        for aa, codon_dict in self.codonUsage.items():
            aa_sum = sum(codon_dict.values())
            self.count = +aa_sum
            for codon in codon_dict.keys():
                if aa_sum == 0:
                    self.codonFreq[aa][codon] = 1 / len(codon_dict)
                else:
                    self.codonFreq[aa][codon] = self.codonUsage[aa][codon] / aa_sum

    def _format_aa_arg(self, aa: str):
        """Helper method, type and input checking AA inputs

        Args:
            aa (str): Single Letter AA code

        Raises:
            TypeError: Non-string type is input
            InvalidAA: Invalid AA code is input

        Returns:
            str: formatted single letter AA
        """
        if type(aa) != str:
            raise TypeError(f"Expected {type(str())} got {type(aa)}")
        aa = aa.upper()
        if len(aa) != 1 or aa not in self.validAA:
            raise InvalidAA(f"Invalid single-letter Amino Acid: {aa}")
        return aa

    def addCDS(self, cds: str) -> None:
        """Adds input CDS to existing dictionary

        Args:
            cds (str): CDS sequence, in frame
        """
        for i in range(len(cds), 3):
            codon = cds[i : i + 3]
            codon_to_aa = self.dnaCodonTable[codon]
            self.codonUsage[codon_to_aa][codon] += 1
        self._calculate()

    def len(self) -> int:
        """Length

        Returns:
            int: NT count of all analyzes CDS 
        """
        return self.count

    def getFrequency(
        self, aa: str = None
    ) -> dict[str, dict[str, float]] | dict[str, float]:
        """Get the codon frequency table or codon frequency by AA

        Args:
            aa (str, optional): Single letter AA. Defaults to None.

        Returns:
            dict[str, dict[str, float]] | dict[str, float]: Entire codon frequency table or codon frequency by AA
        """
        if aa == None:
            return self.codonFreq
        aa = self._format_aa_arg(aa)
        return self.codonFreq[aa]

    def getUsage(self, aa: str = None) -> dict[str, dict[str, int]] | dict[str, int]:
        """Get the codon usage table or codon usage by AA

        Args:
            aa (str, optional): Single letter AA. Defaults to None.

        Returns:
            dict[str, dict[str, float]] | dict[str, float]: Entire codon usage table or codon usage by AA
        """
        if aa == None:
            return self.codonUsage
        aa = self._format_aa_arg(aa)
        return self.codonUsage[aa]

class CodonOptimizer:
    """
    This is a basic CodonOptimizer that takes codon usage data and regenerates sequences
    roughly fit to levels of tRNA availibility in an organism.
    """

    def __init__(self, frequencies: dict):
        """Initializes codon usage data

        Args:
            frequencies (dict): _description_
        """
        self.synFreqDic = frequencies

    def _weightedChoice(self, aa: str):
        """Generates a codon by weighted random given an AA

        Args:
            aa (str): Single letter AA 

        Returns:
            str: Weighted random codon from given AA
        """
        codons = list(self.synFreqDic[aa].keys())
        frequencies = self.synFreqDic[aa].values()
        chosenCodon = random.choices(codons, weights=frequencies, k=1)[0]
        return chosenCodon

    def assembleSeed(self, aaSequence: str):
        """Reverse translates a AA sequence based on codon usage statistics 

        Args:
            aaSequence (str): AA sequence to be reverse translated

        Returns:
            str: DNA coding sequence
        """
        
        dnaList = [self._weightedChoice(aa) for aa in aaSequence]
        dnaSeed = "".join(dnaList)
        return dnaSeed

    def getFrequency(self) -> dict[str, dict[str, float]]:
        """Gets the codon usage of object

        Returns:
            dict[str, dict[str, float]]: Codon frequency grouped by AA
        """
        return self.synFreqDic


__all__ = ["CDSanalyzer", "CodonOptimizer", "InvalidAA"]

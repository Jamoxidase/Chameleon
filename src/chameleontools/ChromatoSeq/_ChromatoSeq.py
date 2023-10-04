import argparse
import random
from Bio.Seq import Seq


class codonOptimizer:
    """
    This is a basic codonOptimizer that takes codon useage data and a sequence
    of aminoAcids as input and recapitulates the overall frequency of the
    given input in the reverse-transcription of the AA seq. This should
    roughly fit to levels of tRNA availibility.
    """

    def __init__(self, synonomous, frequencies):
        self.synFreqDic = {}
        for aa, codons in synonomous.items():
            codFreqList = [(codon, frequencies[codon]) for codon in codons]
            self.synFreqDic[aa] = codFreqList

    def weightedChoice(self, aa):
        codons, frequencies = zip(*self.synFreqDic[aa])
        chosenCodon = random.choices(codons, weights=frequencies, k=1)
        return chosenCodon[0]

    def assembleSeed(self, aaSequence):
        dnaList = [self.weightedChoice(aa) for aa in aaSequence]
        dnaSeed = "".join(dnaList)
        return dnaSeed


class patternConstrainer:
    """
    Initialized via an instance of putativly "optimized" sequence,
    patternConstrainer finds regions of the seed that contain motifs of
    interest (via a 5'->3' sliding window approach) and passes associated
    indicies (start pos. based on specific window size & fixed adapter
    sequences at the begainning and end of each window) to some outside
    func.
    In this case, the outside func returns a 'more ideal' (synonomous variate)
    which patternConstrainer catches and overwrites into the working sequence.

    edge case handling:
    """

    def __init__(self, motifs, codonOptions, codonPref):
        self.motifs = motifs  # A list of motifs to avoid, likely to trigger RMS
        self.codonOptions = codonOptions  # The dictionary of codons mapping to their corresponding amino acids.
        self.codonPref = codonPref

    def motifInSeq(self, sequence):
        """
        remove? if we use a motif checker that checks for a motif at each position, we can add variation to the shuffles actual footprint size
        """
        # Iterate over each motif
        for motif in self.motifs:
            # If the motif is found in the sequence
            if motif in sequence:
                return True  # Return True immediately. This indicates that the sequence contains a motif that we want to avoid.
        return False  # If none of the motifs were found, return False. This indicates that the sequence is safe to use as is.

    def optimizeCodon(self, sequence, start):
        length = len(sequence)
        workingDnaSequence = sequence[start : start + 12]

        bioSeq = Seq(workingDnaSequence)
        workingAaSeq = bioSeq.translate()

        if start < 9:
            headStart = start
        elif 9 <= start:
            headStart = start - 9

        if length <= start + 21:
            tailEnd = length
            tailStart = length
        else:
            tailEnd = start + 21
            tailStart = start + 12

        headAdapt = sequence[
            headStart:start
        ]  # change to handle begainning case (no -6 index)
        tailAdapt = sequence[tailStart:tailEnd]

        myPermuter = Permutations(self.motifs, self.codonOptions, self.codonPref)
        cleanedRegion = myPermuter.rank(
            workingDnaSequence, workingAaSeq, headAdapt, tailAdapt
        )

        cleanCutOff = len(cleanedRegion) - len(tailAdapt)
        newSequence = (
            sequence[:start]
            + cleanedRegion[len(headAdapt) : cleanCutOff]
            + sequence[start + 12 :]
        )

        return newSequence

    def optimizeSequence(self, sequence):
        """
        understand behaviour
        """
        optimizedSequence = sequence
        for i in range(9, len(sequence) - 9, 3):
            if self.motifInSeq(optimizedSequence[i : i + 12]):
                newSequence = self.optimizeCodon(optimizedSequence, i)
                optimizedSequence = newSequence
        return optimizedSequence


class Permutations:
    """
    Generates all possible codon permutations of given window (codon subsequence
    of sequence) & ranks permutations based on rank items as they are generated.
    Each round of generation ends with the appension of scores and permutations
    added to net list as a tuple. These tuples are then sorted via score descrimination-
    only canditates with the best score for first item contenue, same with second item
    and then the 3rd item varies. This can be adjusted but allows for the prioritization
    of sorting methodologies for each specific score items.
    prior implementations had score bleeding issues, this is currently an open end
    """

    def __init__(self, motifs, codonOptions, frequencies):
        self.Motifs = motifs
        self.frequencies = frequencies
        self.codonOptions = codonOptions

    def permute(self, amino_acid_sequence, current_codons=""):
        if not amino_acid_sequence:
            yield current_codons
            return

        current_aa = amino_acid_sequence[0]
        remaining_aas = amino_acid_sequence[1:]

        for codon in self.codonOptions[current_aa]:
            yield from self.permute(remaining_aas, current_codons + codon)

    def rank(self, dnaSeq, aaSeq, headAdapt, tailAdapt):
        """ """
        ranks = []
        permutation_generator = self.permute(aaSeq)

        while True:
            try:
                putative = next(permutation_generator)
                permutation = headAdapt + putative + tailAdapt

                severity, codScore = self.analyze_dna_motifs(permutation)
                adherance = self.gaugeAdherance(
                    dnaSeq, putative
                )  # returns a value between 0 and 1, 0 being less adherant

                ranks.append((severity, adherance, codScore, permutation))
            except StopIteration:
                break  # Stop when all permutations have been generated

        min_first_item = min(ranks, key=lambda x: x[0])[0]
        filter1 = sorted(
            [t for t in ranks if t[0] == min_first_item],
            key=lambda x: x[1],
            reverse=True,
        )  # Sort the remaining tuples based on the second item in descending order

        max_second_item = filter1[0][1]
        filter2 = [
            t for t in filter1 if t[1] == max_second_item
        ]  # Extract tuples with the highest second item values

        max_third_item = max(filter2, key=lambda x: x[2])[2]
        filter3a = [
            t for t in filter2 if 0.95 * max_third_item <= t[2] <= 1.05 * max_third_item
        ]  # Extract tuples with a third item value within 5% of the maximal value
        filter3b = sorted(filter3a, key=lambda x: x[2], reverse=True)

        # add another layer, for now first item
        chosenTup = filter3b[0]
        chosenOne = chosenTup[3]
        return chosenOne

    def analyze_dna_motifs(self, sequence):
        occurrences = {}  # Dictionary to store occurrences of short motifs
        total_occurrences = 0  # Total count of occurrences in the longer sequence

        for short_motif in self.Motifs:
            motif_positions = []  # List to store positions of motif occurrences
            motif_count = 0  # Count of occurrences for each motif
            start_position = 0

            while start_position < len(sequence):
                position = sequence.find(short_motif, start_position)

                if position == -1:
                    break  # No more occurrences found
                motif_positions.append(position)
                motif_count += 1
                total_occurrences += 1
                start_position = position + 1  # Move to the next position

            occurrences[short_motif] = {
                "positions": motif_positions,
                "count": motif_count,
            }

        score = 0.0
        codLen = len(sequence) / 3
        for i in range(0, len(sequence), 3):
            codon = sequence[i : i + 3]
            if codon in self.frequencies:
                score += self.frequencies[codon]
                Score = score / codLen
        return total_occurrences, Score

    def gaugeAdherance(self, refSeq, perm):
        codonRef = []
        codonPerm = []
        adherance = 0
        """
        for i in range(0, len(refSeq), 3): #both refSeq & perm should be of the same length
            codon = refSeq[ i : i + 3 ]
            codonRef.append(codon)
        for i in range(0, len(perm), 3):
            codon = perm[ i : i + 3 ]
            codonPerm.append(codon)
        """
        for i in range(0, len(refSeq), 3):
            refCodon = refSeq[i : i + 3]
            codonRef.append(refCodon)

            perCodon = perm[i : i + 3]
            codonPerm.append(perCodon)

        adherance = 0
        codLen = len(refSeq) / 3
        for tup in zip(codonRef, codonPerm):
            cod1 = tup[0]
            cod2 = tup[1]
            if cod1 == cod2:
                adherance += 1 / codLen
        return adherance  # confers percentage of conserved codons in perm


class MotifChecker:
    """
    helper class, used to fetch final stats in main
    """

    def __init__(self, motifs):
        self.motifs = motifs

    def checkMotifs(self, sequence):
        totalOccurrences = 0  # Total count of occurrences in the longer sequence

        for motif in self.motifs:
            startPosition = 0
            while startPosition < len(sequence):
                position = sequence.find(motif, startPosition)

                if position == -1:
                    break  # No more occurrences found

                totalOccurrences += 1
                startPosition = position + 1  # Move to the next position

        return totalOccurrences


def main():
    # codonOptions is a constant
    codonOptions = {
        "F": ["TTT", "TTC"],
        "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "Y": ["TAT", "TAC"],
        "C": ["TGT", "TGC"],
        "W": ["TGG"],
        "P": ["CCT", "CCC", "CCA", "CCG"],
        "H": ["CAT", "CAC"],
        "Q": ["CAA", "CAG"],
        "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "I": ["ATT", "ATC", "ATA"],
        "M": ["ATG"],
        "T": ["ACT", "ACC", "ACA", "ACG"],
        "N": ["AAT", "AAC"],
        "K": ["AAA", "AAG"],
        "V": ["GTT", "GTC", "GTA", "GTG"],
        "A": ["GCT", "GCC", "GCA", "GCG"],
        "D": ["GAT", "GAC"],
        "E": ["GAA", "GAG"],
        "G": ["GGT", "GGC", "GGA", "GGG"],
        "*": ["TAA", "TAG", "TGA"],
    }

    # Create a parser with appropriate descriptions for the arguments
    parser = argparse.ArgumentParser(description="DNA Processing")
    parser.add_argument("amino_sequence", help="Amino sequence as a string")
    parser.add_argument("frequencies", help="Frequencies as a dictionary string")
    parser.add_argument("motifs", help="Motifs as a list string")

    args = parser.parse_args()

    # Parse the input arguments
    aaSequence = args.amino_sequence

    # Convert the input strings to Python dictionaries and lists
    frequencies = eval(args.frequencies)  # Use eval to convert to a dictionary
    motifs = eval(args.motifs)  # Use eval to convert to a list

    """
    frequencies = {'TAA': 0.526, 'TAG': 0.264, 'TGA': 0.21, 'GCA': 0.182, 'GCC': 0.316, 
                    'GCG': 0.168, 'GCT': 0.334, 'TGC': 0.259, 'TGT': 0.741, 'GAC': 0.233,
                    'GAT': 0.767, 'GAA': 0.741, 'GAG': 0.259, 'TTC': 0.242, 'TTT': 0.758,
                    'GGA': 0.292, 'GGC': 0.184, 'GGG': 0.182, 'GGT': 0.343, 'CAC': 0.351,
                    'CAT': 0.649, 'ATA': 0.128, 'ATC': 0.324, 'ATT': 0.548, 'AAA': 0.791,
                    'AAG': 0.209, 'CTA': 0.126, 'CTC': 0.129, 'CTG': 0.117, 'CTT': 0.099,
                    'TTA': 0.389, 'TTG': 0.14, 'ATG': 1.0, 'AAC': 0.277, 'AAT': 0.723,
                    'CCA': 0.162, 'CCC': 0.411, 'CCG': 0.172, 'CCT': 0.256, 'CAA': 0.695,
                    'CAG': 0.305, 'AGA': 0.19, 'AGG': 0.067, 'CGA': 0.153, 'CGC': 0.218, 
                    'CGG': 0.139, 'CGT': 0.234, 'AGC': 0.142, 'AGT': 0.266, 'TCA': 0.107,
                    'TCC': 0.173, 'TCG': 0.103, 'TCT': 0.211, 'ACA': 0.18, 'ACC': 0.36, 
                    'ACG': 0.139, 'ACT': 0.32, 'GTA': 0.181, 'GTC': 0.227, 'GTG': 0.264, 
                    'GTT': 0.328, 'TGG': 1.0, 'TAC': 0.283, 'TAT': 0.717}
    motifs = ['CACCTGC', 'GTGGACG', 'CCGG', 'CGCG', 'GGCC', 'TGCA', 'CCAGG', 'CCTGG', 'GGACC', 'GGCCC', 'GGGCC', 'GGTCC', 'TCTGA', 'AAATTT', 'AATATT', 'ACATGT', 'ACGCGT', 'ACTAGT', 'AGCGCT', 'AGTACT', 'ATCGAT', 'ATGCAT', 'CAATTG', 'CAGCTG', 'CATATG', 'CCGCGG', 'CCTAGG', 'CGCGCG', 'CGTACG', 'CTCGAG', 'CTGCAG', 'GAATTC', 'GACGTC', 'GAGCTC', 'GATATC', 'GCATGC', 'GCGCGC', 'GCTAGC', 'GGCGCC', 'GGTACC', 'GTATAC', 'GTCGAC', 'GTGCAC', 'TAATTA', 'TACGTA', 'TCCGGA', 'TGCGCA', 'TGTACA', 'TTCGAA', 'TTTAAA', 'GGCC', 'ACTAGT', 'AGTACT', 'ATCGAT', 'ATGCAT', 'CTCGAG', 'CTGCAG', 'GAATTC', 'GCTAGC', 'GGTACC', 'GTCGAC', 'TAATTA', 'TGCGCA', 'TTCGAA', 'CCTCAGG', 'CCTGAGG', 'GCTAAGC', 'GCTTAGC', 'GGTAACC', 'GGTCACC', 'GGTGACC', 'GGTTACC', 'AAAATTTT', 'AAACGTTT', 'AAAGCTTT', 'AAATATTT', 'AACATGTT', 'AACCGGTT', 'AACGCGTT', 'AACTAGTT', 'AATATATT', 'AATCGATT', 'AATGCATT', 'AATTAATT', 'AGATATCT', 'AGTTAACT', 'ATAATTAT', 'ATTATAAT', 'CAGATCTG', 'CAGTACTG', 'CATATATG', 'CATTAATG', 'CCAATTGG', 'CCACGTGG', 'CCAGCTGG', 'CCCATGGG', 'CCCCGGGG', 'CCCGCGGG', 'CGAGCTCG', 'CGATATCG', 'CGCGCGCG', 'CGCTAGCG', 'CTAATTAG', 'CTAGCTAG', 'CTATATAG', 'CTGATCAG', 'CTGGCCAG', 'CTGTACAG', 'CTTATAAG', 'CTTGCAAG', 'CTTTAAAG', 'GAGATCTC', 'GAGTACTC', 'GATATATC', 'GATTAATC', 'GTGATCAC', 'GTTATAAC', 'TACCGGTA', 'TACTAGTA', 'TATCGATA', 'TATTAATA', 'TGATATCA', 'TGTTAACA']
    aaSequence = 'MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK'
    """

    seedDNA = codonOptimizer(codonOptions, frequencies)
    avoidPatterns = patternConstrainer(motifs, codonOptions, frequencies)
    evalMortality = MotifChecker(motifs)

    roundsData = []
    i = 0
    while i < 5:
        seed = seedDNA.assembleSeed(aaSequence)
        optimization = avoidPatterns.optimizeSequence(seed)

        mortality = evalMortality.checkMotifs(optimization)
        roundsData.append((mortality, optimization))
        i += 1
    lowestMortality = min(roundsData, key=lambda x: x[0])
    best = lowestMortality[0]
    bestSeq = lowestMortality[1]

    print(best, "MOTIFS REMAIN: \n", bestSeq)


if __name__ == "__main__":
    main()

from Bio.Seq import Seq

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
    """

    def __init__(self, motifs, codonPref):
        self.motifs = motifs  # A list of motifs to avoid, likely to trigger RMS
        self.codonOptions = { # The dictionary of codons mapping to their corresponding amino acids.
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
        self.codonPref = codonPref

    def motifInSeq(self, sequence):
        """
        Given a sequence and initialized with a list of motifs, evaluates presence or absence of motifs in the sequence.
        In practice this method is used to support the sliding window approach to pattern constraint requirements.
        """
        for motif in self.motifs:
            if motif in sequence:
                return True  # Return True immediately. This indicates that the sequence contains a motif that we want to avoid.
        return False  # If none of the motifs were found, return False. This indicates that the sequence is safe to use as is.

    def optimizeCodon(self, sequence, start):
        '''
        Organizes optimization of the codon sequence at a given start position in a DNA sequence.

        This method takes a DNA sequence and a start position as input. It extracts a 12-base pair
        (12-bp) subsequence starting from the given position and translates it into an amino acid sequence.
        The method then calculates head and tail subsequences based on the start position and a fixed size
        window. The head and tail subsequences are used to integrate the current window of the sliding
        window approach with the rest of the sequence.

        The method employs a Permutations object to rank the 12-bp DNA sequence based on motifs and codon preferences.
        The resulting cleaned region is then inserted back into the original DNA sequence, replacing the original
        12-bp sequence. The final sequence with the optimized codon is returned.

        We chose a 12-bp/ 4 codon window as the longest motif we are screening is 8-bp long, and to resolve that the 
        window must have at least a 4 codon scope. The reason that the window is not larger is to retain stochastic
        reconstruction as much as possible, and to reduce the footprint of the pattern constraint optimization.
        '''
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

        myPermuter = Permutations(self.motifs, self.codonPref)
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

    def optimizeSequence(self, sequence) -> str:
        """
        This method iterates through the input DNA sequence, looking for motifs within 12-base pair (12-bp) subsequences.
        When a motif is found, the method calls the 'optimizeCodon' method to optimize the codon at the motif's start position.
        The optimized DNA sequence is updated, and the process continues until all relevant motifs have been optimized.
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
    """
    
    def __init__(self, motifs, frequencies):
        self.Motifs = motifs
        self.frequencies = {}
        for _,codon_dict in frequencies.items():
            self.frequencies.update(codon_dict)
        self.codonOptions = { # The dictionary of codons mapping to their corresponding amino acids.
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
        

    def permute(self, amino_acid_sequence, current_codons=""):
        '''
        Recusive generator to produce all possible codon permutations for a given amino acid sequence. 
        '''
        if not amino_acid_sequence:
            yield current_codons
            return

        current_aa = amino_acid_sequence[0]
        remaining_aas = amino_acid_sequence[1:]

        for codon in self.codonOptions[current_aa]:
            yield from self.permute(remaining_aas, current_codons + codon)

    def rank(self, dnaSeq, aaSeq, headAdapt, tailAdapt):
        """ 
        Organizes evalutation of each permutation (putative solution) for motifs present, stochasticity, & preferential 
        codons to select optimal solutions. This method uses the head & tail adapter sequences to consider cases where 
        motifs are split by the sliding window.
        """
        ranks = []
        permutation_generator = self.permute(aaSeq)

        for putative in permutation_generator:
            permutation = headAdapt + putative + tailAdapt

            severity, codScore = self.analyze_dna_motifs(permutation)
            adherance = self.gaugeAdherance(
                dnaSeq, putative
            )  # returns a value between 0 and 1, 0 being less adherant

            ranks.append((severity, adherance, codScore, permutation))

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
        '''
        Supports rank(), this method returns score values for motifs present & preferential 
        codon useage to elute optimal solutions.
        '''
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
        '''
        Supports rank by providing score value reflective of stochasticity (origional input is stochastic).
        '''
        codonRef = []
        codonPerm = []
        adherance = 0
        
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
    helper class, used to fetch final stats on motifs present in main
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

__all__ = ["patternConstrainer","Permutations","MotifChecker"]

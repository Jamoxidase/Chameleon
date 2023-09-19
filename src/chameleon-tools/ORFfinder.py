"""
TABI - UCSC iGEM 2023

Modified from ORFfinder Class written for Winter 2022 BME160

BME160 Professor: David L. Bernick

Author: Tyler Gaw
"""

# ORF finder used to parse FastA format genome files.
# Optional parameters:
#     longestGene -> only reports longest genes in an ORF
#     minGene -> only reports genes of this size or larger
#     starts -> define set of valid start codons
#     stop -> define set fo valid stop codons

class ORFfinder:
    """
    ORF finder class that scans through a genome and reports ORFs from all frames
    
    Data is reported in a list of lists in 
    
    format: list[[startPos,stopPos,length,frame]]
    """

    def __init__(self,seq: str,longestGene=False,minGene = 100,starts: set = {"ATG"},stops: set = {"TAA","TAG","TGA"}):
        """
        Usage:
        gene_finder = ORFfinder(seq,**kwargs)
        
        gene_candidates = gene_finder.get_genes()
        """
        assert minGene >= 0
        self.geneCandidates = []
        format = seq.strip().upper().replace(" ",'') 
        flipTable = format[::-1].maketrans("ATCG", "TAGC") # Builds reverse compliment strand
        
        self.fiveThree = format 
        self.threeFive = format[::-1].translate(flipTable) 
        self.validStart = starts
        self.validStops = stops
        self.minGene = minGene
        self.longestGene = longestGene
        
        self._geneFinder()

    def _geneFinder(self,seq,rev,longGene):
        """
        Helper function that collects starts and identifies stops
        """
        start = set()
        for frame in range(0,3):
            start.clear()
            for i in range(frame,len(seq),3):
                codon = seq[i:i+3]

                if codon in self.validStart:
                    start.add(i)

                if codon in self.validStops and len(start) > 0:
                    _frame = -1-frame if rev else frame+1 
                    self._addGene(start,i,_frame,longGene)
                    start = set() 
                    

            if len(start) > 0: 
                _frame = -1-frame if rev else frame+1
                self._addGene(start,len(seq)-3,_frame,longGene)


    def _addGene(self,start,i,frame,longGene):
        """
        Helper Function
        Does all gene-candidate calculations using [start] [stop index] [frame #] [longest gene flag]
        """
        for s in sorted(start):
            length = i - s + 3
            if frame < 0: # Handles reverse strand calcs
                startPos = len(self.fiveThree) - i -2
                stopPos = len(self.fiveThree) - s
            else: # Handles top strand
                startPos = s + 1
                stopPos = i + 3
            if length >= self.minGene:
                self.geneCandidates.append([startPos, stopPos, length, frame])
            if longGene: # LongestGene flag
                return


    def _geneFinder(self):
        """
        Main geneFinder function. Runs hidden _geneFinder() function to parse through both top/bottom strand
        """
        self._geneFinder(self.fiveThree,False,self.longestGene)
        self._geneFinder(self.threeFive,True,self.longestGene)
        
        #sorts geneCandidates list, sorting by gene length
        self.formatted = [x for x in sorted(self.geneCandidates, key=lambda entry: entry[2], reverse=True)] 
        return self.formatted

    def get_genes(self):
        return self.geneCandidates
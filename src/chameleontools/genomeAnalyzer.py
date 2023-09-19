#!/usr/bin/env python3
# Name: Allison Jaballas (acjaball)

##########################################################################
############################## FASTA READER ##############################
##########################################################################

import sys
class FastAreader:
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

########################################################################
############################## NUC PARAMS ##############################
########################################################################

class NucParams:
    '''
    Define objects to store counts of codons and amino acids from input sequence.
    '''
    
    dnaCodonTable = {
    # Second Base
    # T             C             A             G
    #T
    'TTT': 'Phe', 'TCT': 'Ser', 'TAT': 'Tyr', 'TGT': 'Cys',
    'TTC': 'Phe', 'TCC': 'Ser', 'TAC': 'Tyr', 'TGC': 'Cys',
    'TTA': 'Leu', 'TCA': 'Ser', 'TAA': '---', 'TGA': '---',
    'TTG': 'Leu', 'TCG': 'Ser', 'TAG': '---', 'TGG': 'Trp',
    #C 
    'CTT': 'Leu', 'CCT': 'Pro', 'CAT': 'His', 'CGT': 'Arg',
    'CTC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
    'CTA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
    'CTG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
    #A
    'ATT': 'Ile', 'ACT': 'Thr', 'AAT': 'Asn', 'AGT': 'Ser',
    'ATC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
    'ATA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
    'ATG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
    #G
    'GTT': 'Val', 'GCT': 'Ala', 'GAT': 'Asp', 'GGT': 'Gly',
    'GTC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
    'GTA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
    'GTG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
    }

    def __init__ (self, inString=''):
        ''' Set up codon composition, amino acid composition, and nucleotide composition dictionaries for future use.'''
        self.inString = inString.upper()
        
        self.codonComp = {codon: 0 for codon in self.dnaCodonTable.keys()} #set up empty (0) dict for codon counts

        self.aaComp = {
            'Phe': 0, 'Leu': 0, 'Ser': 0, 'Tyr': 0, '---': 0, 'Cys': 0, 'Trp': 0, 'Pro': 0,
            'His': 0, 'Gln': 0, 'Arg': 0, 'Ile': 0, 'Met': 0, 'Thr': 0, 'Asn': 0, 'Lys': 0,
            'Val': 0, 'Ala': 0, 'Asp': 0, 'Glu': 0, 'Gly': 0
            }
        
        self.nucComp = { 'A': 0, 'C': 0, 'G' : 0, 'T' : 0, 'U' : 0, 'N' : 0 } 

        self.addSequence(inString)
        
    def addSequence (self, inSeq):
        '''Add nucleotide count totals and decode codon information to appropriate dictionaries.'''
        inSeq = self.inString.join(inSeq.split()).upper()

        # count totals of each nucleotide from file input
        for nuc in self.nucComp.keys():
            self.nucComp[nuc] += inSeq.count(nuc)

        # to decode codon -> AA
        AA = [ inSeq[n:n+3] for n in range(0, len(inSeq), 3) ]
        
        # add to codonComp + aaComp dictionaries!
        for codon in AA:
            if codon in self.codonComp: 
                self.codonComp[codon] += 1
                aa = self.dnaCodonTable[codon]
                self.aaComp[aa] += 1
       
    def nucComposition(self):
        '''Return dictionary counts of valid nucleotides from file input.'''
        return self.nucComp
    def codonComposition(self):
        '''Return dictionary counts of codons from file input.'''
        return self.codonComp
    def nucCount(self):
        '''Return sum of all the valid nucleotides based on nucComposition.'''
        return sum(self.nucComp.values())

#!/usr/bin/env python3
# Name: Allison Jaballas (acjaball)

#############################################################################
############################## GENOME ANALYZER ##############################
#############################################################################

class genomeAnalyzer:
    '''
    GenomeAnalyzer calls the methods from the classes imported from sequenceAnalyzer to calculate the
    file input's sequence length, GC content, and relative codon usages. This program prepares the
    summaries and final display of all the data.

    Dictionaries are as follows:

    totalCodonsPerAA = {codon:#count}
    codonFractionUsed = {codon:fraction}

    '''
    def __init__ (self, filename = None):
        '''Constructor: saves attribute filename.'''
        self.calculations(filename)
        
    def calculations (self,fastafile):
        '''Calculate the file input's sequence length and relative codon usages.'''
        myReader = FastAreader(fastafile)
        myNuc = NucParams()
        for head,seq in myReader.readFasta():
            myNuc.addSequence(seq)
            
        seqLength = myNuc.nucCount()/1000000
        #print('sequence length = {:.2f} Mb'.format(seqLength))
        #print()

        codonCount = myNuc.codonComposition() #shorthand version
        
        # calculate fraction used
        # set up dictionaries
        # totalCodonsPerAA gets the total # of codons for each aa
        codonFractionUsed = {}
        totalCodonsPerAA = {aa: sum(count for codon, count in myNuc.codonComp.items() if myNuc.dnaCodonTable[codon] == aa)for aa in myNuc.aaComp}
        
        codonTot = 0

        for codon, count in myNuc.codonComp.items(): # calculate fraction
            aa = myNuc.dnaCodonTable[codon] #the aa of each codon
                
            if totalCodonsPerAA[aa] > 0:
                codonFractionUsed[codon] = count/totalCodonsPerAA[aa] # builds the dict codonFractionUsed

            codonTot += count
                                    
        for codons in sorted(codonCount.keys()): # pull together all information and prepare to print
            codCount = codonCount[codons]
            fracCount = codonFractionUsed[codons]
            aa = myNuc.dnaCodonTable[codons] if myNuc.dnaCodonTable[codons] != "---" else "End"
            aaCount = myNuc.aaComp[aa] if aa != 'End' else myNuc.aaComp['---']

            # codon freq calc: #codons/totalAA
            codonUsage = (codCount/codonTot) * 1000

            # print everything together
            print('{:5} {:5} {:5d} {:5.2f} {:5.2f}'.format(aa, codons, codCount, codonUsage, fracCount))
   
genomeAnalyzer("microcystis_blast_contigs_CDS.fasta")

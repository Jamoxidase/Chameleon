from typing import Iterator
import argparse,sys

class FastAreader:
    """
    Define objects to read FastA files.

    Args:
    fname (str): file name (optional), default is None (STDIN)

    Usage:
    thisReader = FastAreader ('testTiny.fa')
    for head, seq in thisReader.readFasta():
        print (head,seq)
    """

    def __init__(self, fname=None):
        """
        Constructor: saves attribute fname
        Args:
        fname (str): file name (optional), default is None (STDIN)
        """
        self.fname = fname

    def doOpen(self):
        """
        Handle file opens, allowing STDIN.

        Returns:
        file: either sys.stdin or file handler for the input file
        """
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        """
        Read an entire FastA record and return the sequence header/sequence

        Yields:
        tuple: header and sequence as a tuple
        """
        header = ''
        sequence = ''

        with self.doOpen() as fileH:
            header = ''
            sequence = '\n'

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = '\n'
                else:
                    sequence += line
        yield header, sequence

class ParseStealth(set):
    '''
    Takes stealth output data as a text file and processes a set of motifs.
    '''
    degenerateBases = { # Degenerate base pair dictionary
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'M': ['A', 'C'],
            'K': ['G', 'T'],
            'S': ['C', 'G'],
            'W': ['A', 'T'],
            'H': ['A', 'C', 'T'],
            'B': ['C', 'G', 'T'],
            'V': ['A', 'C', 'G'],
            'D': ['A', 'G', 'T'],
            'N': ['A', 'C', 'G', 'T']
            }

    def __init__(self, file_path: str):
        '''
        Constructor: Inherits from set\n
        args : <file_path> - File path to stealth output\n
        Output: Set containing expanded motifs
        '''
        super().__init__()
        self._readIn(file_path)
    
    def _readIn(self,file_path: str) -> set:
        '''
        helper func - from file reads all motifs into set
        '''
        file_p = open(file_path,'r') if type(file_path) == str else sys.stdin
        with file_p as fd:
            fd.readline()
            while line := fd.readline().strip().split():
                for expanded_motif in self._permute(line[0]):
                    self.add(expanded_motif)
    
    def _permute(self,string, cur_idx=0, cur_perm="") -> Iterator[str]:
        '''
        helper func - generates permutations of all degenerate sequences
        '''
        if cur_idx == len(string):
            yield cur_perm
            if cur_perm == string:
                return
        else:
            cur_let = string[cur_idx]
            if cur_let in self.degenerateBases:
                subs = self.degenerateBases[cur_let]
                for sub in subs:
                    yield from self._permute(string,cur_idx + 1,cur_perm + sub)
            else:
                yield from self._permute(string,cur_idx + 1,cur_perm + cur_let)

class PalindromeParseStealth(ParseStealth):
    
    def __init__(self, file_path: str) -> set:
        '''
        Constructor: Inherits from ParseStealth --> Set\n
        args : <file_path> - File path to stealth output\n
        Output: Set containing expanded RC palindrome motifs
        '''
        super().__init__(file_path)
        
    def _readIn(self,file_path: str) -> set: 
        '''
        helper func overload - from file reads palindromes into set
        '''
        file_p = open(file_path,'r') if type(file_path) == str else sys.stdin
        with file_p as fd:
            fd.readline()
            while line := fd.readline().strip().split():
                if line[-1] == 'Palindrome':
                    for expanded_motif in self._permute(line[0]):
                        if self._revComp(expanded_motif):
                            self.add(expanded_motif)

    def _revComp(self,seq: str) -> bool:
        '''
        helper function - boolean check on reverse compliment
        '''
        t = {'A':'T','C':'G','T':'A','G':'C'}
        return ''.join([t[i] for i in seq[::-1]]) == seq
    
def writeFile(motif: set, outfile: str, sortBool: bool):
    '''
    writes stealth output to DNAworks format file
    '''
    file_p = open(outfile,"w") if type(outfile) == str else outfile
    with file_p as fd:
        if sortBool:    
            'explicit sort by length first, then alphabetical'
            for seq in sorted(list(motif),key= lambda x : (len(x),x)):
                fd.write(seq+" [temp]\n")
        else:
            for seq in list(motif):
                fd.write(seq+" [temp]\n")

def main():
    '''
    Commandline parse, executes file write
    '''
    parser = argparse.ArgumentParser(description= "Reads in Stealth file, outputs motifs in DNAworks compatible file",usage= f"{sys.argv[0]} -i <input file | default= stdin> -o <outfile | default stdout> -[optional]\n\t-s [sorted | default= False]\n\t-p [RC paindromes only | default= False]")
    parser.add_argument("--infile",'-i',default=sys.stdin,type=str,action='store',help='input file directory')
    parser.add_argument("--outfile",'-o',default=sys.stdout,type=str,action='store',help='output file')
    parser.add_argument("--sorted",'-s',default=False,action='store_true',help='sort output')
    parser.add_argument("--palindrome",'-p',default=False,action='store_true',help='palindrome output')
    args = parser.parse_args()
    conserved = PalindromeParseStealth(args.infile) if args.palindrome else ParseStealth(args.infile)
    writeFile(conserved,args.outfile,args.sorted)

if __name__ == "__main__":
    main()








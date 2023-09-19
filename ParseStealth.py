"""
TABI - UCSC iGEM 2023

Classes to parse Stealth output and isolate restricted motifs

Author: Tyler Gaw
"""

# Takes in a file path to Stealth output and generates a set of expanded degenerate motifs
# If no file_path is specified, reads from STDIN
# Usage Example:
#     motif_set_1 = ParseStealth(file_path_1)
#     motif_set_2 = PalindromeParseStealth(file_path_2)
#     conserved_motifs = motif_set_1.intersect(motif_set_2) -> set

from typing import Iterator
import sys

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
            while line := fd.readline():
                if not line.startswith("N ="):
                    continue
                line = fd.readline().strip().split()
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
            
            while line := fd.readline():
                
                if not line.startswith("N ="): # Skip to first motif
                    continue
                
                line = fd.readline().strip().split()
                if line[-1] == 'Palindrome':
                    for expanded_motif in self._permute(line[0]):
                        if self._revComp(expanded_motif):
                            self.add(expanded_motif)

    def _revComp(self,seq: str) -> bool:
        '''
        helper function - boolean check on reverse compliment
        '''
        t = {'A':'T','C':'G','T':'A','G':'C'}
        if len(seq) % 2 == 0:
            return ''.join([t[i] for i in seq[::-1]]) == seq
        else:
            return "".join([t[i] for i in seq[:len(seq)//2:-1]]) == seq[:len(seq)//2]
    
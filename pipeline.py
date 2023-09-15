#!/usr/bin/env python3
"""
TABI - UCSC iGEM 2023

Commandline interface for parsing Stealth output

Standalone program for ParseStealth.py

Author: Tyler Gaw
"""

# Usage: ./pipeline.py -i <input file | default= stdin> -o <outfile | default stdout> -[optional]
#         -s [sorted | default= False]
#         -p [RC paindromes only | default= False]
#         -d DNAworks compatible output


import argparse,sys
from ParseStealth import ParseStealth,PalindromeParseStealth

def writeFile(motif: set, outfile: str, sortBool: bool, header: str, dnaWorks: bool):
    '''
    writes stealth output to DNAworks format file
    '''
    file_p = open(outfile,"w") if type(outfile) == str else sys.stdout
    with file_p as fd:
        fd.write(f"Run Options: {header}\n")
        if sortBool:    
            'explicit sort by length first, then alphabetical'
            for i,seq in enumerate(sorted(list(motif),key= lambda x : (len(x),x))):
                descriptor = f" re{i}" if dnaWorks else f" [temp]"
                fd.write(f"{descriptor}\t{seq}\n")
        else:
            for seq in list(motif):
                descriptor = f" re{i}" if dnaWorks else f" [temp]"
                fd.write(f"{descriptor}\t{seq}\n")
    

def main():
    '''
    Commandline parse, executes file write
    '''
    parser = argparse.ArgumentParser(description= "Reads in Stealth file, outputs motifs in formatted file",usage= f"{sys.argv[0]} -i <input file | default= stdin> -o <outfile | default stdout> -[optional]\n\t-s [sorted | default= False]\n\t-p [RC paindromes only | default= False]\n\t-d DNAworks compatible output")
    parser.add_argument("--infile",'-i',default=None,type=str,action='store',help='input file directory')
    parser.add_argument("--outfile",'-o',default=None,type=str,action='store',help='output file')
    parser.add_argument("--sorted",'-s',default=False,action='store_true',help='sort output')
    parser.add_argument("--palindrome",'-p',default=False,action='store_true',help='palindrome output')
    parser.add_argument("--dnaWorks",'-d',default=False,action='store_true',help='DNAWorks compatible output')
    args = parser.parse_args()
    conserved = PalindromeParseStealth(args.infile) if args.palindrome else ParseStealth(args.infile)
    writeFile(conserved,args.outfile,args.sorted," ".join(sys.argv),args.dnaWorks)


if __name__ == "__main__":
    main()








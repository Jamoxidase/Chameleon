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


import argparse, sys, re
from chameleontools.StealthParser import ParseStealth


def writeFile(motif: set, outfile: str, sortBool: bool, header: str, dnaWorks: bool):
    """
    writes stealth output to DNAworks format file
    """
    file_p = open(outfile, "w") if type(outfile) == str else sys.stdout
    with file_p as fd:
        fd.write(f"Run Options: {header}\n")
        if sortBool:
            "explicit sort by length first, then alphabetical"
            for i, seq in enumerate(sorted(list(motif), key=lambda x: (len(x), x))):
                descriptor = f" re{i}" if dnaWorks else f" [temp]"
                fd.write(f"{descriptor}\t{seq}\n")
        else:
            for seq in list(motif):
                descriptor = f" re{i}" if dnaWorks else f" [temp]"
                fd.write(f"{descriptor}\t{seq}\n")


def main():
    """
    Commandline parse, executes file write
    """
    script_name = re.split(r"\\|/", sys.argv[0])[-1]
    parser = argparse.ArgumentParser(
        description="Reads in Stealth file, outputs motifs in formatted file",
        usage=f"{script_name} -i <input file> -o <outfile | default stdout> -[optional]\n\t-s [sorted | default= False]\n\t-p [RC paindromes only | default= False]\n\t-d DNAworks compatible output",
    )
    parser.add_argument(
        "--infile", "-i", type=str, action="store", help="input file", required=True
    )
    parser.add_argument(
        "--outfile",
        "-o",
        default=sys.stdout,
        type=str,
        action="store",
        help="output file",
    )
    parser.add_argument(
        "--sorted", "-s", default=False, action="store_true", help="sort output"
    )
    parser.add_argument(
        "--palindrome",
        "-p",
        default=False,
        action="store_true",
        help="palindrome output",
    )
    parser.add_argument(
        "--dnaWorks",
        "-d",
        default=False,
        action="store_true",
        help="DNAWorks compatible output",
    )
    args = parser.parse_args()
    conserved = ParseStealth(args.infile, args.palindrome)
    header = f"{script_name} {' '.join(sys.argv[1:])}"
    writeFile(conserved, args.outfile, args.sorted, header, args.dnaWorks)
    # import os
    # parser.add_argument("--file",'-f', nargs="?", default=".", help="Path to file or directory")
    # args = parser.parse_args()

    # # Get the absolute path of the file or directory specified
    # target_path = os.path.abspath(args.file)

    # # Access the working directory
    # current_directory = os.getcwd()

    # print(f"target path {target_path}, CWD {current_directory}")

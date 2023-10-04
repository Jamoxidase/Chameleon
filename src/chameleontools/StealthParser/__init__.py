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

from ._StealthParser import *

"""
TABI - UCSC iGEM 2023

FastA Reader class to read FastA files

Modified from FastAreader written by David L. Bernick

Author: David L. Bernick
"""

# Args:
#     fname (str): file name (optional), default is None (STDIN)

#     Usage:
#     thisReader = FastAreader ('testTiny.fa')
#     for head, seq in thisReader.readFasta():
#         print (head,seq)

from ._FastAreader import *
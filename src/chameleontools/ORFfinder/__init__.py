"""
TABI - UCSC iGEM 2023
Modified from ORFfinder Class written for Winter 2022 BME160
BME160 Professor: David L. Bernick
"""

# ORF finder used to parse FastA format genome files.
# Optional parameters:
#     longestGene -> only reports longest genes in an ORF
#     minGene -> only reports genes of this size or larger
#     starts -> define set of valid start codons
#     stop -> define set fo valid stop codons

from ._ORFfinder import *
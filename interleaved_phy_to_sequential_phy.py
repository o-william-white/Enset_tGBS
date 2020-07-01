import sys
from Bio import AlignIO

inFile = sys.argv[1]
outFile = sys.argv[2]

input_handle = open(inFile, "rU")
output_handle = open(outFile, "w")

alignments = AlignIO.parse(input_handle, "phylip")

AlignIO.write(alignments, output_handle, "phylip-sequential")

output_handle.close()
input_handle.close()


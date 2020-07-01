import sys
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped

inFile = sys.argv[1]
outFile = sys.argv[2]

input_handle = open(inFile, "rU")
output_handle = open(outFile, "w")

alignments = AlignIO.parse(input_handle, "phylip", alphabet=Gapped(alphabet=Gapped(IUPAC.ambiguous_dna)))

AlignIO.write(alignments, output_handle, "nexus")

output_handle.close()
input_handle.close()


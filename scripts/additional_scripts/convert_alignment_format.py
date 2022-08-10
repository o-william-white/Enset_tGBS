import sys
import re
import os
from Bio import AlignIO

# note the stacks "#" line at the end of th alignment file needs to be removed

# directory with "populations.all.phylip" and "populations.var.phylip" input files
# output will be added to the same directory
inFile = sys.argv[1]

# if input specified without a final "/", add one
if not re.search("/$", inFile):
	inFile = inFile + "/"

# check input files present in dir
assert os.path.exists(inFile+"populations.all.phylip") and os.path.exists(inFile+"populations.var.phylip"), "Input files not present"

# convert files

# all sites to sequential phylip
print("Converting to all sites phylip to sequential phylip")
AlignIO.convert(inFile+"populations.all.phylip", "phylip", inFile+"populations.all.seq.phylip", "phylip-sequential")

# all sites to fasta
print("Converting to all sites phylip to fasta")
AlignIO.convert(inFile+"populations.all.phylip", "phylip", inFile+"populations.all.fasta", "fasta")

# all sites to nexus
print("Converting to all sites phylip to nexus")
AlignIO.convert(inFile+"populations.all.phylip", "phylip", inFile+"populations.all.nexus", "nexus", "DNA")

# var sites to nexus
print("Converting to variable sites phylip to nexus")
AlignIO.convert(inFile+"populations.var.phylip", "phylip-relaxed", inFile+"populations.var.nexus", "nexus", "DNA")

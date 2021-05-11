
from Bio import SeqIO

input_file  = "convert_file_formats_80_single_snp/populations.all.seq.phylip"
id_file     = "popmap_ids_without_semi_domesticated.txt"
output_file = "convert_file_formats_80_single_snp/populations.all.seq.without.semi.domesticated.phylip"

with open(id_file) as id_handle:
    wanted = set(line.rstrip("\n").split(None, 1)[0] for line in id_handle)
print("Found %i unique identifiers in %s" % (len(wanted), id_file))

records = (r for r in SeqIO.parse(input_file, "phylip-sequential") if r.id in wanted)
count = SeqIO.write(records, output_file, "phylip-sequential")
print("Saved %i records from %s to %s" % (count, input_file, output_file))
if count < len(wanted):
    print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))

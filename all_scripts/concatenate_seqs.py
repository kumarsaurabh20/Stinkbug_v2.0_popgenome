import sys
import re
import os
from Bio import SeqIO

chromosome = ["Backbone", "contig", "scaffold"]
input = sys.argv[1]
base = os.path.splitext(input)[0]
output = base + "mod.psmcfa"
f = open(output, 'w')

for each in chromosome:
	sequence = ''
	for s in SeqIO.parse(input, "fasta"):
		if re.search(each, s.id):
			sequence += s.seq.strip()
	f.write(">{}\n".format(each))
	f.write("{}\n".format(sequence))

f.close()

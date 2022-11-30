"""
The script takes a fasta file and outputs the sequence ID and and sequence length in two columns.

"""

from __future__ import print_function
from Bio import SeqIO
import sys,getopt

__author__= 'Kumar'
__desc__ = 'This script outputs a 2 column file with sequence ID in one cloumn and its length in other column'

fasta_file = ""
sep = ""
block = 1

try:
	myopts, args = getopt.getopt(sys.argv[1:], "f:p:s:b:")
	for o, a in myopts:
		if o == "-f":
			fasta_file = str(a)
		elif o == "-s":
			sep = str(a)
		elif o == "-b":
			block = int(a)
except getopt.GetoptError as e:
	print(str(e))
	print("Usage: %s -f <fasta file> -s <header delimiter> -b <header block>")
	print("Defaults: -s '|' -b 1")
	sys.exit(2)
	#Delimiter could be any separator symbol which is used in the sequence header
	#block is the number of column which is separated by a delimiter			

try:
	for s in SeqIO.parse(fasta_file, "fasta"):
		#temp = s.id.split("")[0]
		temp = str(s.id)
		print("%s\t%i"%(temp,len(s.seq)))
except IOError:
	print("Usage: %s -f <fasta file> -s <header delimiter> -b <header block>" %(sys.argv[0]))
	print("Defaults: -s '|' -b 1")
	sys.exit(2)

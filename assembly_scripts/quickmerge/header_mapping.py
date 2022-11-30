#!/usr/bin/env python2
from __future__ import (print_function, division)
from Bio import SeqIO
import sys, os, textwrap, getopt

__author__ = 'Kumar'

"""
Mapping file:

scaffold1|size2271171   Scaffold00001
scaffold2|size2035591   Scaffold00002
scaffold3|size2019171   Scaffold00003
scaffold4|size1894588   Scaffold00004
scaffold5|size1456919   Scaffold00005

"""
fasta_file = ""
map_file = ""

try:
	myopts, args = getopt.getopt(sys.argv[1:],"f:m:")
	for o, a in myopts:
		if o == '-f':
			fasta_file = str(a)
		elif o == '-m':
			map_file = str(a)   
except getopt.GetoptError as e:
	print(str(e))
	print("Usage:: %s -f <fasta file> -m <mapping file>" % sys.argv[0])
	sys.exit(2)   

old_header = {}
try:
	file = open(map_file, "r")
	for line in file:
		temp = line.split("\t")
		old_header[temp[0]] = temp[1].strip()
	file.close
except IOError:
	print("Usage:: %s -f <fasta file> -m <mapping file>" % sys.argv[0])
	sys.exit(2)

try:
	for s in SeqIO.parse(fasta_file, "fasta"):
		for key, value in old_header.iteritems():
			if s.id == key:
				print(">%s\n%s"%(value, str(s.seq)))
except IOError:
	print("Usage:: %s -f <fasta file> -m <mapping file>" % sys.argv[0])
	sys.exit(2)

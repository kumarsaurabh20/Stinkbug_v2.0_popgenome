from Bio import SeqIO
import os
import gzip

class ReadFa:

	def __init__(self, fastq):
		self.fastq = fastq


	def fa_check(self):
		# name_tuple = os.path.splitext(os.path.basename(self.bam))
		if os.path.exists(self.fastq):
			name_tuple = os.path.splitext(os.path.basename(self.fastq))
			if name_tuple[1].lower() is "fasta":
				print("file is a FASTA file!")
				return False
			elif name_tuple[1].lower() is "fastq":
				print("file is a FASTQ file!")
				return True
			else:
				print("Unknown file extension!!")
				exit()
		else:
			print("BAM file was not found!")
			return False

	
	def count_reads(self):
		temp = []	
		file = gzip.open(self.fastq,"rt")
		for s in SeqIO.parse(file, "fastq"):
			temp.append(s.id)
		return len(temp)

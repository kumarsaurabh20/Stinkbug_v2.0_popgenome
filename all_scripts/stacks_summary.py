import os
import re
import sys
import pandas as pd
from bam_summary import Samstats
from matplotlib import pyplot as plt

"""
This script first parses stacks process_radtag.log file and returns all the details(including stacks version, 
command and sequencing run or global summary) and also a per sample retained reads summary. Additionally it also 
takes pstacks log file and returns exploratory graphs about the quality of sample alignment. 

Author: Kumar Saurabh Singh (ks575@exeter.ac.uk)
"""

sys.path.append(".")

rx_dict = {'version': re.compile(r'(process_radtags\sv.+)\n'), 'command': re.compile(r'(process_radtags\s-p.+)\n'), 'run_summary': re.compile(r'(.+\.[fastq]+\.gz)\t(\d{0,})\t(\d{0,})\t(\d{0,})\t(\d{0,})\t(\d{0,})\n'), 'sample_summary': re.compile(r'([A-Z]+)\t([A-Za-z0-9]+)\t(\d{0,})\t(\d{0,})\t(\d{0,})\t(\d{0,})\n')}


def _parse_line(line):
	"""
	Do a regex search against all defined regexes and
	return the key and match result of the first matching regex

	"""
	for key, rx in rx_dict.items():
		match = rx.search(line)
		if match:
			return key, match
		# if there are no matches
	return None, None


def parse_file(filepath):

	version = ''
	command = ''
	rundata = []  # create an empty list to collect the data
	sampledata = []
	# open the file and read through it line by line
	with open(filepath, 'r') as file_object:
		line = file_object.readline()
		while line:
			# at each line check for a match with a regex
			key, match = _parse_line(line)

			# extract stacks version name
			if key == 'version':
				version = match.group()

			# extract command used to process barcodes
			if key == 'command':
				command = match.group()

			# identify a sequencing run summary 
			if key == 'run_summary':
				# extract sequencing run summary
				run_summary = match.group()
			
				# read each line of the table until a blank line
				while line.strip():
					print(line)
					# extract number and value
					readfile, ret_reads, lowQ, ambiguous_bc, ambiguous_radtag, total = line.strip().split('\t')
					total = total.strip()
					# create a dictionary containing this row of data
					# print(readfile)
					# print(total)
					runsum = {
						'filename': readfile,
						'retained_reads': int(ret_reads),
						'low_quality': int(lowQ),
						'barcode': int(ambiguous_bc),
						'radtag': int(ambiguous_radtag),
						'total': int(total)
						}
				# append the dictionary to the data list
					rundata.append(runsum)
					line = file_object.readline()

			if key == 'sample_summary':
				sample_summary = match.group()
				line = file_object.readline()
				
				while line.strip():
					print(line)
					barcode, sample, total, drop_rt, drop_qual, retained = line.strip().split('\t')
					retained = retained.strip()

					samplesum = {
							'barcode': barcode,
							'sample': sample,
							'total': int(total),
							'drop_radtag': int(drop_rt),
							'drop_qual': int(drop_qual),
							'retained': int(retained)
						}
				
					sampledata.append(samplesum)
					line = file_object.readline()	

			line = file_object.readline()

		# create a pandas DataFrame from the list of dicts
		data = pd.DataFrame(sampledata, columns = ['barcode', 'sample', 'total', 'drop_radtag', 'drop_qual', 'retained'])
		coveragelist = estimate_cov_for_all_samples(data)
		data['bam_coverage'] = coveragelist
		data['map_percentage'] = (data.bam_coverage / data.retained) * 100
		export_df(data, "tsv")
		# plotdata = data.iloc[:, [1, 7]]
		# print(plotdata)
		# plotdata.set_index('sample', inplace=True)
		# plotdata.plot(kind="bar")
	
	return version, command, rundata, data


def plot_summary(data):
	plt.style.use('ggplot')
	plotdata = data.loc[:,['sample', 'map_percentage']]
	plotdata.set_index('sample', inplace=True)
	plotdata.plot(kind="bar")
	fig = plotdata.plot(kind='bar', figsize=(70, 30), fontsize=20).get_figure()
	fig.savefig('bugs_summary_plot.pdf')


def read_summary(summaryfile):
	data = pd.read_csv(summaryfile, sep="\t")
	coveragelist = estimate_cov_for_all_samples(data)
	data['bam_coverage'] = coveragelist
	data['map_percentage'] = (data.bam_coverage / data.tot_cov) * 100
	export_df(data, "tsv")
	return data


def export_df(datafile, type):
	if type == "csv":
		datafile.to_csv("bugs.summary.csv", index=False, float_format='%.2f')
	elif type == "tsv":
		datafile.to_csv("bugs.summary.tsv", index=False, float_format='%.2f', sep="\t")
	elif type == "xls":
		datafile.to_excel("bugs.summary.xlsx")
	else:
		print("File type is not provided!")


def estimate_mapped_reads(bam):
	new_bam = Samstats(bam)
	if new_bam.sam_check():
		print(new_bam.mapped_reads())
		return new_bam.mapped_reads()
	else:
		print("BAM file path not found!")


def estimate_cov_for_all_samples(pddata):
	cov = []
	samples = pddata.iloc[:,0].tolist()
	print(samples)
	#print(f'all samples in the rad-seq data \n {samples}')
	for sample in samples:
		path = "./bugs/"
		filename = path + str(sample) + ".bam"
		cov.append(estimate_mapped_reads(filename))
	return cov


if __name__ == '__main__':
	#filepath = 'bugs2.process_radtags.Files.log'
	filepath = 'process_radtags.raw_sequences.log'
	summary = 'bugs_reads_counts.txt'
	version, command, runsummary, data = parse_file(filepath)
	# print(version)
	# print(command)
	# print(summary)
	#data = read_summary(summary)
	print(data)
	##plot_summary(data)
	# print(estimate_mapped_reads("./bams/GOLB1.bam"))
	# all_bams = os.listdir("./bugs")
	# for a in all_bams:
	#	if a.endswith(".bam"):
	#		print(a)
	#	else:
	#		pass

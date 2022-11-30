import os, sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from bam_summary import Samstats
import gzip

def estimate_mapped_reads(bam):
	new_bam = Samstats(bam)
	if new_bam.sam_check():
		print(new_bam.mapped_reads())
		return new_bam.mapped_reads()
	else:
		print("BAM file path not found!")

sys.path.append(".")
data = pd.read_csv("read_counts_bugs.txt", sep="\t")
print(data)
count = []
samples = data.iloc[:,0].tolist()
print(f'all samples in the rad-seq data \n {samples}')
for sample in samples:
	path = "./bugs/"
	filename = path + str(sample) + ".bam"
	count.append(estimate_mapped_reads(filename))
data['bam_coverage'] = count
data['map_percentage'] = (data.bam_coverage / data.tot_cov) * 100
data.to_csv("bugs.summary.tsv", index=False, float_format='%.2f', sep="\t")
plt.style.use('ggplot')
plotdata = data.loc[:,['sample', 'map_percentage']]
plotdata.set_index('sample', inplace=True)
plotdata.plot(kind="bar")
fig = plotdata.plot(kind='bar', figsize=(70, 30), fontsize=20).get_figure()
fig.savefig('bugs_new_summary_plot.pdf')

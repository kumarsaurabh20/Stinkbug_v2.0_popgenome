import pandas as pd


##df[df.column.str.contains("word")]

set1 = pd.read_csv("popmap_set1.txt", sep="\t", header=None)
print(set1)
set2 = pd.read_csv("popmap_set2.txt", sep="\t", header=None)
merge = pd.read_csv("bugs.merged.filtered.txt", sep="\t", header=None)

frames = [set1, set2]
sets = pd.concat(frames)
print(merge)
print(sets)
newfile = open("filtered_sets.txt", "a+")
for i,r in merge.iterrows():
	for x,y in sets.iterrows():
		#print(f"{y}")
		if str(r[0]) == y[0]:
			newfile.write(f"{y[0]}\t{y[1]}\n")
		else:
			pass
newfile.close()
			 

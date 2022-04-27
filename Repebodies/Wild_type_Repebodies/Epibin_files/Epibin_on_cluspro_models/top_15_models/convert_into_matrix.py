#!/usr/bin/python
import pandas as pd

user_input = input("Enter output file name:\n")
output_file_name = user_input.split(".csv")[0].strip()
IFH1 = open(output_file_name+".csv", "r")
lines = IFH1.readlines()
line_count = 0
overlap_dict = {}
for line in lines:
	if line_count > 0:
		line = line.split(",")
		Ab_1 = line[0].strip()
		Ab_2 = line[1].strip()
		score = line[2].strip()
		if Ab_1 not in overlap_dict:
			overlap_dict[Ab_1] = {}
			if Ab_2 not in overlap_dict[Ab_1]:
				overlap_dict[Ab_1][Ab_2] = score
		else:
			if Ab_2 not in overlap_dict[Ab_1]:
				overlap_dict[Ab_1][Ab_2] = score

	line_count = line_count+1


data = overlap_dict
df = pd.DataFrame.from_dict(data, dtype="object")
df.to_csv(output_file_name+'_bins.csv', header = True)
print(df)

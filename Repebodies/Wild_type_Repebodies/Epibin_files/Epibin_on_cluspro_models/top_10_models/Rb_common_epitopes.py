#/usr/bin/env/python
#from __future__ import division
import os, sys
import csv


def common_epitopes(file1, file2):
	
	file1 = file1.strip()
	IFH1 = open(file1, "r")
	repebody_1 = file1.split("_epitopes.csv")[0]
	file2 = file2.strip()
	IFH2 = open(file2, "r")
	repebody_2 = file2.split("_epitopes.csv")[0]
	print(repebody_2)
	print(repebody_1)
	
	output_file_1 = open(repebody_1+"_"+repebody_2+"_score.csv", "w")
	output_file_2 = open(repebody_1+"_"+repebody_2+"_common_epitopes.csv", "w")
	output_file_2.write("%s,%s,Common Epitopes,\n"%(repebody_1, repebody_2))
	csv_writer = csv.writer(output_file_1)
	lines2 = IFH2.readlines()
	top_row = []
	top_row.append(" ")
	for line2 in lines2:
		line2 = line2.strip("\n")
		line2 = line2.split(",")
		line2[0] = line2[0].strip()
		if line2[0] not in top_row:
			top_row.append(line2[0])
	csv_writer.writerow(top_row)

	lines1 = IFH1.readlines()
	for line1 in lines1:
		row = []
		epitope_score = 0
		Ag_epitope_model1 = set()
		line1 = line1.strip("\n")
		#print(line1)
		line1 = line1.split(",")
		print(line1)
		print(len(line1))
		total_epitopes_1 = len(line1)-1
		print(total_epitopes_1)
		row.append(line1[0])
		count_1 = 0
		for residues_model1 in line1:
			if count_1 != 0: #ensures only epitopes and not model names are included in the list
				if residues_model1 != '':
					residues_model1 = residues_model1.strip("\r")
					print(residues_model1)
					if residues_model1 not in Ag_epitope_model1:
						Ag_epitope_model1.add(residues_model1)
			count_1 = count_1+1

		for line2 in lines2:
			Ag_epitope_model2 = set()
			line2 = line2.strip("\n")
			line2 = line2.split(",")
			total_epitopes_2 = len(line2)-1


			count_2 = 0

			for residues_model2 in line2:
					
				if count_2 !=0:
					if residues_model2 != '':
					
						residues_model2 = residues_model2.strip("\r")
						if residues_model2 not in Ag_epitope_model2:
							Ag_epitope_model2.add(residues_model2)
				count_2 = count_2+1			

			common_epitopes = Ag_epitope_model1.intersection(Ag_epitope_model2)
			total_epitopes = Ag_epitope_model1.union(Ag_epitope_model2)
			epitope_score = len(common_epitopes)
			output_file_2.write("%s,%s,%s\n"%(line1[0], line2[0], common_epitopes))
			numerator = epitope_score 
			denominator = len(total_epitopes)
			normalized_epitope_score = float(numerator)/float(denominator)
			normalized_epitope_score = format(normalized_epitope_score, '.3f')
			row.append(normalized_epitope_score)
		csv_writer.writerow(row)
	IFH1.close()
	IFH2.close()
	output_file_1.close()
	output_file_2.close()


IFH3 = open("list_of_epitope_files", "r")
file_names = IFH3.readlines()
list1 = []
list2 = []

for file_name in file_names:
	file_name = file_name.strip()
	if file_name not in list1:
		list1.append(file_name)
	if file_name not in list2:
		list2.append(file_name)
count1 = 0
for file1 in list1:
	count2 = 0
	for file2 in list2:
		if count2 > count1:
			common_epitopes(file1, file2)
		elif count2 == count1:
			common_epitopes(file1, file2)
		count2 = count2 + 1
	count1 = count1 + 1
	print(count1)
	print(count2)






				



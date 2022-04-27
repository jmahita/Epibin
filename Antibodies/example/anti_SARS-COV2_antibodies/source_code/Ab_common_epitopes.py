#!usr/bin/python
import os
import sys
import csv


def common_epitopes(file1, file2):
	
	file1 = file1.strip()
	IFH1 = open(file1, "r") 
	epitopes_1 = csv.reader(IFH1, delimiter = ",")
	list_epitopes_1 = list(epitopes_1)
	print(epitopes_1)
	print(list_epitopes_1)
	antibody_1 = file1.split("_epi.csv")[0]
	file2 = file2.strip()
	IFH2= open(file2, "r")
	epitopes_2 = csv.reader(IFH2, delimiter = ",")
	antibody_2 = file2.split("_epi.csv")[0]
	print(antibody_2)
	print(antibody_1)
	
	output_file_1 = open(antibody_1+"_"+antibody_2+"_score.csv", "w")
	csv_writer_1 = csv.writer(output_file_1, delimiter =",")
	output_file_2 = open(antibody_1+"_"+antibody_2+"_common_epitopes.csv", "w")
	csv_writer_2 = csv.writer(output_file_2, delimiter = ",")
	of2_row = []
	of2_row.append(antibody_1)
	of2_row.append(antibody_2)
	of2_row.append("Common_Epitopes")
	csv_writer_2.writerow(of2_row)
	top_row = []
	top_row.append(" ")
	for row2a in epitopes_2:
		row2a[0] = row2a[0].strip()

		if row2a[0] not in top_row:
			top_row.append(row2a[0])
	csv_writer_1.writerow(top_row)
	IFH2.close()
	IFH3 = open(file2, "r")
	epitopes_3 = csv.reader(IFH3, delimiter = ",")
	list_epitopes_2 = list(epitopes_3)
	for row1 in list_epitopes_1:
		row = []
		epitope_score = 0
		Ag_epitope_model1 = set()
		total_epitopes_1 = len(row1) - 1
		row.append(row1[0])
		count_1 = 0
		for residues_model1 in row1:
			if count_1 != 0: #ensures only epitopes and not model names are included in the list
				if residues_model1 != '':

					if residues_model1 not in Ag_epitope_model1:
						Ag_epitope_model1.add(residues_model1)

			count_1 = count_1 + 1

		for row2b in list_epitopes_2:

			Ag_epitope_model2 = set()

			total_epitopes_2 = len(row2b) - 1

			count_2 = 0

			for residues_model2 in row2b:
					
				if count_2 != 0:
					if residues_model2 != '':

						if residues_model2 not in Ag_epitope_model2:
							Ag_epitope_model2.add(residues_model2)

				count_2 = count_2 + 1			

			common_epitopes = Ag_epitope_model1.intersection(Ag_epitope_model2)
			epitope_score = len(common_epitopes)
			total_epitopes = Ag_epitope_model1.union(Ag_epitope_model2)

			of2_row1 = []
			of2_row1.append(row1[0])
			of2_row1.append(row2b[0])
			if epitope_score > 0:

				for epi_res in common_epitopes:
					of2_row1.append(epi_res)
			else:
				of2_row1.append("No common_epitopes")

			csv_writer_2.writerow(of2_row1)
			
			numerator = epitope_score 

			denominator = len(total_epitopes)

			normalized_epitope_score = float(numerator)/float(denominator)
			normalized_epitope_score = format(normalized_epitope_score, '.3f')

			row.append(normalized_epitope_score)

			print(row1[0], row2b[0], normalized_epitope_score)
		csv_writer_1.writerow(row)
	IFH1.close()

	IFH3.close()
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
#count2 = 0
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








				



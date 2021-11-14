#!/usr/bin/env/python
import os, sys
import csv

def number_of_models(threshold, file_name):
	input_file = open(file_name, "r")
	repebody_pair = file_name.split("_")
	lines1 = csv.reader(input_file)
	threshold = threshold.strip()
	count_row = 0
	total_rows = 0
	print(repebody_pair[0])
	print(repebody_pair[1])
	for line1 in lines1:
		element_count = 0
		for score in line1:
			score = score.strip("''")
			print(score)
			if element_count > 0 and score[0] != "f":
				if float(score) >= float(threshold):
					count_row = count_row + 1
					print(count_row)
					break
			element_count = element_count + 1
			print(element_count)
		total_rows = total_rows + 1

	print("this is the count")
	print(count_row)
	print("total_no_of_models")
	print(int(total_rows)-1)
	
	normalized_value_row = format((float(count_row)/(int(total_rows)-1))*100, ".2f")
	return(normalized_value_row)
	input_file.close()

def averaged(threshold, file_name):
	input_file_1 = open(file_name, "r")
	transpose_file = file_name.split(".csv")[0]
	input_file_2 = open(transpose_file+"_T.csv", "r")
	repebody_pair = file_name.split("_")
	normalized_value_row = number_of_models(threshold, file_name)
	normalized_value_column = number_of_models(threshold, file_name)
	averaged_value = (float(normalized_value_row) + float(normalized_value_column))/2
	print(averaged_value)
	epibin_score = (100 - averaged_value)/100
	return(repebody_pair[0], repebody_pair[1], epibin_score)
	input_file_1.close()
	input_file_2.close()



user_input_1 = input("Enter threshold value (between 0 and 1):\n")
threshold_EOS = user_input_1.strip()
user_input_1 = threshold_EOS.strip(".")
IFH2 = open("list_common_epitope_scores", "r")
output_file = open("repebody_epibin_EOS_"+user_input_1+".csv", "w")
csv_writer = csv.writer(output_file)
row1 = []
row1.append("Repebody 1")
row1.append("Repebody 2")
row1.append("Epibin Score")
csv_writer.writerow(row1)
lines = IFH2.readlines()
for line in lines:
	line = line.strip("\n")
	row2 = averaged(threshold_EOS, line)
	csv_writer.writerow(row2)
IFH2.close()
output_file.close()



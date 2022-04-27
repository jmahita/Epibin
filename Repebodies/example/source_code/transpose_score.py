import os, sys
import csv
import pandas as pd

IFH1 = open("list_common_epitope_scores" , "r")
lines = IFH1.readlines()
for file in lines:
	file = file.strip("\n")
	file_name = file.split(".csv")[0]
	pd.read_csv(file, header=None).T.to_csv(file_name+'_T.csv', header=False, index=False)
IFH1.close()
#!/usr/bin/python

import os, sys
Ag_epitope = set()
repebody_clone = input("Enter name of repebody clone:\n")
repebody_clone = repebody_clone.strip()
IFH1 = open("list_det_cont", "r")
file_names = IFH1.readlines()
OFH1 = open(repebody_clone+"_epitopes.csv", "w")
for name in file_names:
	Ag_epitope = set()
	name = name.strip("\n")
	IFH2 = open(name, "r")
	lines = IFH2.readlines()
	count = 0
	for line in lines:
		line.strip("\n")
		line = line.split(",")
		print(line)
		if count > 0:
			#print line[0]+line[1]
			if line[1] != "\n":
				if line[0] not in Ag_epitope:
					Ag_epitope.add(str(line[0]+" "+line[1]))
		count = count+1

	print(Ag_epitope)
	OFH1.write("%s,"%(name))
	for i in Ag_epitope:
		OFH1.write("%s,"%(i))
	Ag_epitope.clear()
	IFH2.close()
	OFH1.write("\n")
IFH1.close()


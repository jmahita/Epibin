#!/usr/bin/python

import os, sys
import csv
from epitope import *
from Bio import PDB
from Bio.PDB import Polypeptide
from Bio.PDB.Polypeptide import *

antibody_id = raw_input("Enter name of antibody:\n")
antibody_id = antibody_id.strip()
H_chain = raw_input("Chain ID of heavy chain:")
H_chain = H_chain.strip()
L_chain = raw_input("Chain ID of light chain:")
L_chain = L_chain.strip()
Ag_chain = raw_input("Chain ID of antigen chain:")
Ag_chain = Ag_chain.strip()
Cutoff = raw_input("Enter cutoff value:")
Cutoff = int(Cutoff.strip())
csvfile = open(antibody_id+"_epi.csv", "w")
csvwriter = csv.writer(csvfile, delimiter = ",")
IFH1 = open("list_docked_models", "r")
docked_models = IFH1.readlines()
for model in docked_models:
	model = model.strip("\n")
	epi_res = epitope(model,Ag_chain, H_chain, L_chain, Cutoff)
	epi_res_nos = []
	epi_res_nos.append(model)
	for residue in epi_res:
		res_id = residue.get_id()
		res_no = res_id[1]
		res_no = str(res_no)
		res_name = residue.get_resname()
		one_letter = three_to_one(res_name)
		epi_res_id = res_no+one_letter
		print(epi_res_id)
		epi_res_nos.append(epi_res_id)

	print(epi_res_nos)
	#row = " ".join(epi_res_nos)
	#print(row)
	csvwriter.writerow(epi_res_nos)

IFH1.close()
csvfile.close()



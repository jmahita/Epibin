import os, sys, tarfile, csv
import numpy
from glob import glob
import math
from Bio import PDB
from Bio.PDB import *
from contact_2 import *

aa_categories = {1 :['GLY', 'ALA', 'VAL', 'LEU', 'MET', 'ILE'],
				2 : ['SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN'],
				3 : ['LYS', 'ARG', 'HIS'],
				4 : ['ASP', 'GLU'],
				5 : ['PHE', 'TYR', 'TRP']}
"""				
def extract_chain_length(dock_model,chain_ID):
	res = 0
	non_res = 0
	pdb_io = PDB.PDBIO()
	pdb_parser = PDB.PDBParser()
	pdbfile = dock_model
	structure = pdb_parser.get_structure(" ", pdbfile)
	for model1 in structure:
		model_id = model1.get_id()	
		for chain1 in model1:
			chain_id_1 = chain1.get_id()
			print(chain_id_1) 
			if chain_id_1 == chain_ID:
				residues1 = chain1.get_list()
				#print(residues1)
				for residue in residues1:
					residue_id = residue.get_id()
					if residue_id[0] == " ":
						print(residue_id[1])
						res = res+1
					else:
						non_res = non_res+1


	return(res)
"""
#output = extract_chain_length("3mac.pdb", "H")
#print(output)

def Ab_Ag_contact_analysis(dock_list,Ag_chain_ID, Ab_chain_H, Contact_cut):	
	Ab_contact_frequency = {}
	Ab_contact_distribution = {}
	Ag_contact_frequency = {}
	Ag_contact_distribution = {}
	dock_models = open(dock_list, "r").readlines()
	file_count = 0
	OFH1a = open("Ab-chHL_contact_frequency.csv", "w")
	#OFH1b = open("Ab-chL_contact_frequency.csv", "w")
	OFH2 = open("Ag_contact_frequency.csv", "w")
	OFH3a = open("Ab-chHL_aatype_distribution.csv","w")
	#OFH3b = open("Ab-chL_aatype_distribution.csv","w")
	OFH4 = open("Ag_aatype_distribution.csv" ,"w")
	OFH5a = open("Ab-chHL_contact_distribution.csv", "w")
	#OFH5b = open("Ab-chL_contact_distribution.csv", "w")
	OFH6 = open("Ag_contact_distribution.csv", "w")
	for dock_model in dock_models:
		dock_model = dock_model.strip("\n")
		base_name = dock_model.strip(".pdb")
		#print(dock_model)
		output=get_contacts(dock_model, Ag_chain_ID, Ab_chain_H, Contact_cut)
		Ab_contacts = {}
		Ag_contacts = {}
		file_count = file_count+1
	#get_contacts(dock_model, "A", "H", "L", "5")
	#print(output)
		for line1 in output:
			line2_count = 0
			for line2 in line1:
				#print(line2)
				Ab_res1a = line2[3].strip()
				Ab_res1b = line2[5].strip()
				Ab_chainid = line2[5].strip()
				Ab_resno = Ab_res1a+Ab_res1b
				Ab_resname = line2[4].strip()
				Ag_resno = str(line2[0])
				Ag_resname = line2[1].strip()
				Ag_atom = line2[2]

				#######CONTACT DICTIONARY FOR ANTIBODY RESIDUES###########
				if Ab_resno in Ab_contacts:
					if Ag_resno in Ab_contacts[Ab_resno]:
						Ab_contacts[Ab_resno][Ag_resno].append(Ag_atom) #Number of Ab_resname values indicates how many atoms of the Ag residue makes contact with the specific Ab_residue
					else:
						Ab_contacts[Ab_resno][Ag_resno] = []
						Ab_contacts[Ab_resno][Ag_resno].append(Ag_atom)
				else:
					Ab_contacts[Ab_resno] = {}
					Ab_contacts[Ab_resno][Ag_resno] = []
					Ab_contacts[Ab_resno][Ag_resno].append(Ag_atom)
				

				#########CONTACT DICTIONARY FOR ANTIGEN RESIDUES############
				if Ag_resno in Ag_contacts:
					if Ab_resno in Ag_contacts[Ag_resno]:
						Ag_contacts[Ag_resno][Ab_resno].append(Ag_atom)
					else:
						Ag_contacts[Ag_resno][Ab_resno] = []
						Ag_contacts[Ag_resno][Ab_resno].append(Ag_atom)
				else:
					Ag_contacts[Ag_resno] = {}
					Ag_contacts[Ag_resno][Ab_resno] = []
					Ag_contacts[Ag_resno][Ab_resno].append(Ag_atom)

				line2_count = line2_count + 1
		#print(Ab_contacts)
		#print(Ag_contacts)
         
		###########CALCULATION OF THE FREQUENCIES OF EACH ANTIBODY RESIDUE BEHAVING AS A CONTACT OVERALL AND OF BEING AN ANTIGEN RESIDUE-SPECIFIC CONTACT #################
		for Ab_contact in Ab_contacts:
			if Ab_contact in Ab_contact_frequency:
				Ab_contact_frequency[Ab_contact] = Ab_contact_frequency[Ab_contact] + 1
				for Ag_contact in Ab_contacts[Ab_contact]:
					if Ag_contact in Ab_contact_distribution[Ab_contact]:
						Ab_contact_distribution[Ab_contact][Ag_contact] = Ab_contact_distribution[Ab_contact][Ag_contact]+1
					else:
						Ab_contact_distribution[Ab_contact][Ag_contact] = 0
						Ab_contact_distribution[Ab_contact][Ag_contact] = Ab_contact_distribution[Ab_contact][Ag_contact]+1

			else:
				Ab_contact_frequency[Ab_contact] = 0
				Ab_contact_frequency[Ab_contact] = Ab_contact_frequency[Ab_contact] + 1
				Ab_contact_distribution[Ab_contact] = {}
				for Ag_contact in Ab_contacts[Ab_contact]:
					Ab_contact_distribution[Ab_contact][Ag_contact] = 0
					Ab_contact_distribution[Ab_contact][Ag_contact] = Ab_contact_distribution[Ab_contact][Ag_contact]+1

		###########CALCULATION OF THE FREQUENCY OF EACH ANTIGEN RESIDUE BEING A CONTACT OVERALL AS WELL AS WITH A SPECIFIC ANTIBODY RESIDUE#################
		for Ag_contact in Ag_contacts:
			if Ag_contact in Ag_contact_frequency:
				Ag_contact_frequency[Ag_contact] = Ag_contact_frequency[Ag_contact] + 1
				for Ab_contact in Ag_contacts[Ag_contact]:
					if Ab_contact in Ag_contact_distribution[Ag_contact]:
						Ag_contact_distribution[Ag_contact][Ab_contact] = Ag_contact_distribution[Ag_contact][Ab_contact]+1
					else:
						Ag_contact_distribution[Ag_contact][Ab_contact] = 0
						Ag_contact_distribution[Ag_contact][Ab_contact] = Ag_contact_distribution[Ag_contact][Ab_contact]+1
			else:
				Ag_contact_frequency[Ag_contact] = 0
				Ag_contact_frequency[Ag_contact] = Ag_contact_frequency[Ag_contact] + 1
				Ag_contact_distribution[Ag_contact] = {}
				for Ab_contact in Ag_contacts[Ag_contact]:
					Ag_contact_distribution[Ag_contact][Ab_contact] = 0
					Ag_contact_distribution[Ag_contact][Ab_contact] = Ag_contact_distribution[Ag_contact][Ab_contact]+1
	
	#print(file_count)
	print(Ab_contact_distribution)	
	
	########GENERATION OF DICTIONARIES OF RESIDUE NUMBER-RESIDUE NAME FOR ANTIGEN AND ANTIBODY###########################
	Antigen_residues = extract_residues(dock_model,Ag_chain_ID)
	Antigen_res_dict = {}
	Ab_res_dict = {}
	for antigen_res in Antigen_residues:
		antigen_resid = antigen_res.get_id()
		antigen_resname = antigen_res.get_resname()
		antigen_resno = str(antigen_resid[1])+str(antigen_resid[2]).strip()
		if antigen_resno not in Antigen_res_dict:
			Antigen_res_dict[antigen_resno] = antigen_resname
	print(Antigen_res_dict)

	#############CALCULATION OF PERCENT NORMALIZED FREQUENCY OF EACH RESIDUE, IN ANTIBODY CHAIN H, OF BEING A CONTACT#########################
	#############CALCULATION OF AATYPE DISTRIBUTION OF ANTIGEN CONTACTS MADE BY EACH RESIDUE IN ANTIBODY CHAIN H#########################
	

	OFH1a.write("Ab_Hchain_ResNo,Frequency of being a contact,Number of Models,Normalized Frequency,Percent Normalized Frequency,\n")
	OFH3a.write("Ab_Hchain_ResNo, Non-Polar Aliphatic, Polar Uncharged, Positively Charged, Negatively Charged, Aromatic,\n")
	OFH5a.write("Ab_Hchain_ResNo,Ag residue making contact,Frequency of this contact,Normalized Frequency of this contact,Percent Normalized Frequency of Contact,\n")
	
	Ab_chainH_res = extract_residues(dock_model,Ab_chain_H)
	for resH in Ab_chainH_res:
		chains = Selection.unfold_entities(resH, "C")
		for chain1 in chains: 
			chainH_id = chain1.get_id()
		chainH_id = chainH_id
		resH_id = resH.get_id()
		resH_name = resH.get_resname()
		resH_no = str(resH_id[1])+(str(resH_id[2]).strip())+str(chainH_id)
		#print(resH_no)
		if resH_no not in Ab_res_dict:
			Ab_res_dict[resH_no] = resH_name

		if resH_no in Ab_contact_frequency.keys():
			normalized_frequency_H = float(Ab_contact_frequency[resH_no])/float(file_count)
			percent_normalized_frequency_H = normalized_frequency_H * 100
			#print(resH_no, Ab_contact_frequency[resH_no], file_count, normalized_frequency_H,percent_normalized_frequency_H)
			OFH1a.write("%s,%s,%s,%s,%s,\n"%(resH_no, Ab_contact_frequency[resH_no], file_count, normalized_frequency_H,percent_normalized_frequency_H))
		else:
			print(0,0)
			OFH1a.write("%s,0,0,0,0,\n"%(resH_no))
		

		if resH_no in Ab_contact_distribution.keys():
			print(resH_no)
			print(resH_name)
			print(Ab_contact_distribution[resH_no])
			aatype_distribution_H = {}
			OFH3a.write("%s,"%(resH_no))
			for Ag_residue in Ab_contact_distribution[resH_no]:
				normalized_value_H = float(Ab_contact_distribution[resH_no][Ag_residue])/float(file_count)
				percent_normalized_value_H = normalized_value_H * 100
				print(resH_no,Ag_residue,Ab_contact_distribution[resH_no][Ag_residue],normalized_value_H,percent_normalized_value_H)
				OFH5a.write("%s,%s,%s,%s,%s,\n"%(resH_no,Ag_residue,Ab_contact_distribution[resH_no][Ag_residue],normalized_value_H,percent_normalized_value_H))
				three_letter_aacode = Antigen_res_dict[Ag_residue]
				for aa_group in aa_categories:
					if three_letter_aacode in aa_categories[aa_group]:
						if aa_group in aatype_distribution_H:
							aatype_distribution_H[aa_group] = aatype_distribution_H[aa_group] + Ab_contact_distribution[resH_no][Ag_residue]
						else:
							aatype_distribution_H[aa_group] = 0
							aatype_distribution_H[aa_group] = aatype_distribution_H[aa_group] + Ab_contact_distribution[resH_no][Ag_residue]

			print(aatype_distribution_H)
		
			

			for aa_group in aa_categories:
				if aa_group in aatype_distribution_H:
					normalized_aatype_distribution_H = float(aatype_distribution_H[aa_group])/float(file_count)
					percent_normalized_aatype_distribution_H = normalized_aatype_distribution_H * 100
					OFH3a.write("%s,"%(percent_normalized_aatype_distribution_H))
					print(aa_group)
					print(percent_normalized_aatype_distribution_H)
				else:
					print(0)
					OFH3a.write("%s,"%(0))
			OFH3a.write("\n")
		else:
			print(resH_no,0,0,0,0)
			OFH3a.write("%s,0,0,0,0,0,\n"%(resH_no))
			OFH5a.write("%s,0,0,0,0,\n"%(resH_no))
			
	#print(resH_no)

	#############CALCULATION OF PERCENT NORMALIZED FREQUENCY OF EACH RESIDUE, IN ANTIBODY CHAIN L, OF BEING A CONTACT#########################
	#############CALCULATION OF AATYPE DISTRIBUTION OF ANTIGEN CONTACTS MADE BY EACH RESIDUE IN ANTIBODY CHAIN L#########################
	"""
	OFH1b.write("Ab_LChain_ResNo, Frequency of being a contact, Number of Models, Normalized Frequency, Percent Normalized Frequency\n")
	OFH3b.write("Ab_Lchain_ResNo, Non-Polar Aliphatic, Polar Uncharged, Positively Charged, Negatively Charged, Aromatic,\n")
	OFH5b.write("Ab_Lchain_ResNo,Ag residue making contact,Frequency of this contact,Normalized Frequency of this contact,Percent Normalized Frequency of Contact,\n")
	
	Ab_chainL_res = extract_residues(dock_model,Ab_chain_L)

	for resL in Ab_chainL_res:
		chains = Selection.unfold_entities(resL, "C")
		for chain2 in chains: 
			chainL_id = chain2.get_id()
		chainL_id = chainL_id
		resL_id = resL.get_id()
		resL_name = resL.get_resname()
		resL_no = str(resL_id[1])+(str(resL_id[2]).strip())+str(chainL_id)
		print(resL_no)
		if resL_no not in Ab_res_dict:
			Ab_res_dict[resL_no] = resL_name
		if resL_no in Ab_contact_frequency.keys():
			normalized_frequency_L = float(Ab_contact_frequency[resL_no])/float(file_count)
			percent_normalized_frequency_L = normalized_frequency_L * 100
			#print(resL_no, Ab_contact_frequency[resL_no], file_count, normalized_frequency_L,percent_normalized_frequency_L)
			print(resL_no, Ab_contact_frequency[resL_no], file_count, normalized_frequency_L,percent_normalized_frequency_L)
			OFH1b.write("%s,%s,%s,%s,%s,\n"%(resL_no, Ab_contact_frequency[resL_no], file_count, normalized_frequency_L,percent_normalized_frequency_L))
		else:
			print(0,0)
			OFH1b.write("%s,0,0,0,0,\n"%(resL_no))

		if resL_no in Ab_contact_distribution.keys():
			print(Ab_contact_distribution[resL_no])
			aatype_distribution_L = {}
			OFH3b.write("%s,"%(resL_no))
			for Ag_residue in Ab_contact_distribution[resL_no]:
				normalized_value_L = float(Ab_contact_distribution[resL_no][Ag_residue])/float(file_count)
				percent_normalized_value_L = normalized_value_L * 100
				print(resL_no,Ag_residue,Ab_contact_distribution[resL_no][Ag_residue],normalized_value_L,percent_normalized_value_L)
				OFH5b.write("%s,%s,%s,%s,%s,\n"%(resL_no,Ag_residue,Ab_contact_distribution[resL_no][Ag_residue],normalized_value_L,percent_normalized_value_L))
				three_letter_aacode = Antigen_res_dict[Ag_residue]
				for aa_group in aa_categories:
					if three_letter_aacode in aa_categories[aa_group]:
						if aa_group in aatype_distribution_L:
							aatype_distribution_L[aa_group] = aatype_distribution_L[aa_group] + Ab_contact_distribution[resL_no][Ag_residue]
						else:
							aatype_distribution_L[aa_group] = 0
							aatype_distribution_L[aa_group] = aatype_distribution_L[aa_group] + Ab_contact_distribution[resL_no][Ag_residue]
			print(aatype_distribution_L)
			for aa_group in aa_categories:
				if aa_group in aatype_distribution_L:
					normalized_aatype_distribution_L = float(aatype_distribution_L[aa_group])/float(file_count)
					percent_normalized_aatype_distribution_L = normalized_aatype_distribution_L * 100
					OFH3b.write("%s,"%(percent_normalized_aatype_distribution_L))
					print(aa_group)
					print(percent_normalized_aatype_distribution_L)
				else:
					print(0)
					OFH3b.write("%s,"%(0))
			OFH3b.write("\n")
		else:
			print(resL_no,0,0,0,0)
			OFH3b.write("%s,0,0,0,0,0,\n"%(resL_no))
			OFH5b.write("%s,0,0,0,0,\n"%(resL_no))
	"""
	#############CALCULATION OF PERCENT NORMALIZED FREQUENCY OF EACH RESIDUE, IN ANTIGEN CHAIN, OF BEING A CONTACT#########################
	#############CALCULATION OF AATYPE DISTRIBUTION OF ANTIBODY CONTACTS MADE BY EACH RESIDUE IN ANTIGEN CHAIN#########################
	OFH2.write("Ag_ResNo, Frequency of being a contact, Number of Models, Normalized Frequency, Percent Normalized Frequency,\n")
	OFH4.write("Ag_ResNo, Non-Polar Aliphatic, Polar Uncharged, Positively Charged, Negatively Charged, Aromatic,\n")
	OFH6.write("Ag_ResNo,Ab residue making contact,Frequency of this contact,Normalized Frequency of this contact,Percent Normalized Frequency of Contact,\n")

	
	for residue in extract_residues(dock_model,Ag_chain_ID):
		resAg_id = residue.get_id()
		resAg_name = residue.get_resname()
		resAg = str(resAg_id[1])+str(resAg_id[2]).strip()
		#print(antigen_resno)
		if resAg in Ag_contact_frequency.keys():
			normalized_frequency_Ag = float(Ag_contact_frequency[resAg])/float(file_count)
			percent_normalized_frequency_Ag = normalized_frequency_Ag * 100
			print(resAg, Ag_contact_frequency[resAg], file_count, normalized_frequency_Ag,percent_normalized_frequency_Ag)
			OFH2.write("%s,%s,%s,%s,%s,\n"%(resAg, Ag_contact_frequency[resAg], file_count, normalized_frequency_Ag,percent_normalized_frequency_Ag))
		else:
			print(0,0)
			OFH2.write("%s,0,0,0,0,\n"%(resAg))

		if resAg in Ag_contact_distribution.keys():
			print(Ag_contact_distribution[resAg])
			aatype_distribution_Ag = {}
			OFH4.write("%s,"%(resAg))
			for Ab_residue in Ag_contact_distribution[resAg]:
				normalized_value_Ag = float(Ag_contact_distribution[resAg][Ab_residue])/float(file_count)
				percent_normalized_value_Ag = normalized_value_Ag * 100
				print(resAg,Ab_residue,Ag_contact_distribution[resAg][Ab_residue],normalized_value_Ag,percent_normalized_value_Ag)
				OFH6.write("%s,%s,%s,%s,%s,\n"%(resAg,Ab_residue,Ag_contact_distribution[resAg][Ab_residue],normalized_value_Ag,percent_normalized_value_Ag))
				three_letter_aacode = Ab_res_dict[Ab_residue]
				for aa_group in aa_categories:
					if three_letter_aacode in aa_categories[aa_group]:
						if aa_group in aatype_distribution_Ag:
							aatype_distribution_Ag[aa_group] = aatype_distribution_Ag[aa_group] + Ag_contact_distribution[resAg][Ab_residue]
						else:
							aatype_distribution_Ag[aa_group] = 0
							aatype_distribution_Ag[aa_group] = aatype_distribution_Ag[aa_group] + Ag_contact_distribution[resAg][Ab_residue]
			print(aatype_distribution_Ag)
			for aa_group in aa_categories:
				if aa_group in aatype_distribution_Ag:
					normalized_aatype_distribution_Ag = float(aatype_distribution_Ag[aa_group])/float(file_count)
					percent_normalized_aatype_distribution_Ag = normalized_aatype_distribution_Ag * 100
					OFH4.write("%s,"%(percent_normalized_aatype_distribution_Ag))
					print(aa_group)
					print(percent_normalized_aatype_distribution_Ag)
				else:
					print(0,0)
					OFH4.write("%s,"%(0))
			OFH4.write("\n")
		else:
			print(resAg,0,0,0,0)
			OFH4.write("%s,0,0,0,0,0,\n"%(resAg))	
			OFH6.write("%s,0,0,0,0,\n"%(resAg))	



	OFH1a.close()
	#OFH1b.close()
	OFH2.close()
	OFH3a.close()
	#OFH3b.close()
	OFH4.close()
	OFH5a.close()
	#OFH5b.close()
	OFH6.close()
	print(Ab_res_dict)
	print(Ag_contact_distribution)

Ab_Ag_contact_analysis("dock_list", "B", "A", "6")


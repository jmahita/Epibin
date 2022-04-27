import os, sys, tarfile, csv
import numpy
from glob import glob
import math
from Bio import PDB
from Bio.PDB import *

def extract_residues(dock_model,chain_ID):
	
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
		
	return(residues1)
								

def extract_Ag_atoms(residue):
	resid = residue.get_id()
	res_no_1a = resid[1]
	res_no_1b = resid[2]
	res_no_1 = str(res_no_1a)+str(res_no_1b)
	res_name_1 = residue.get_resname()
	if res_name_1 != "HOH":
		atoms1 = residue.get_list()
		return(atoms1)

def extract_Ab_atoms(dock_model, chain_ID):
	pdb_io = PDB.PDBIO()
	pdb_parser = PDB.PDBParser()
	pdbfile = dock_model
	structure = pdb_parser.get_structure(" ", pdbfile)
	for model1 in structure:
		for chain in model1:
			if chain.get_id() == chain_ID:
				Ab_chain_atoms = Selection.unfold_entities(chain, "A")
				return(Ab_chain_atoms)

def get_contacts(dock_model, Ag_chain_ID, Ab_chain_HL, Contact_cut):
	Ab_chainH_atoms = extract_Ab_atoms(dock_model,Ab_chain_HL)
	#Ab_chainL_atoms = extract_Ab_atoms(dock_model,Ab_chain_L)
	Ag_reslist = extract_residues(dock_model,Ag_chain_ID)
	contact_details = {}
	Contact_cut = int(Contact_cut)
	for res in Ag_reslist:
		res_id = res.get_id()
		Ag_res_name = res.get_resname()
		Ag_res_no = res_id[1]
		if res_id[0] == " ":
			#print(res)
			Ag_chain_atom_list = extract_Ag_atoms(res)
			for atom in Ag_chain_atom_list:
				atom_center = atom.get_coord()
				ns1 = NeighborSearch(Ab_chainH_atoms)
				#ns2 = NeighborSearch(Ab_chainL_atoms)
				H_atoms = ns1.search(atom_center, Contact_cut)
				H_residue_contacts = Selection.unfold_entities(H_atoms, "R")
				#L_atoms = ns2.search(atom_center, Contact_cut)
				#L_residue_contacts = Selection.unfold_entities(L_atoms, "R")
				Ab_residue_contacts = H_residue_contacts
				for Ab_residue in Ab_residue_contacts:
					Ab_chain = Selection.unfold_entities(Ab_residue, "C")
					for chain in Ab_chain:
						Ab_chain_id = chain.get_id()
					Ab_chain_id = Ab_chain_id
					Ab_res_name = Ab_residue.get_resname()
					Ab_res_id = Ab_residue.get_id()
					if Ab_res_id[2] != " ":
						Ab_res_no_1a = Ab_res_id[1]
						Ab_res_no_1b = Ab_res_id[2]
						Ab_res_no = str(Ab_res_no_1a)+str(Ab_res_no_1b)
					else:
						Ab_res_no = str(Ab_res_id[1])

					Ab_contact_info = [Ag_res_no, Ag_res_name, atom, Ab_res_no, Ab_res_name, Ab_chain_id]

					if res not in contact_details.keys():
						contact_details[res] = [Ab_contact_info]
					else:
						if Ab_contact_info not in contact_details[res]:
							contact_details[res].append(Ab_contact_info)
	#print(contact_details)		
	Ag_Ab_contact_details = []					
	for Ag_residue in contact_details:
		#return(contact_details[Ag_residue])
		#print(Ag_residue,contact_details[Ag_residue])
		Ag_Ab_contact_details.append(contact_details[Ag_residue])	
	return(Ag_Ab_contact_details)

def VH_VL_interface_contacts(dock_model, Ab_chain_H, Ab_chain_L, Contact_cut):
	VH_VL_contact_dict = {}
	Ab_chainH_residues = extract_residues(dock_model, Ab_chain_H)
	for H_residue in Ab_chainH_residues:
		H_resname = H_residue.get_resname()
		H_resid = H_residue.get_id()
		if H_resid[2] != " ":
			H_resno_1a = H_resid[1]
			H_resno_1b = H_resid[2]
			H_resno = str(H_resno_1a)+str(H_resno_1b)
		else:
			H_resno = str(H_resid[1])
		if H_resid[0] != "W":
			Ab_chainH_atoms = extract_Ag_atoms(H_residue)
			#print(H_residue)
			#print(Ab_chainH_atoms)
	#Ab_chainH_atoms = extract_Ab_atoms(dock_model,Ab_chain_H)
		Ab_chainL_atoms = extract_Ab_atoms(dock_model,Ab_chain_L)
		Contact_cut = int(Contact_cut)
		for H_atom in Ab_chainH_atoms:
			atom_H_center = H_atom.get_coord()
			ns1 = NeighborSearch(Ab_chainL_atoms)
			L_atom_contacts = ns1.search(atom_H_center,Contact_cut)
			L_residues = Selection.unfold_entities(L_atom_contacts, "R")
			#H_residue = Selection.unfold_entities(H_atom, "R")
		
			for L_residue in L_residues:
				L_resname = L_residue.get_resname()
				L_resid = L_residue.get_id()
				if L_resid[2] != " ":
					L_resno_1a = L_resid[1]
					L_resno_1b = L_resid[2]
					L_resno = str(L_resno_1a)+str(L_resno_1b)
				else:
					L_resno = str(L_resid[1])
				VH_VL_contacts = [Ab_chain_H, H_resno, H_resname, H_atom, Ab_chain_L, L_resno, L_resname]
				if H_residue not in VH_VL_contact_dict.keys():
					VH_VL_contact_dict[H_residue] = [VH_VL_contacts]
				else:
					if VH_VL_contacts not in VH_VL_contact_dict[H_residue]:
						VH_VL_contact_dict[H_residue].append(VH_VL_contacts)

	VH_VL_contact_details = []
	for H_residue in VH_VL_contact_dict:
		VH_VL_contact_details.append(VH_VL_contact_dict[H_residue])
	#print(VH_VL_contact_details)
	return(VH_VL_contact_details)

def write_contacts(dock_list,Ag_chain_ID, Ab_chain_HL, Contact_cut):	
	dock_models = open("dock_list", "r").readlines()
	#Contact_cut = int(Cutoff)
	#list = []
	for dock_model in dock_models:
		dock_model = dock_model.strip("\n")
		base_name = dock_model.strip(".pdb")
		print(dock_model)
	#contact_details = {}
		output=get_contacts(dock_model, Ag_chain_ID, Ab_chain_HL, Contact_cut)
		#output_VH_VL=VH_VL_interface_contacts(dock_model, Ab_chain_H, Ab_chain_L, Contact_cut)
	#get_contacts(dock_model, "A", "H", "L", "5")
	#print(output)
		with open(base_name+"_det_contacts.csv", "w") as csvfile:
		#csvfile.write("%s,%s,%s,%s,%s"%("Agresnumber","Agresname","Agatom","Abresnumber","Abresname","Abchain"))
		#heading = "Agresnumber Agresname Agatom Abresnumber Abresname Abchain"
			rowwriter= csv.writer(csvfile, delimiter = ',')
		#rowwriter.writeheader()
			for line1 in output:
				for line2 in line1:
					rowwriter.writerow(line2)
					print(line2)

		csvfile.close()

		#with open(base_name+"_VH_VL_interface.csv", "w") as csvfile:
			#rowwriter= csv.writer(csvfile, delimiter = ',')
			#for line3 in output_VH_VL:
				#for line4 in line3:
					#rowwriter.writerow(line4)
					#print(line4)

		#csvfile.close()
#write_contacts("dock_list", "A", "H","L", "5")
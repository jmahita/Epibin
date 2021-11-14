import os, sys, tarfile, csv
import numpy
from glob import glob
import math
import Bio
from Bio import PDB
from Bio.PDB import *

def extract_atoms(model, chain_ID):
	pdb_io = PDB.PDBIO()
	pdb_parser = PDB.PDBParser()
	pdbfile = model
	structure = pdb_parser.get_structure(" ", pdbfile)
	for model1 in structure:
		for chain in model1:
			if chain.get_id() == chain_ID:
				#atom_list = []
				#for residue in chain:
					#res_id = residue.get_id()
					#if res_id[0] == " ":
				atoms = Selection.unfold_entities(chain, "A")
						#atom_list = atom_list+atoms
				#print(atom_list)
				return(atoms)

#extract_atoms("6azz.pdb" , "E")

def epitope(model, Ag_chain_ID, Ab_chain_H, Ab_chain_L, Contact_cut):
	Ab_chainH_atoms = extract_atoms(model,Ab_chain_H)
	Ab_chainL_atoms = extract_atoms(model,Ab_chain_L)
	Ag_atoms = extract_atoms(model, Ag_chain_ID)
	epitope = []
	Contact_cut = int(Contact_cut)
	Ab_atoms = Ab_chainH_atoms+Ab_chainL_atoms
	for Ab_atom in Ab_atoms:
		Ab_atom_id = Ab_atom.get_id()
		if Ab_atom_id != "H":
			Ab_atom_center = Ab_atom.get_coord()
			ns = NeighborSearch(Ag_atoms)
			epitope_atoms = ns.search(Ab_atom_center, Contact_cut)
			Ag_residues = Selection.unfold_entities(epitope_atoms , "R")
			for Ag_residue in Ag_residues:
				Ag_resid = Ag_residue.get_id()
				if Ag_resid[0] == " ":
					if Ag_residue not in epitope:
						epitope.append(Ag_residue)
	#for residue in epitope:
		#res_atoms = residue.get_atoms()
	#print(res_atoms)
	return(epitope)

def paratope(model, Ag_chain_ID, Ab_chain_H, Ab_chain_L, Contact_cut):
	Ab_chainH_atoms = extract_atoms(model,Ab_chain_H)
	Ab_chainL_atoms = extract_atoms(model,Ab_chain_L)
	Ag_atoms = extract_atoms(model, Ag_chain_ID)
	paratope = []
	Contact_cut = int(Contact_cut)
	Ab_atoms = Ab_chainH_atoms+Ab_chainL_atoms
	for Ag_atom in Ag_atoms:
		Ag_atom_id = Ag_atom.get_id()
		if Ag_atom_id != "H":
			Ag_atom_center = Ag_atom.get_coord()
			ns = NeighborSearch(Ab_atoms)
			paratope_atoms = ns.search(Ag_atom_center, Contact_cut)
			Ab_residues = Selection.unfold_entities(paratope_atoms, "R")
			for Ab_residue in Ab_residues:
				Ab_resid = Ab_residue.get_id()
				if Ab_resid[0] == " ":
					if Ab_residue not in paratope:
						paratope.append(Ab_residue)
	return(paratope)

#epitope_residues = epitope("6azz.pdb" ,"E", "F", "D" , "4" )
#print(epitope_residues)
def write_epitope(model, Ag_chain_ID, Ab_chain_H, Ab_chain_L, Contact_cut,output_path):
	parser = PDBParser()
	structure = parser.get_structure('Epitope' , model)
	model_name = model.split(".pdb")[0]
	#output_name = model_name.split("_")[0]
	epitope_residues = epitope(model, Ag_chain_ID, Ab_chain_H, Ab_chain_L, Contact_cut)
	epitope_dict = {}
	for epitope_residue in epitope_residues:
	#print(epitope_residue)
		epitope_resid = epitope_residue.get_id()
		epitope_full_id = epitope_residue.get_full_id()
		epitope_chain = epitope_full_id[2]
		epitope_resno = epitope_resid[1]
		epitope_resname = epitope_residue.get_resname()
		if epitope_resno not in epitope_dict.keys():
			epitope_dict[epitope_resno] = epitope_chain
	#print(epitope_resname
	#print(epitope_dict)

	class epitopeselect(Select):
		def accept_residue(self, residue):
			residue_id = residue.get_id()
			residue_full_id = residue.get_full_id()
			residue_name = residue.get_resname()
		#for epitope_residue in epitope_residues:
			#epitope_resid = epitope_residue.get_id()
			#epitope_resname = epitope_residue.get_resname()
			#if residue_id[1] == epitope_resid[1] and residue_name == epitope_resname:
		#residue_name = residue.get_resname()
			#print(epitope_residues)
			if residue_id[1] in epitope_dict.keys():
				if residue_full_id[2] == epitope_dict[residue_id[1]]:
					return True
			else:
				return False
	io = PDBIO()
	io.set_structure(structure)
	io.save(output_path+"/"+model_name+'_epi_res.pdb',epitopeselect())

def write_paratope(model, Ag_chain_ID, Ab_chain_H, Ab_chain_L, Contact_cut, output_path):
	parser = PDBParser()
	structure = parser.get_structure('Paratope' , model)
	model_name = model.split(".pdb")[0]
	#output_name = model_name.split("_")[0]
	paratope_residues = paratope(model, Ag_chain_ID, Ab_chain_H, Ab_chain_L, Contact_cut)
	paratope_dict = {}
	for paratope_residue in paratope_residues:
	#print(epitope_residue)
		paratope_resid = paratope_residue.get_id()
		paratope_full_id = paratope_residue.get_full_id()
		paratope_chain = paratope_full_id[2]
		paratope_resno = paratope_resid[1]
		paratope_resname = paratope_residue.get_resname()
		if paratope_resno not in paratope_dict.keys():
			paratope_dict[paratope_resno] = paratope_chain
	#print(epitope_resname
	#print(paratope_dict)

	class paratopeselect(Select):
		def accept_residue(self, residue):
			residue_id = residue.get_id()
			residue_full_id = residue.get_full_id()
			residue_name = residue.get_resname()
		#for epitope_residue in epitope_residues:
			#epitope_resid = epitope_residue.get_id()
			#epitope_resname = epitope_residue.get_resname()
			#if residue_id[1] == epitope_resid[1] and residue_name == epitope_resname:
		#residue_name = residue.get_resname()
			#print(epitope_residues)
			if residue_id[1] in paratope_dict.keys():
				if residue_full_id[2] == paratope_dict[residue_id[1]]:
					return True
			else:
				return False
	io = PDBIO()
	io.set_structure(structure)
	io.save(output_path+"/"+model_name+'_para_res.pdb',paratopeselect())

"""
for residue in epitope_residues:
	print(residue)
	#atoms = Selection.unfold_entities(residue, "A")
	#residue_id = residue.get_id()
	#for atom in atoms:
		#atom_coord = atom.get_coord()
		#print(residue)
		#print(atom_coord)
"""		

#write_paratope("6azz.pdb" ,"D", "F", "E" , "4")

import os, sys, tarfile, csv
import numpy
from glob import glob
import math


def identify_contacts(Cutoff, dock_list, Ag_chain, Ab_chain1, Ab_chain2):
        Contact_cut = int(Cutoff)
        AA3 = {"ALA": "A", "GLU": "E", "GLN": "Q", "ASP": "D", "ASN": "N",
        "LEU": "L", "GLY": "G", "LYS": "K", "SER": "S", "VAL": "V",
        "ARG": "R", "THR": "T", "PRO": "P", "ILE": "I", "MET": "M",
        "PHE": "F", "TYR": "Y", "CYS": "C", "TRP": "W", "HIS": "H",
        "CYX": "C", "CYD": "C", "MSE": "M", "HIP": "H", "HSC": "H"}
        #if os.path.isfile(dock_list): #list of docked models
        #    dock_models = dock_models.strip("\n")
        #else:
        #    sys.stderr.write("There is no such file: %s%s"%(dock_list, os.linesep))
        #    sys.exit()
        dock_models = open(dock_list, "r").readlines()
        for dock_model in dock_models:
            dock_model = dock_model.strip("\n")
            base_name = dock_model.strip(".pdb")
        #for dock in dock_model: #name of docked model
            IFH1 = open(dock_model).readlines() #open PDB file and read each line of the file
            #print IFH1
            Ab_pdb = []
            Ag_pdb =[]
            for line in IFH1:
                if line[:4] == "ATOM":
                    atom_name = line[12:16].strip()
                    chain_id = line[21]
                    res_num = line[22:26].strip()
                    res_name = line[17:20]
                    #line = line.split(" ")
                    #print line 
                
                    
                    if chain_id == Ab_chain1 or chain_id == Ab_chain2:
                        Ab_pdb.append(line)
                        #print "atom name=%s,chain_id =%s,res_num=%s, res_name=%s"%(atom_name, chain_id,res_num,res_name)
                    else:
                        Ag_pdb.append(line)
                        #print "atom name=%s,chain_id =%s,res_num=%s, res_name=%s"%(atom_name, chain_id,res_num,res_name)
            #print Ag_pdb
            #print Ab_pdb
            lig_in_contact = {}
            contact_details = {}
            lig_pos = []
            for l in Ag_pdb:
                if l[:4] == "ATOM":
                    atom_name_1 = l[12:16].strip()
                    chain_id_1 = l[21]
                    res_num_1 = int(l[22:26].strip())
                    res_name_1 = l[17:20]
                    res_name_1_1 = AA3[res_name_1]
                    
                    x1 = float(l[30:38].strip())
                    y1 = float(l[39:46].strip())
                    z1 = float(l[47:54].strip())
                    if atom_name_1 != "H":
                        #print atom_name_1
                        #print "%s,%s,%s\n"%(x1,y1,z1)
      
                        for r in Ab_pdb:
                            if r[:4] == "ATOM":
                                atom_name_2 = r[12:16].strip()
                                atom_name_2a = r[12:]
                                chain_id_2 = r[21]
                                res_num_2 = int(r[22:26].strip())
                                res_name_2 = r[17:20]
                                res_name_2_1 = AA3[res_name_2]
                                
                                x2 = float(r[30:38].strip())
                                y2 = float(r[39:46].strip())
                                z2 = float(r[47:54].strip())
                                if atom_name_2 != "H":
                                    #print atom_name_2
                                    #print "%s,%s,%s\n"%(x2,y2,z2)
                                #filter(lambda x: x.atom_name.count("H") == 0, Ab_pdb):
                                    #dist = numpy.linalg.norm(l.crd-r.crd)
                                    dist = math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
                                    #print dist
                                    #elm = [dist, res_num_2, res_name_2, chain_id_2]
                                    #print elm
                                    if res_num_1 == 166 and res_num_2 == 301:
                                        if atom_name_1 == "O" and atom_name_2 == "O":
                                            #print "Distance between Tyr and Met is %s\n"%(dist)
                                            #print "Atomic coordinates of O is %s, %s, %s\n"%(x1, y1, z1)
                                            print("Atomic coordinates of O is %s, %s, %s\n"%(x2, y2, z2))
                                    
                                    if dist < Contact_cut:
                                        #print "%s,%s\n"%(dist, Contact_cut)
                                        elm = [res_num_2, res_name_2_1, chain_id_2]
                                        contact_info = [res_num_1, res_name_1_1, atom_name_1, res_num_2, res_name_2_1, atom_name_2, dist]
                                        #print elm
                                        #lig_in_contact(res_num_1) = [elm]
                                        #print dist
                                        #if not lig_in_contact.has_key(res_num_1):
                                        if res_num_1 not in lig_in_contact.keys():
                                            lig_in_contact[res_num_1] = [elm]
                                            contact_details[res_num_1] = [contact_info]
                                            lig_pos.append(res_num_1)
                                            #print lig_in_contact[res_num_1]
                                            #print "Ag position:%s, Ab residues:%s\n"%(lig_pos[res_num_1],lig_in_contact[res_num_1])
                                            break
                                        else:
                                            if elm not in lig_in_contact[res_num_1]:
                                                lig_in_contact[res_num_1].append(elm)
                                            if contact_info not in contact_details[res_num_1]:
                                            	contact_details[res_num_1].append(contact_info)
                                                #print lig_in_contact[res_num_1]
    
            output_1 = csv.writer(open(base_name+"_contacts.csv", "wb"))
            output_2 = open(base_name+".det_cont.csv", "w")
            output_2.write("Ag res number, Ag res name, Atom name 1, Rb res number, Rb res name, Atom name 2, Distance(Angstroms)\n")
            for p in lig_pos:
                  #print "Model_no: %s,Res no on Ag:%s,Ab residues in contact:%s\n"%(dock_model, p, lig_in_contact[p])
                  final_result = [p, ";".join(map(lambda x: x[1], lig_in_contact[p]))]
                  #print final_result_1
                  output_1.writerow(final_result)
                  output_2.write("%s,\n"%(p))
                  for detail in contact_details[p]:
                  	#print("%s,%s,%s,%s,%s,%s,%s,%s\n")%(p,detail[0], detail[1], detail[2],detail[3],detail[4],detail[5],detail[6])
                  	output_2.write("%s,%s,%s,%s,%s,%s,%s\n"%(detail[0], detail[1], detail[2],detail[3],detail[4],detail[5],detail[6]))
             

identify_contacts(6,"list_docked_models", "D", "A", "A")

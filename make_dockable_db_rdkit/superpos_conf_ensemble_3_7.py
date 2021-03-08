#!/usr/bin/env python3

#make sure superposition in planar rings works properly, see also around 360


#adpoted => can read multimol2 file, changes solvation table if multiple molecule entries were generated
#problem: if ring is not flat, mol2db gets confused => have to be seperated (independend, if ring is used as anchor or if ring is attached to anchor)
#an issue: mol2db only allows label with 9 digits -> where to add my mutiple_ring_identifier?
#=> had to sacrify first two digits, add ring identifier at end

#make sure that the first atom in the mol2 file is part of the rigid fragment, otherwise mol2db might get confused


#looses molecules that don't have a ring!

#can not handle duplicates in file names <------------------------------

import os, sys, Mol2Writer
import chemistry3 as chemistry
from rdkit import Chem
from rdkit.Chem import AllChem
import math


#-------------------------------
def calc_rmsd (mol, atom_list,conf_list):
	ref_conf = mol.GetConformer(conf_list[0])
	ref_coords = ref_conf.GetPositions()
	sum_dist_sq = 0.0

	for conf_id in conf_list[1:] :# frist one (0) the referecne
		conf = mol.GetConformer(conf_id)
		coords = conf.GetPositions()
		for atom_id in atom_list:
			print (ref_coords[atom_id])
			print (coords[atom_id])

			x_diff = ref_coords[atom_id][0] - coords[atom_id][0]
			y_diff = ref_coords[atom_id][1] - coords[atom_id][1]
			z_diff = ref_coords[atom_id][2] - coords[atom_id][2]

			distance_sq = math.pow(x_diff,2) + math.pow(y_diff,2) + math.pow(z_diff,2)
			sum_dist_sq = sum_dist_sq + distance_sq
		rmsd = math.sqrt(sum_dist_sq)/len(atom_list)
		print (rmsd, conf_id)

#---------------------
def gen_new_title(title, counter):
    new_title = title[0:4] + title[6:] + str(counter).rjust(2)
    new_title = new_title.replace(' ', '0')
    print (title, new_title, '********')
    return new_title

#-------------------------------
def read_solv_table(solv_table_file):
	dict = {}
	for i in solv_table_file.readlines():
		if i[0:1] != ' ':
			end = i.find(' ')
			identifier = i[:end]
			dict[identifier] = [i,[]]
		else:
			dict[identifier][1].append(i)
	
	return dict
		

		
	
#---------------------
def find_atoms_in_plane_extend_ring(mol,debug):  

#get all (non H) substituents of planr ring atoms that are not part of a ring, if not planar, do not search, as their position might change depending on ring conformation
#should not be necessary as we say exactly which atom matches which one
	ri = mol.GetRingInfo()
	#print(ri.AtomRings())
	atom_list = []
	for ring in range (0, len(ri.AtomRings())): #convert tulple to list so that it can be modified
		atom_list.append([[]]) #list should also hold information about planarity
		for atmidx in ri.AtomRings()[ring]:
			atom_list[ring][0].append(atmidx)
		#print (atom_list)

	for i in range(0, len(atom_list)):
		#need to create indepent copy of list
		atom_idx_list = atom_list[i][0][:]
		#print (atom_idx_list, 'atom_idx_list')
		if chemistry.check_planar(mol, atom_idx_list,debug):
	      	  atom_list[i].append(True)
	      	  for atom in mol.GetAtoms():
	      	          #print (atom.GetIdx())
	      	          #print (atom_idx_list,  'atom_idx_list 2') 
	      	          if atom.GetIdx() in atom_idx_list: #atom is part of this ring system, check for neighbours
	      	      	    #print (atom.GetIdx(), 'checked')
	      	      	    for neighbour in atom.GetNeighbors() :
	      	      	    	#print (neighbour.GetIdx(), atom_idx_list)
	      	      	    	if neighbour.GetIdx() not in atom_idx_list and neighbour.GetAtomicNum() != 1: #do not add H
	      	      	    		#atom_list[i][0].append(neighbour.GetIdx())
	      	      	    		continue
	      	      	    		#print ('added',neighbour.GetSmarts(),neighbour.GetAtomicNum())

		else:
	      	  atom_list[i].append(False) #ring is not planar



	#print (atom_list)
	count = len(atom_list)
	return atom_list, count
#---------------------
def fuse_rings(mol,debug):  

#only check if ring is planar
#should not be necessary to extend rings as we say exactly which atom matches which one
#I think there is no need to check if the ring is planar, we can just fuse them anyway as we do an rmsd check when superposing. In case there is no good match, they will just be alinged in a different template
	ri = mol.GetRingInfo()
	#print(ri.AtomRings())
	atom_list = []
	if debug:
		print('number rings: ', len(ri.AtomRings())) 
	for ring in range (0, len(ri.AtomRings())): #convert tulple to list so that it can be modified
		atom_list.append([[]]) #list should also hold information about planarity
		for atmidx in ri.AtomRings()[ring]:
			atom_list[ring][0].append(atmidx)
		#print (atom_list)

	#for i in range(0, len(atom_list)):
		#need to create indepent copy of list
		#atom_idx_list = atom_list[i][0][:]
		#print (atom_idx_list, 'atom_idx_list')
		#if chemistry.check_planar(mol, atom_idx_list,debug):
	      	  #atom_list[i].append(True)
		#else:
	      	  ##atom_list[i].append(False) #ring is not planar

	for i in range(0, len(atom_list)):
		atom_list[i].append(False) #means ring is not planar, but for now we just need this to keep the data format consistent
	
	#check if two planar rings are fused together, this check is no longer carried out, we can fuse any rings

	remove_from_list = []
	print (atom_list)


	#find out if at least one atom is in common between two rings, if yes, fuse
	for ring_1 in range(0,len(atom_list)):
		print (ring_1, '. ring')
		#if atom_list[i][1]: #planar ring
		for ring_2 in range(1,len(atom_list)):
			#print (j, 'j **********')
			#if atom_list[j][1] and i!=j and j not in remove_from_list: #planar ring, don't compare identical rings
			if  ring_1 !=ring_2 and ring_2 not in remove_from_list: #don't compare identical rings
				#check if they have atoms in common => fused
				for atom_idx in atom_list[ring_2][0]:
					#print (atom_idx, 'atom_idx')
					if atom_idx in atom_list[ring_1][0]: #fused
						#fuse both rings
							
						for atom_idex_2 in atom_list[ring_1][0]:
							if atom_idex_2 not in atom_list[ring_2][0]:
								atom_list[ring_2][0].append(atom_idex_2)
						if debug:
							print (atom_list[ring_2][0], "joined list")
						remove_from_list.append(ring_1)
						#print ('remove form list', remove_from_list)
						if debug:
							print('fused rings')
						break #no need to ceck remaining atoms of this list, is already fused and joined
	#the remove list can contain the same ring several times => clean up
	remove_from_list = list(dict.fromkeys(remove_from_list))
	if debug:
		print(remove_from_list, 'remove fused rings') 

	for i in range(0,len(remove_from_list)):
		del atom_list[remove_from_list[i]-i] #index will change
			



	print (atom_list)
	count = len(atom_list)
	return atom_list, count#---------------------

#----------------------------
def super_impose_only_on_planar_rings(ringlist,mol, count):
	#for now, I'll only extend planar ring systems
	pred = OEPartPredAtom(ringlist)
	refmol = OECreateOEMol()
	
	old_list = ringlist
	ring_count = 0

	match_list = []
	for j in old_list:
		match_list.append(j)
	
	for i in xrange(1,count+1):
  		pred.SelectPart(i)
  		OESubsetMol(refmol, mol, pred) 
  		if chemistry.check_planar(refmol):
			#print "keep this ring"
			#print ringlist
  		        ring_count = ring_count + 1
			
  		else:
			#print "reject this ring"
			#delete this ring from ring_list
  		        for atom_number in range(0,len(match_list)):
				#print match_list[atom_number]
  		                if match_list[atom_number] == i:
  		                        match_list[atom_number] = 0
   		                	#print 'test'
   		                #adjust ring_numbers
  		                elif match_list[atom_number] > i:
  		                        match_list[atom_number] = match_list[atom_number] -1
   		                	
			

	#print match_list
	return match_list, ring_count

   		                	
	
#---------------------
def OrderMol2db(mol, atom_list):
#for example: if newOrder is [3,2,0,1], then atom 3 in the original molecule will be atom 0 in the new one
	sorted_atom_list = []
	index_list = []
	#tmp_mol = Chem.Mol()
	for id in atom_list:
		sorted_atom_list.append(id)
		index_list.append(id)
	for atom in mol.GetAtoms():
		idx = atom.GetIdx()
		if idx not in sorted_atom_list:
			sorted_atom_list.append(idx)
			index_list.append(idx)

	tmp_mol = Chem.rdmolops.RenumberAtoms(mol, sorted_atom_list)

	return tmp_mol, index_list

#---------------------
def  get_matches(mol, patlist):
	#generate empty list
	pattern_matches = []

	match = 0
	for pat,smarts,description in patlist:
		for matchbase in pat.Match(mol,True) :
			#print description
			#generate empty list
			pattern_match = []
			for index in range(1, mol.NumAtoms() + 1):
				pattern_match.append(0)

			match = match + 1
			#print match
			for matchpair in matchbase.GetAtoms():
				#print matchpair.target.GetIdx(), matchpair.target.GetName()
				pattern_match[matchpair.target.GetIdx()] = match
			pattern_matches.append(pattern_match)
      


	
	#for i in pattern_matches:
		#print i

	return pattern_matches, match
#---------------------
def superpos(file_name_prefix, mol_list, amsol_dir, only_planar_rings=False, pattern_file =None,debug=False):
	"""superposed conformational ensembles on rings or smarts patterns
	ARGUMENTS:
	      - file_name_prefix: prefix of file name
	      - mol_list: list with ensembles of molecules
	      - amsol_dir: directory under which xx*cav_CYC lies (bad name for variable)
	      - only_planar_rings: superimpose only on planar rings (optional)
	      - zinc_title: (optional) zinc titles have to be treated diffrently, because they are already "pre-adjusted" (they are not unique when you download them)
	      - pattern_file: (optional) file with smarts for pattern matching
	      - debug: (optinal) gives more output

	    RETURNS:

	      None (output will be saved as files)
	"""


#rot=DoubleArray(9)
#trans=DoubleArray(3)


	if pattern_file != None:
		pattern_matching = True
		smarts_file = open(pattern_file, 'r')
		rms_cut_off = 0.001 # fewer atoms in pattern matching -> must be more accurate, otherwise mol2db screws up
		print ("pattern matching")
	else:
		pattern_matching = False
		print ("ring matching")
		rms_cut_off = 0.015 #used to be 0.01 with OE kits. Does 0.025 really work for mol2db?






	solv_dict = {}
	output_dict = {}

	#pattern matching is currently not working

	if pattern_matching:
		#some options from the original match.py file
		usa=1; asym=1; exph=0; max=0; smartsfile=None; verbose=0;
		qmolfile=None; kekule=0; atom_expr=None; bond_expr=None;
		smartslist=None; qmollist=None; nowarn=0;

		smarts_list=chemistry.get_smartslist_file(smarts_file)

		patlist= chemistry.get_patlist(smarts_list, atom_expr, bond_expr,max)



	#out_file = file_name_prefix + '_om_mult_rings.mol2'
	out_file_tmp = file_name_prefix + '_mult_rings_tmp.mol2'
	solv_table_file_name = amsol_dir + '/' + file_name_prefix + '_amsol_cav_CYC/' + file_name_prefix +'_CYC.solv'
	#print solv_table_file_name
	new_solv_table_file_name = amsol_dir + '/' + file_name_prefix + '_amsol_cav_CYC/' + file_name_prefix +'_mult_rings_CYC.solv'
	solv_table_file = open(solv_table_file_name, 'r')
	#new_solv_talbe_file = open(new_solv_table_file_name, 'w')
	#read solvation_table, this might have to be changed, if tables get too long to be held in memory
	solv_dict = read_solv_table(solv_table_file)
	#print solv_dict.keys()

	for mol in mol_list: #loop over all molecule esembles in the list
		title =  mol.GetProp("_Name")
		print (title)
		output_dict[title]=[]
		#go only on, if solvation energy calculation was sucessfull
		if title in solv_dict:
			if not pattern_matching:
				#get the atoms that have to be superimposed
				if only_planar_rings: #currently not working
					match_list, count = super_impose_only_on_planar_rings(ringlist,mol, count)
					#print count
				else:
					match_list, count = fuse_rings(mol,debug)
			else: #currently not working
				#print title
				match_list, count = get_matches(mol, patlist)
				#print match_list
				refmol = OECreateOEMol()
				super_pos_counter = 0
			#print count, ' count'

			cids = [x.GetId() for x in mol.GetConformers()] #get all conformer ids

			#if there is only one conformer, I don't need to align them

			if len(cids) == 1:
				count = 1 #no need to loop over all rings of only one conformer

			super_pos_counter = 0


			for i in range(0,count): #loop over ring systems
				print (match_list, 'match_list')
				print ('ring id:', i)
				saved= 0

#write out amsol charges if this easily possible

				#align molecules on best reference molecule				
				found_ref = [] 
				for conf_id in cids:								
					if conf_id not in found_ref:
						found_ref.append(conf_id)
						align_list = [conf_id] #align onto current conformer
						for cid_2 in cids: #find other conformers to align
							if cid_2 not in found_ref:
								align_list.append(cid_2) #align these conformers in this round

						rms_values = []
						AllChem.AlignMolConformers(mol,  confIds=align_list, atomIds=match_list[i][0], RMSlist=rms_values)
						if debug:						
							print ('cids:', cids, 'len rms_values', len(rms_values), rms_values)
							print ('align_list', align_list)
						found_this_round = [conf_id]
						for i_rms in range (0, len(rms_values)): #find conformes that aligned well
							#first entry in rms_values list corresponds to alignement of 2nd entry in align list on first entry in align list and so on
							if rms_values[i_rms] <= rms_cut_off:
								found_ref.append(align_list[i_rms+1])
								found_this_round.append(align_list[i_rms+1])
						print ('found this round:', found_this_round)
						#write out aligned molecules
						super_pos_counter = super_pos_counter + 1
						print ("super_pos_counter", super_pos_counter)
						#for zinc, you have to adjust the titles before => don't destroy this adjustment here
						new_title = gen_new_title(title, super_pos_counter)
						#mol.SetProp("_Name", new_title)
						print (mol.GetProp("_Name") , "name >----------")
	 					#order atoms for mol2db
						tmp_mol, index_list = OrderMol2db(mol, match_list[i][0])
						tmp_mol.SetProp("_Name", title)

						print (index_list, "index_list")

						#reorder solv table	
						#new_first_line = solv_dict[title][0].replace(title, new_title)
						new_solv_table_file_name = amsol_dir + '/' + file_name_prefix + '_amsol_cav_CYC/' + new_title +'_mult_rings_CYC.solv'
						new_solv_table_file = open(new_solv_table_file_name, 'w')
						new_solv_table_file.write(solv_dict[title][0])
						for index in index_list: #sort solvation table according to atom indexes
							new_solv_table_file.write(solv_dict[title][1][index])
						new_solv_table_file.close()


						#problem with mol2db: aligned is not accurate enough, does not find identical atoms
						#get coordinates from ref conformer
						#use them for matched conformers

						if len(found_this_round) > 1: #no need to resort if only one molecule in the list
							coord_list = []
							for atom_idx in match_list[i][0]:
								pt1 = mol.GetConformer(found_this_round[0]).GetAtomPosition(atom_idx)
								coord_list.append(pt1)

						out_file = new_title + '_om_mult_rings.mol2'
						output_dict[title].append(new_title)

						for cid_3 in found_this_round:
							#saved = saved + 1
							if len(found_this_round) > 1: #update coordinates, no need to do if only one molecule in the list
								for pos in range (0, len(coord_list)): #atoms are already sorted
									tmp_mol.GetConformer(cid_3).SetAtomPosition(pos, coord_list[pos])

							#get amsol charges to write into mol2 file
							for atom in tmp_mol.GetAtoms():
								pos_unsorted_solv_table = index_list[atom.GetIdx()]
								#print (pos_unsorted_solv_table, atom.GetIdx())
								charge = solv_dict[title][1][pos_unsorted_solv_table].split()[0].strip()
								#print (charge)
								atom.SetProp('_GasteigerCharge', charge )

							Mol2Writer.MultiMolToMol2File([tmp_mol], out_file, confId=cid_3, addHs=True, append=True, addCharges=False)

								


				#print ('saved: ', saved, 'ring id:', i)
				#saved = 0

		
				

	#new_solv_talbe_file.close()
	solv_table_file.close()
	print ("done superpose_conf_ensemble!")
	return output_dict


#!/usr/bin/env python3


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
def gen_new_title(title, counter, zinc_label):
	if zinc_label:
		#it is an already adjusted title
		old_number = int(title[4]) * 10 #this assumes that there are a maximum of nine versions of the same molecules in Zinc
		new_number = old_number + counter 
		new_title = 'MFC' + str(new_number).rjust(2) + title[5:]
	else:
		new_title = 'MFC' + title[5:] + str(counter).rjust(2)

	new_title = new_title.replace(' ', '0')
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
def find_atoms_in_plane(mol,debug):  

#only check if ring is planar
#should not be necessary to extend rings as we say exactly which atom matches which one
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
		else:
	      	  atom_list[i].append(False) #ring is not planar



	#print (atom_list)
	count = len(atom_list)
	return atom_list, count#---------------------
def super_impose_only_on_planar_rings(ringlist,mol, count):
	#for now, I'll only extend planar ring systems
	pred = OEPartPredAtom(ringlist)
	refmol = OECreateOEMol()
	
	old_list = ringlist
	ring_count = 0

	new_list = []
	for j in old_list:
		new_list.append(j)
	
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
  		        for atom_number in range(0,len(new_list)):
				#print new_list[atom_number]
  		                if new_list[atom_number] == i:
  		                        new_list[atom_number] = 0
   		                	#print 'test'
   		                #adjust ring_numbers
  		                elif new_list[atom_number] > i:
  		                        new_list[atom_number] = new_list[atom_number] -1
   		                	
			

	#print new_list
	return new_list, ring_count

   		                	
	
#---------------------
def OrderMol2db(mol2, pred):
	index_list = []
	atom_list_match = []
	for atom in mol2.GetAtoms(pred):
	#print atom.GetIdx()
		atom_list_match.append(atom)
		index_list.append(atom.GetIdx())
	for atom in mol2.GetAtoms():
		if atom not in atom_list_match:
			atom_list_match.append(atom)
			index_list.append(atom.GetIdx())

	mol2.OrderAtoms(atom_list_match)
	return mol2, index_list

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
def superpos(file_name_prefix, mol_list, amsol_dir, only_planar_rings=False, zinc_title=False, pattern_file =None,debug=False):
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
		rms_cut_off = 0.01

	#zinc titles have to be treated diffrently, because they are already "pre-adjusted" (they are not unique when you download them
	if zinc_title:
		print ("processing zinc files")
	else:
		print ("processing non-zinc files")





	solv_dict = {}

	tmp_mol_list = []
	tmp2_mol_list = []

	#pattern matching is currently not working

	if pattern_matching:
		#some options from the original match.py file
		usa=1; asym=1; exph=0; max=0; smartsfile=None; verbose=0;
		qmolfile=None; kekule=0; atom_expr=None; bond_expr=None;
		smartslist=None; qmollist=None; nowarn=0;

		smarts_list=chemistry.get_smartslist_file(smarts_file)

		patlist= chemistry.get_patlist(smarts_list, atom_expr, bond_expr,max)



	out_file = file_name_prefix + '_mult_rings.mol2'
	solv_table_file_name = amsol_dir + '/' + file_name_prefix + '_amsol_cav_CYC/' + file_name_prefix +'_CYC.solv'
	#print solv_table_file_name
	new_solv_table_file_name = amsol_dir + '/' + file_name_prefix + '_amsol_cav_CYC/' + file_name_prefix +'_mult_rings_CYC.solv'
	solv_table_file = open(solv_table_file_name, 'r')
	new_solv_talbe_file = open(new_solv_table_file_name, 'w')
	#read solvation_table, this might have to be changed, if tables get too long to be held in memory
	solv_dict = read_solv_table(solv_table_file)
	#print solv_dict.keys()

	for mol in mol_list: #loop over all molecule esembles in the list
		title =  mol.GetProp("_Name")
		print (title)
		#go only on, if solvation energy calculation was sucessfull
		if title in solv_dict:
			if not pattern_matching:
				#get the atoms that have to be superimposed
				if only_planar_rings: #currently not working
					new_list, count = super_impose_only_on_planar_rings(ringlist,mol, count)
					#print count
				else:
					new_list, count = find_atoms_in_plane(mol,debug)
			else: #currently not working
				#print title
				match_list, count = get_matches(mol, patlist)
				#print match_list
				refmol = OECreateOEMol()
				super_pos_counter = 0
			#print count, ' count'
			for i in range(0,count): #loop over ring systmes
				print (new_list, 'new_list')
				ref_mol_list = []
				empty_rms = []

#here, i need to check if I get a good superpostion of all molecules on the current reference, if not add current molecule as reference

				Mol2Writer.MultiMolToMol2File([mol], '_before.mol2', confId=None, addHs=True)

				#weights_list1 = []
				#weights_list2 = []
				#for r in new_list[i][0]:
					#weights_list1.append(0.10)	
					#weights_list2.append(10)	

				#I think the following aligns 1,3 to 1 and so on!
				#print (AllChem.AlignMolConformers(mol, atomIds=new_list[i][0], RMSlist=empty_rms))
				#print (empty_rms)
				#calc_rmsd (mol, new_list[i][0],[0,1,2,3,4,5,6,7,8,9])
				#Mol2Writer.MultiMolToMol2File([mol], str(i) + '_om1_weights1.mol2', confId=None, addHs=True)

				#print (AllChem.AlignMolConformers(mol, confIds=[0, 1, 2],atomIds=new_list[i][0], RMSlist=empty_rms,maxIters=5000,weights=weights_list2))
				#print (empty_rms)

				#Mol2Writer.MultiMolToMol2File([mol], '_om1_weights2.mol2', confId=None, addHs=True)



				confs = mol.GetConformers()
				for aid in new_list[i][0]:
			        	mpos = 0
			        	print ('new atom')
			        	for k, conf in enumerate(confs):
			        		if (k == 0):
			        			mpos = list(conf.GetAtomPosition(aid))
			        			continue
			        		else:
			        			pos = list(conf.GetAtomPosition(aid))
				
			        			#print ('a', mpos)
			        			#print('b', pos)
			        			#lstFeq(mpos, pos, .5)

				Mol2Writer.MultiMolToMol2File([mol], '_om1.mol2', confId=None, addHs=True)
				#print(Chem.MolToMolBlock(mol),file=open('foo.mol','w+'))

				print (new_list[i][0], 'atoms for matching')
				weights_list = []
				for r in new_list[i][0]:
					weights_list.append(10)						

				ms = []
				cids = [x.GetId() for x in mol.GetConformers()]
				for cid in cids:
					newmol = Chem.Mol(mol)
					for ocid in cids:
						if ocid == cid:
							continue
						newmol.RemoveConformer(ocid)
					ms.append(newmol)
				mol2_text = ''
				for new_mol in ms:
					mol2_text = mol2_text + Mol2Writer.MolToMol2Block(new_mol,addCharges=False)
				out_file = open('f_om2_new.mol2', 'w')
				out_file.write(mol2_text)
				out_file.close()




				print (AllChem.AlignMolConformers(mol, confIds=[3, 2, 3],atomIds=new_list[i][0],RMSlist=empty_rms,maxIters=5000))
				print (empty_rms)
				Mol2Writer.MultiMolToMol2File([mol], '_om2.mol2', confId=None, addHs=True)
				#print(Chem.MolToMolBlock(mol),file=open('foo2.mol','w+'))

					
				#ref_mol_list.append(refmol)


#to do: need to attach conformer to ref_mol_list, not sure if conf has to be converted to mol, or maybe I just have to save the conformer ID?
#I think confId should work
#go on after if not planar below
#HTtp://rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html
#would about fused rings? They end up to be two rings. Would be good if we could just save the aromatic part in case mixed fused? They are easy to find: atoms do occur in two lits => can just be combined (and if only armatci onces should be kept, that can be done in a similar way by first testing if the rings a planar)
#I think the alignment works, but the rmsd values are wrong, the change also dependent on the weight => have added a function
#already with some planar rings I can not get very low rmsds => do they also need to be seperated to not confuse mol2db? => test

				 #find out, how many reference molecules are needed (because some rings are not planar)
				#this has also to be done for pattern matching because not all patterns are planar					
				if not pattern_matching:
					print (new_list, new_list[i][1])
					planar = new_list[i][1]
				else:
					planar = False # always check if pattern matching
				print (planar, '**********')
				if not planar:
						   #print mol.GetTitle()
						   #while OEReadMolecule(ifs2, mol2):
						   for mol2 in tmp_mol_list:
	   					        thismol = OECreateOEMol()
	   					        OESubsetMol (thismol, mol2, pred)
							#print 'Check reference molecules'
	   					        found = False
	   					        for refmol in ref_mol_list:
	   					        	if not found:
									#rms=OERMSD(refmol,thismol,0,1,1,rot,trans)
	   					        		rms=OERMSD(refmol,thismol,False,True,True,rot,trans)

									#print 'rmsd %5.2f\n' % rms
	   					        		if rms <= rms_cut_off:
	   					        			found = True
										#print found
							#further reference molecules needed?
	   					        if not found:
	   					        	refmol = OECreateOEMol()
	   					        	OESubsetMol(refmol, mol2, pred)
	   					        	ref_mol_list.append(refmol)

				for refmol in ref_mol_list:
						   super_pos_counter = super_pos_counter + 1
						   #print super_pos_counter
						  
						   #ifs2.close()
						   #ifs2 = oemolistream('tmp.mol2')
						   #while OEReadMolecule(ifs2, mol2):
						   for mol2 in tmp_mol_list:
							#print 'here', mol2.GetTitle()
	   					        thismol = OECreateOEMol()
	   					        OESubsetMol (thismol, mol2, pred)
	   					        smi = OECreateSmiString(thismol)
	   					        #print smi
	   					        #rms=OERMSD(refmol,thismol,0,1,1,rot,trans)
	   					        rms=OERMSD(refmol,thismol,False,True,True,rot,trans)
	   					        #print rms
	   					        if (rms <= rms_cut_off): 
								#print 'rmsd %5.2f\n' % rms
	   					                OERotate(mol2,rot)
	   					                OETranslate(mol2,trans)

								#for zinc, you have to adjust the titles before => don't destroy this adjustment here
	   					                new_title = gen_new_title(title, super_pos_counter, zinc_title)

	   					                mol2.SetTitle(new_title)

								#order atoms for mol2db
	   					                mol2, index_list = OrderMol2db(mol2, pred)					        
								
	   					                OEWriteMolecule(ofs_mult, mol2)
								#print mol2.GetTitle()
						  #identifier must be unique
						   #new_solv_talbe_file.write('MFC' + solv_dict[title][0][4:12] + str(super_pos_counter) + solv_dict[title][0][12:])	
						   #print index_list
						   new_solv_talbe_file.write(new_title + solv_dict[title][0][12:])
						   for index in index_list:
						   	new_solv_talbe_file.write(solv_dict[title][1][index])
				



			if not ifs_mult.IsValid():
					go_on = False
					#Here I would lose the last molecule if it is rigid

					#ifs = oemolistream('tmp2.mol2')
					last_mol = OEGraphMol()
					#OEReadMolecule(ifs, last_mol)
					#if the last molecule is an ensemble => tmp2_list is empty
					if len(tmp2_mol_list) > 0:
						last_mol = tmp2_mol_list[0]
						title = last_mol.GetTitle()
						new_title = gen_new_title(title, 1, zinc_title)
						last_mol.SetTitle(new_title)	

						#order atoms for mol2db
						last_mol, index_list = OrderMol2db(last_mol, pred)				
						OEWriteMolecule(ofs_mult, last_mol)
						new_solv_talbe_file.write(new_title + solv_dict[title][0][12:])	
						for index in index_list:
					 		new_solv_talbe_file.write(solv_dict[title][1][index])

				
				

	new_solv_talbe_file.close()
	solv_table_file.close()
	command ='nice bzip2 ' + file_name
	#os.system(command)
	command = 'nice bzip2 ' +  file_name[:-5] + '_mult_rings.mol2'
	#os.system(command)


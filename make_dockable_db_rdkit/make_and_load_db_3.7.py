#!/usr/bin/env python3

#to do: errors are not caught properly if mol2db failes, can be tested if cmpd id does not start with MFC
#mol2 for error catching is not read properly, could simply treat it as text file, blocks are read correctly=> use them, or use PyBel  (python wrapper for OpenBable)?

#it's better to call this script from cl_make_and_load_db.py because then a check is performed if all necessary files are present and some directories are created

#in a few tests, this seems to fail regualary with molecules with more than 7 rotatable bonds, can this be fixed? (fails at mol2db stage)


#use $TMPDIR on cluster

import sys, os, os.path,time,subprocess, signal, my_mysql3 as mysql
import chemistry3 as chemistry
import Mol2Writer
import superpos_conf_ensemble_3_7 as  superpos_conf_ensemble
import amsol_functions

from rdkit import Chem
from rdkit.Chem import AllChem

from argparse import ArgumentParser


description = "Script to make dockable db and load it into MySQL db"

parser = ArgumentParser(description=description)

parser.add_argument("-u", "--username", type=str, action="store", dest="username", help="username for Mysql DB",required=False)
parser.add_argument("-p", "--password", type=str, action="store", dest="password",  help="password for Mysql DB",required=False)
parser.add_argument("-i", "--input_file", type=str, action="store", dest="in_file", help="input file (smiles)",required=True)
parser.add_argument("-c", "--cluster",  action="store_true", dest="cluster", default=False, help="running on cluster?",required=False)
parser.add_argument("-up", "--updatedb",  action="store_true", dest="update_db", default=False, help="Update MySQL DB, default none",required=False)
parser.add_argument("-uf", "--update_failed",  action="store_true", dest="update_failed", default=False, help="Record failed molecues in MySQL DB, default none",required=False)
parser.add_argument("-pro", "--project", type=str, action="store", dest="project",  help="name of project in Mysql DB",required=False)
parser.add_argument("-s", "--script_path", type=str, action="store", dest="script_path",  help="path to required scripts",required=True)
parser.add_argument("-smarts", "--smarts_file", type=str, action="store", dest="smarts_file",  default = None, help="file with smarts for smarts machting",required=False)
parser.add_argument("-d", "--debug",  action="store_true", dest="debug", default=False, help="Debug, default none",required=False)
parser.add_argument("-r", "--rmsd",  type=float, action="store", dest="rmsd", default=1.0, help="rmsd cut off for rotamers, default 1.0",required=False)



args = parser.parse_args()
args_dict = vars(args)

#print (args_dict)
locals().update(args_dict) #generates local variables from key / value pairs

#print (smarts_file)



if update_db or update_failed:
	if (not username) or (not password) or (not project):
		print ('need username, password and project in order to update DB')
		sys.exit(1)

if script_path[-1] != '/':
	script_path = script_path + '/'	




#----------------------------
def reformat_mol2_file(format_file):
	#OElib now produces a format that is not compatible with mol2db -> change
	correct_file = open(format_file,'r')
	correct_file_text = correct_file.read()
        #print correct_file_text
	correct_file_text = correct_file_text.replace('NO_CHARGES', 'NO_CHARGES\n')
	correct_file.close()
	correct_file = open(format_file,'w')
	correct_file.write(correct_file_text)
        #print correct_file_text
	correct_file.close()
#----------------------------
def solv_remove(remove_id, solv_file):
	command = 'rm temp.solv'
	os.system(command)
	temp_solv = open('temp.solv', 'w')
	print (solv_file)
	org_solv_file = open(solv_file, 'r')
	for line in org_solv_file.readlines():
		#print line
                if line[0:3] == 'MFC':
                	title = line[0:12]
                	#print '----------- check ', title
                	if title in removed_id:
                		delete = True
                	else:
                        	delete = False
                if not delete:
                	temp_solv.write(line)
	temp_solv.close()
	org_solv_file.close()
#----------------------------


#check if input file exists -> is an issue with array jobs because they always expect consecetutive numbering which might not always be the case

if not os.path.exists(in_file):
	print (in_file, 'does not exist')
	sys.exit()

if update_db or update_failed:

	conn=mysql.connect2server(password, username,'purchasable')  
	cursor = conn.cursor ()
			


	#get table names

	command = "select stereoisomer_smiles from purchasable.projects where project = '" + project + "'"
	print (command)
	cursor.execute(command)
	results=cursor.fetchall()
	if len(results) == 0:
		print ('project does not exist')
		sys.exit(1)

	stereo_table = results[0][0]

	#get stereo_id,blob, error from table
	command = "show fields from " + stereo_table
	cursor.execute(command)
	results=cursor.fetchall()
	for result in results:
		print (result[3], result[0])
		if result[3] == 'PRI':
			stereo_field = result[0]
		elif result[0].find('blob') != -1 :
			blob_field = result[0]
		elif result[0].find('error') != -1 :
			error_field = result[0]
	print (stereo_field, blob_field, error_field)
	#error_field = blob_field


cur_dir = os.getcwd()


#make some needed directoryies
if not os.path.exists('omega_mult'):
	os.makedirs('omega_mult')
if not os.path.exists('solv'):
	os.makedirs('solv')
if not os.path.exists('errors'):
	os.makedirs('errors')	
if not os.path.exists('input_smiles'):
	os.makedirs('input_smiles')	

#read smiles file
ifs = Chem.SmilesMolSupplier(in_file,delimiter='\t',titleLine=False)
for mol in ifs:

	name = mol.GetProp("_Name")


	#make subdir
	dir_name = name + '_work_dir'
	print (dir_name)

	if os.path.exists(dir_name):
		print ('working dir already exists:', in_file)
	else:
		os.makedirs(dir_name)
		os.chdir(dir_name)
		#os.makedirs('amsol2')
		os.makedirs('omega_mult')
		if debug:
			smi_file = open(name + '.smi', 'w')
			smi_file.write(Chem.MolToSmiles(mol) + '\t' + name)
			smi_file.close()
 


		#smarts file
		if smarts_file != None:
			command = 'cp ' + smarts_file + ' omega_mult'
			print (command)
			os.system(command)
		#command = 'cp  ' + cur_dir + '/omega.* .'
		#os.system(command)

		#amsol has a problem if the frist for atoms contain a nitril group, same for C#C, 
		patt = Chem.MolFromSmarts('*#*')
		hit_ats = list(mol.GetSubstructMatch(patt))
		print (hit_ats)
		at_id_counter = 0
		for at_id in hit_ats:
			if at_id <=4:
				at_id_counter = at_id_counter + 1
		if debug:
			print (at_id_counter, ", if > 2 => nitril group among first four atoms")
		if at_id_counter >=2: #there is a nitrile group among the first 
			#at this stage, there a no H atoms in the molecule
			new_order = []
			for atom in mol.GetAtoms():
				if atom.GetIdx() not in hit_ats:
					new_order.append(atom.GetIdx())
			for at_id in hit_ats:
				new_order.append(at_id)
			mol = Chem.RenumberAtoms(mol, new_order)
			mol.SetProp("_Name",name)


		

	
		#omega
		single_conf_mol = chemistry.gen_confs_single(mol,  1,debug=debug,rmsd_threshold=rmsd)
		single_conf_mol_list=[single_conf_mol]
		if debug: #save mol2 file
			Mol2Writer.MultiMolToMol2File(single_conf_mol_list, './' + name + '_os.mol2', confId=None, addHs=True)
		#print (single_conf)


		#amsol
		failed_amsol_list = amsol_functions.calc_amsol(single_conf_mol_list, 'CYC',name+'.smi',debug=debug) 


		#mark failed molecules, amsol
		for name in failed_amsol_list:
			if update_db or update_failed:
				name = name.replace('MFC', '')
				command = 'update ' + stereo_table + ' set ' + error_field + " = 'amsol' where " + stereo_field + ' =' + name
				print (command)
				cursor.execute(command)
				conn.commit()
			else:
				print (name, 'AMSOL failed')
		if len(failed_amsol_list) > 0:
			#molecule failed, we only treat one molecule at a time, no need to go on
			os.chdir('../')
			continue


	

		#omega mult
		print ('calc conformations')
		mult_conf_mol = chemistry.gen_confs_single(mol, 1000,debug=debug,rmsd_threshold=rmsd)
		mult_conf_mol_list = [mult_conf_mol]
		if debug: #save mol2 file
			Mol2Writer.MultiMolToMol2File(mult_conf_mol_list, 'omega_mult/' + name + '_om.mol2', confId=None, addHs=True)
		




	

	
		#multio....
		os.chdir('omega_mult')
		output_dict = superpos_conf_ensemble.superpos(name,mult_conf_mol_list,'..', smarts_file,debug=debug)
	
		for title in output_dict.keys():
			#create name file
			name_file = open('name.txt','w')
			name_line = 'name.txt 1 ' + title + ' ../' + in_file + ' | NO_LONG_NAME\n'
			name_file.write(name_line)
			name_file.close()

			for new_mol in output_dict[title]:
				command = '/scratch/software/anaconda3/envs/py2/bin/python /scratch/software/DOCK-3.7.2rc1/ligand/mol2db2/mol2db2.py -m ' + new_mol + '_om_mult_rings.mol2 -s ../' + name + '_amsol_cav_CYC/' + new_mol + '_mult_rings_CYC.solv -o ' + new_mol + '.db.gz'
				print (command)
				os.system(command)

		os.chdir('../..')
		print (os.getcwd())

			


print ('make_and_load_dockabe.py done')

			

		

		
		
		
	
	
	

	
	
	


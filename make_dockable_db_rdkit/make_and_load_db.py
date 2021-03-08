#!/usr/bin/env python3

#to do: errors are not caught properly if mol2db failes, can be tested if cmpd id does not start with MFC
#mol2 for error catching is not read properly, could simply treat it as text file, blocks are read correctly=> use them, or use PyBel  (python wrapper for OpenBable)?

#it's better to call this script from cl_make_and_load_db.py because then a check is performed if all necessary files are present and some directories are created

#in a few tests, this seems to fail regualary with molecules with more than 7 rotatable bonds, can this be fixed? (fails at mol2db stage)


#use $TMPDIR on cluster

import sys, os, os.path,time,subprocess, signal, my_mysql3 as mysql
import chemistry3 as chemistry
import Mol2Writer
import superpos_conf_ensemble
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


#make subdir
dir_name = in_file[:-4] + '_work_dir'
print (dir_name)

if os.path.exists(dir_name):
	print ('working dir already exists:', in_file)
else:
	os.makedirs(dir_name)
	os.chdir(dir_name)
	#os.makedirs('amsol2')
	os.makedirs('omega_mult')
	if cluster:
		command = 'cp ' + cur_dir + '/' + in_file + ' .'
	else:
		command = 'mv ../' + in_file + ' .'
	print (command)
	os.system(command)
	#smarts file
	if smarts_file != None:
		command = 'cp ' + smarts_file + ' omega_mult'
		print (command)
		os.system(command)
	#command = 'cp  ' + cur_dir + '/omega.* .'
	#os.system(command)
	
	#omega
	single_conf_mol_list = chemistry.gen_confs(in_file[:-4] + '.smi', 1,debug=debug)
	if debug: #save mol2 file
		Mol2Writer.MultiMolToMol2File(single_conf_mol_list, './' + in_file[:-4] + '_os.mol2', confId=None, addHs=True)
	#print (single_conf)

	#command = script_path + 'omega_script.py ' + in_file[:-4] + '.smi' + in_file[:-4] + '_os.mol2 1 single.log'

	#print (command)
	#os.system(command)
	

	#omega mult
	print ('calc conformations')
	mult_conf_mol_list = chemistry.gen_confs(in_file[:-4] + '.smi', 1000,debug=debug)
	if debug: #save mol2 file
		Mol2Writer.MultiMolToMol2File(mult_conf_mol_list, 'omega_mult/' + in_file[:-4] + '_om.mol2', confId=None, addHs=True)
		

	#amsol
	failed_amsol_list = amsol_functions.calc_amsol(single_conf_mol_list, 'CYC',in_file,debug=debug) 


	
	#mark failed molecules
	for name in failed_amsol_list:
				if update_db or update_failed:
					name = name.replace('MFC', '')
					command = 'update ' + stereo_table + ' set ' + error_field + " = 'amsol' where " + stereo_field + ' =' + name
					print (command)
					cursor.execute(command)
				else:
					print (name, 'failed')

	if update_db or update_failed:
		conn.commit()
		#cursor.close ()
		#conn.close
	
	#multio....
	os.chdir('omega_mult')
	superpos_conf_ensemble.superpos(in_file[:-4],mult_conf_mol_list,'..', smarts_file,debug=debug)

	#reformat_mol2_file(in_file[:-4]+'_om_mult_rings.mol2')

	#mol2db
	os.system('cp ' + cur_dir + '/inhier .')
	omega_file = in_file[:-4] + '_om_mult_rings.mol2'
	db_file = omega_file[:-5] + '.db'
	command = 'ln -s ' + omega_file + ' temp.mol2'
	print (command)
	os.system(command)
	solv_file = '../' + in_file[:-4] + '_amsol_cav_CYC/' + in_file[:-4] + '_mult_rings_CYC.solv'
	command = 'ln -s ' + solv_file +  ' temp.solv'
	print (command)
	os.system(command)

	moldb_command = ['mol2db', 'inhier']

	
	#sometimes mol2db hangs -> catch this
	mol2db_loop_finished = False #rerun mol2db until it finishes normally or no molecules are left
	removed_id = []
	attempts = 0
	while (not mol2db_loop_finished): 
		attempts = attempts + 1
		mol2db_log = in_file[:-4] + '_om_mult_rings.log'
		log_file = open(mol2db_log, 'w')
		print (moldb_command)
		proc = subprocess.Popen(moldb_command,stdout=log_file) #start process

		#let it try for ten minutes
		minutes = 1
		error = False
		mol2db_finished = False
		while (not mol2db_finished) and (minutes <= 6):
			time.sleep(10) #short sleep, enough for short jobs
			print (proc.poll(), 'proc.poll()')
			if proc.poll() == 0: #job finished normally
				print ('finished normally')
				mol2db_finished = True
				error = False
			elif proc.poll() is not None: #The jobs has finished, but there was an error
				#do something
				#print proc.poll()
				print ('there was a problem with mol2db')
				mol2db_finished = True
				if proc.poll() == -11 and attempts > 1:
					error = True # this used to be False because "there is a segmentation fault because no molecules are left => quit" but this is not correct
				else:
					error = True

			elif proc.poll() is None: #job is still running 
				time.sleep(30) 
				minutes = minutes + 1
				print ('mol2db is still running')
		if not mol2db_finished:
		     #kill job	
		     print ('mol2db is hanging -> will now be killed')	
		     print (proc.pid)
		     print (signal.SIGUSR1)
		     os.kill(proc.pid, signal.SIGUSR1)
		log_file.close()
		if error or not mol2db_finished: # mol2db did not run on all molecules => check log file
			log_file = open( mol2db_log , 'r')
			log_file_lines = log_file.readlines()
			log_file.close()
			#find last line which contains id number
			worked_id = ''
			for i in reversed(log_file_lines):
				#print i
				if len(i) == 74: #line with molecule id
					worked_id = 'MFC' + i.split()[0]
					#print worked_id
					break # found last working id

			if worked_id == '' and attempts > 0:
				#the first molecule caused the problem
				org_solv_file = open(solv_file, 'r')
				first_line = org_solv_file.readline()
				org_solv_file.close()
				#print (first_line)
				failed = first_line[0:12]
				#print (failed, 'failed id')
				removed_id.append(failed)

			#clean mol2 file
			command = 'rm temp.mol2'
			os.system(command)
			#print omega_file
			print ('make sure that this rescue procedure does not change the order of the atoms in the molecules!')
			saved_mols = []
			omega_mols = chemistry.Mol2Supplier(omega_file,sanitize=False) #if not sanitize => can not read all molecules
			delete_following = False
			for mol in omega_mols:

				#print (list( mol.GetPropNames()))

				title = mol.GetProp("_Name") #this does not work, why? If I read the c code correctly, this should be set
				#print (title)
				if title == worked_id:
					delete_following = True
				if title in removed_id:
					#skip
					continue #we do not want these molecules
				elif delete_following and title != worked_id:
					delete_following = False #found the id number
					removed_id.append(title)
				else: 
					saved_mols.append(mol)
					#print title
			Mol2Writer.MultiMolToMol2File(saved_mols, 'temp.mol2', confId=None, addHs=True)


			print (removed_id)
			#reformat_mol2_file('temp.mol2')

			#adopt solv file
			solv_remove(removed_id, solv_file)
		if (not mol2db_finished) or error:
			mol2db_loop_finished = False # repeat mol2db after molis have been removed
		else:
			mol2db_loop_finished = True	
		
	print ('failed in mol2db', removed_id)
	#update db
	if update_db or update_failed:
		for stereo_id in removed_id:
			stereo_id = stereo_id.replace('MFC','')
			stereo_id = stereo_id[0:-2] #last two digits are for ring systems
			command = 'update ' + stereo_table + ' set ' + error_field + " = 'mol2db' where " + stereo_field + ' =' + stereo_id 
			print (command)
			cursor.execute(command)	
			conn.commit()		



	command = 'mv temp.db ' + in_file[:-4] + '_om_mult_rings.db'
	os.system(command)

	#check that log file is okay
	log_file = open( in_file[:-4] + '_om_mult_rings.log' , 'r')
	log_file_text = log_file.readlines()
	log_file.close()
	b_fin = True
	b_err = False
	l_err = []
	for s_line in log_file_text:
		if 'Finished' in s_line:
			b_fin = False
		if 'error' in s_line and not 'Exiting' in s_line:
			b_err = True
			l_err.append(s_line.split()[-1])
	if b_err:
		b_fin = False		
		for s_name in l_err:
			if update_db or update_failed: #for sure only do something if update  is set
				command = 'update ' + stereo_table + ' set ' + error_field + " = 'failed' where " + stereo_field + ' =' + s_name[:-2]
				print (command)
				cursor.execute(command)
				conn.commit()		
			else:
				print ('MFC' + s_name,  'failed')
				
	if b_fin:		
		#did not work
		command = 'mv ../' + in_file  + ' ' + cur_dir + '/errors'
		print (command)
		os.system(command)
		
		
	else:
		#everything went fine
		if update_db:
			cursor.close ()  #close connection before opened new one in update_db section
			conn.close
		
		
		if debug:
			command = 'mv *mult*mol2* ' + cur_dir + '/omega_mult'
			print (command)
			os.system(command)
			command = 'mv ../*CYC/*mult*solv ' + cur_dir + '/solv'
			print (command)
			os.system(command)
		if update_db:
	
			#load file
			#header
			#command ='cat ../../header.txt *db > ../new.db'
			command = 'mv *db ../new.db'
			print (command)
			os.system(command)
			os.chdir('..')

			command = script_path + 'mysql_load_dockable.py ' + username + ' ' + password +  ' ' + stereo_table + ' ' + stereo_field + ' ' + blob_field +  ' ' + 'False'
			print (command)
			upload_success =  os.system(command)                   
			if upload_success == 0:
				print ('upload succesful')
				if not debug:
					if cluster:
						command = 'cp ' + cur_dir + '/' + in_file + ' ' + cur_dir + '/input_smiles'
					else:
						command = 'mv *smi ../input_smiles'
					print (command)
					os.system(command)
					os.chdir('..')
					if not cluster:					
						command = 'rm -rf ' + dir_name
						print (command)
						os.system(command)

			else:
				print ('could not upload molecules')
		elif not debug:
			command = 'mv *db ' + cur_dir + '/new_dbs/'
			print (command)
			os.system(command)
			command = 'bzip2 ' + cur_dir + '/new_dbs/' +  in_file[:-4] + '_om_mult_rings.db'
			os.system(command)
			os.chdir('..')
			if not cluster:
				command = 'mv *smi ../input_smiles'
			else:
				command = 'mv ' + cur_dir + '/' + in_file + ' ' + cur_dir + '/input_smiles'
			print (command)
			os.system(command)
			os.chdir('..')
			if not cluster:
				print (command)
				os.system(command)

print ('make_and_load_dockabe.py done')

			

		

		
		
		
	
	
	

	
	
	


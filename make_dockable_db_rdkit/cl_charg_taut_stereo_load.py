#!/usr/bin/env python

import os, sys
import my_mysql3 as mysql
from argparse import ArgumentParser


#-----------

#*******************************************************************************
# MAIN help
#**********
#


description = "Script to charge and tautomerzie molecules on cluster and load them to db into charged or tautomerzied smiles tables"
#usage = "drugpred.py [options]"

parser = ArgumentParser(description=description)

parser.add_argument("-d", "--directory", type=str, action="store", dest="path", help="directory with charger.py and tautomerizer.py (absolute path)",required=True)
parser.add_argument("-u", "--user",  type=str, action="store", dest="username", default ='webuser', help="username",required=False)
parser.add_argument("-p", "--password",  type=str, action="store", dest="password", help="password", default = '', required=False)
parser.add_argument("-n", "--project_name",  type=str, action="store", dest="project", help="name of project in db",required=False)
parser.add_argument("-t", "--initial",  type=str, action="store", dest="initial", help="initial",required=False)
parser.add_argument("-l", "--load_to_db",  action="store_true", dest="load", default=False, help="load to db when done (can be slow on cluster)",required=False)
parser.add_argument("-j", "--job_to_do",  action="store", dest="job", default='both', help="run charger or tautomerizer",choices=['charger', 'tautomerizer','stereoisomerizer'],required=True)


args = parser.parse_args()
args_dict = vars(args)

#print (args_dict)
locals().update(args_dict) #generates local variables from key / value pairs

#check that everything is entered
failed = False
if load:
	if username == 'webuser' or password == '':
		print ('need username and password to load compounds')
		failed = True
	if not project:
		print ('need project to upload compounds')
		failed = True
	if not initial:
		print ('need initial to upload compounds')
		failed = True

	if failed:
		sys.exit()

#check db

conn=mysql.connect2server(password, username,'purchasable')  
cursor = conn.cursor ()


if load:
	#check project
	command = 'select unique_compounds, supplier_info_table, supplier_table from projects where project = "' + project + '"'
	#print (command)
	cursor.execute(command)
	results=cursor.fetchall()
	if len(results) == 0:
		print ('project not valid')
		conn.close()
		sys.exit(1)
conn.close()

path_cpd_upload = path.replace('make_dockable_db_rdkit', 'compound_upload_rdkit', )

#********************

files = os.listdir('./')

clean_list = []
max_file_number = 0

#remove all none .smi from files
for i in files:
	if (i[-4:] == '.smi'):
		clean_list.append(i)

files = clean_list



files.sort()



first_bit = "#$ -S /bin/tcsh\n#$ -cwd\nworkdir=/$SCRATCH/$SLURM_ARRAY_TASK_ID\nmkdir -p $workdir\ncd $SLURM_SUBMIT_DIR\n"
counter = 0

for file in files:


		counter = counter + 1

		#print 'submit job on cluster'

		file_name = str(counter) + '_start.bin'
		start_file = open(file_name, 'w')
		start_file.write(first_bit)

		start_file.write('echo $workdir\n')

		start_file.write('cp ' + file + ' $workdir\n')

		start_file.write('cd $workdir\n')

		if job == 'charger':

			start_file.write(path_cpd_upload + '/charger.py  ')
			start_file.write('-i ' + file + ' -o ' + file [:-4] + '_charged.smi -p 7 -r 2\n')
			start_file.write('cp -R ' +  file[:-4] + '_charged.smi'  + ' $SLURM_SUBMIT_DIR\n')

		elif job == 'stereoisomerizer':

			start_file.write(path + '/stereoisomerizer.py  ')
			start_file.write('-i ' + file + ' -o ' + file [:-4] + '_stereo.smi\n')
			start_file.write('cp -R ' +  file[:-4] + '_stereo.smi'  + ' $SLURM_SUBMIT_DIR\n')
		else:
			start_file.write(path_cpd_upload + '/tautomerizer.py  ')
			start_file.write('-i ' + file + ' -o ' + file[:-4] + '_taut.smi\n')
			start_file.write('cp -R ' + file[:-4] + '_taut.smi'  + ' $SLURM_SUBMIT_DIR\n')

	
		if load:
			if job == 'charger':
				start_file.write(path + '/mysql_load_prot_states.py ')
				start_file.write(username + ' ' + password + ' ' +  file[:-4] + '_charged.smi ' + project + ' ' + initial + '\n')


			elif job == 'stereoisomerizer':
				start_file.write(path + '/mysql_load_stereoisomers.py ')
				start_file.write(username + ' ' + password + ' ' +  file[:-4] + '_stereo.smi ' + project + ' ' + initial + '\n')
			else:
				start_file.write(path + '/mysql_load_taut_states.py ')
				start_file.write(username + ' ' + password + ' ' +  file[:-4] + '_taut.smi ' + project + ' ' + initial + '\n')



		start_file.close()

		os.system('chmod 744 '  + file_name)


file_name = 'charge_taut_load_array_start.bin'
start_file = open(file_name, 'w')
start_file.write("#!/bin/bash\n")
start_file.write('#SBATCH --time=0-03:00:00\n')
start_file.write('#SBATCH --partition normal\n')
start_file.write('${SLURM_ARRAY_TASK_ID}_start.bin\n')
start_file.close()
os.system('chmod 744 '  + file_name)

command = 'sbatch --array=1-' + str(counter) + ' ' + file_name  
print (command)
os.system(command)



print ('done!!!!!!!!!!!!!')





	
	

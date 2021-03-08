#!/usr/bin/env python

import os, sys
import my_mysql3 as mysql
from argparse import ArgumentParser

try:
	path = os.environ['MakeDble']
except:
	print ("Environment variable MakeDble not set")
	sys.exit()


#-----------

#*******************************************************************************
# MAIN help
#**********
#


description = "Script to generate 3D conformations on cluster"
#usage = "drugpred.py [options]"

parser = ArgumentParser(description=description)

parser.add_argument("-r", "--rmsd_threshold",  type=float, action="store", dest="rmsd_threshold", default=1.0, help="rmsd threshold, default 1.0",required=False)
parser.add_argument("-e", "--energy_diff",  type=float, action="store", dest="energy_diff", default=10.0, help="max energy difference between lowest and highest energy for conformers",required=False)
parser.add_argument("-m", "--max_confs", type=int, action="store", dest="max_confs", default=1, help="maximum number of conformations to generate",required=False)
parser.add_argument("-d", "--debug",  action="store_true", dest="debug", default=False, help="Debug, default none",required=False)
parser.add_argument("-t", "--time", type=int, action="store", dest="time", default=3, help="max time (in hours) per node",required=False)





args = parser.parse_args()
args_dict = vars(args)

locals().update(args_dict) #generates local variables from key / value pairs


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

		start_file.write(path + '/gen_conf.py ' + '-i ' + file + ' -o ' + file[:-4] + '.sdf' + ' -r ' + str(rmsd_threshold) + ' -e ' + str(energy_diff) + ' -m ' + str(max_confs))
		if debug:
			start_file.write(' d')
		start_file.write('\n')		

		start_file.write('cp -R ' +  file[:-4] + '.sdf'  + ' $SLURM_SUBMIT_DIR\n')

		start_file.close()

		os.system('chmod 744 '  + file_name)


file_name = 'conf_gen_start.bin'
start_file = open(file_name, 'w')
start_file.write("#!/bin/bash\n")
start_file.write('#SBATCH --time=0-' + str(time) + ':00:00\n')
start_file.write('#SBATCH --partition normal\n')
start_file.write('${SLURM_ARRAY_TASK_ID}_start.bin\n')
start_file.close()
os.system('chmod 744 '  + file_name)

command = 'sbatch --array=1-' + str(counter) + ' ' + file_name  
print (command)
os.system(command)



print ('done!!!!!!!!!!!!!')





	
	

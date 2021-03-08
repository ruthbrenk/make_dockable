#!/usr/bin/env python3



import sys, chemistry3 as chemistry
from rdkit import Chem
from rdkit.Chem import AllChem
from argparse import ArgumentParser

#rmsd_threshold = 1 #this is what we used in omega.mult
#energy_diff = 10.0


description = "Script to generate conformations"

parser = ArgumentParser(description=description)

parser.add_argument("-i", "--input_smiles", type=str, action="store", dest="smi_input_file", help="smiles file with input molecules",required=True)
parser.add_argument("-o", "--output_sdf", type=str, action="store", dest="sdf_output_file",  help="output SDF",required=True)
parser.add_argument("-r", "--rmsd_threshold",  type=float, action="store", dest="rmsd_threshold", default=1.0, help="rmsd threshold, default 1.0",required=False)
parser.add_argument("-e", "--energy_diff",  type=float, action="store", dest="energy_diff", default=10.0, help="max energy difference between lowest and highest energy for conformers",required=False)
parser.add_argument("-m", "--max_confs", type=int, action="store", dest="max_confs", default=1, help="maximum number of conformations to generate",required=False)
parser.add_argument("-d", "--debug",  action="store_true", dest="debug", default=False, help="Debug, default none",required=False)


args = parser.parse_args()
args_dict = vars(args)

#print (args_dict)
locals().update(args_dict) #generates local variables from key / value pairs


junk,file_type = smi_input_file.split('.')
print ("type", file_type)
if file_type == 'smi':
	ifs = Chem.SmilesMolSupplier(smi_input_file,delimiter='\t',titleLine=False)
elif file_type == 'sdf':
	ifs = Chem.SDMolSupplier(smi_input_file)
	print ('reading SDF')


elif file_type == 'mol2': #this does not seem to be a well supported format

	ifs = chemistry.Mol2Supplier(smi_input_file)

else:
	print ('unknown file type')
	sys.exit()

#input_smiles = Chem.SmilesMolSupplier(smi_input_file, delimiter='\t', titleLine=False)
writer = Chem.SDWriter(sdf_output_file)


for mol in ifs:
	print (mol)
	mol3D = chemistry.gen_confs_single(mol, max_confs,debug=debug,rmsd_threshold = rmsd_threshold, energy_diff = energy_diff)
	if mol3D: # do not go on if molecule has failed
		for i in range (0, len(mol3D.GetConformers())):
			writer.write(mol3D, confId = i)
writer.close()

print ("Done!")


#make sure that partial charge is correct in output file
#need to modify the the name with ring id, I believe, but maybe I am wrong


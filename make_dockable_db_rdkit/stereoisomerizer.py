#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem import AllChem
import chemistry3 as chemistry
import my_mysql3 as mysql
import sys, string
from argparse import ArgumentParser

#I have tried to also count possible EZ stereoisomers, but this is not straightforward, as this gets not assigned by RDKIT


#*******************************************************************************
# MAIN help
#**********
#


description = "Script to generate stereoisomers"

parser = ArgumentParser(description=description)

parser.add_argument("-i", "--input_file", type=str, action="store", dest="infile", help="input file (smiles)",required=True)
parser.add_argument("-o", "--output_file", type=str, action="store", dest="output_file",  help="output file (smiles)",required=True)
parser.add_argument("-m", "--max_centres",  type=int, action="store", dest="maxcenters", default=4, help="Don't generate stereoisomers for molecules that have more stereo centers than specified by maxcenters",required=False)
parser.add_argument("-f", "--flip_all",  action="store_false", dest="flip_only_unasigned", default=True, help="if set stereoisomers for all stereo centers in the molecule will be generated, also if stereoisomery is already specified",required=False)



args = parser.parse_args()
args_dict = vars(args)

#print (args_dict)
locals().update(args_dict) #generates local variables from key / value pairs


junk,type = infile.split('.')
if type == 'smi':
	ifs = Chem.SmilesMolSupplier(infile,delimiter='\t',titleLine=False)
else:
	print ('file type not supported (but you could modify the script to test if it works)')
	sys.exit()


out_file = open(output_file, 'w')



#--------------------------------------------------
def FindMolStereoBonds(mol):
	for bond in mol.GetBonds():
		print (bond.GetStereo())

	centers = 0

	return centers


#------------------------------------------------

for mol in ifs:
	#print (mol.GetProp("_Name"))
	#opts = StereoEnumerationOptions(tryEmbedding=True,unique=True) #Because the molecule is constrained, not all of those isomers can actually exist. We can check that, but this is computationally quite expensive, should save time later
	#maxisomers is only used if not all centers are assigned
	
	#opts = StereoEnumerationOptions(unique=True,maxIsomers=maxcenters*2,onlyUnassigned=flip_only_unasigned,tryEmbedding=True)
	#embedding is really very slow, it's 50% faster, if we test it later and just kick out the ones that failed
	
	opts = StereoEnumerationOptions(unique=True,maxIsomers=maxcenters*2,onlyUnassigned=flip_only_unasigned)


	isomers = tuple(EnumerateStereoisomers(mol,options=opts))


	n_chiral = len(Chem.FindMolChiralCenters(mol,includeUnassigned=False))
	#print (Chem.FindMolChiralCenters(mol))
	n_chiral_all = len(Chem.FindMolChiralCenters(mol,includeUnassigned=True))
	#print (n_chiral, n_chiral_all, maxcenters)



	if len(isomers) == 1:
		s_type ='specified'  
	elif (n_chiral_all - n_chiral) > maxcenters:
		s_type = 'exceeded max centres, random selection'
	else:
		s_type = 'guessed'
	

	for smi in sorted(Chem.MolToSmiles(x) for x in isomers):
		#print(smi, mol.GetProp("_Name"), s_type)
		if s_type != 'specified': #make sure that stereoisomer can actually exist, we will loose a few
			m = Chem.MolFromSmiles(smi)
			m = Chem.AddHs(m)
			ids=AllChem.EmbedMultipleConfs(m, numConfs=1)
			#print ('>>>>',len(ids),ids)
			if len(ids) > 0: #valid stereoisomer was generated
				out_file.write(smi + '\t' + mol.GetProp("_Name") + '\t' + s_type + '\n')
		else: #always keep specified
			out_file.write(smi + '\t' + mol.GetProp("_Name") + '\t' + s_type + '\n')


#make sure to not lose molecules with nothing to do
print ('stereoisomerizer done')





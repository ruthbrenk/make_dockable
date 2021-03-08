#!/usr/bin/env python3

from rdkit import Chem

import sys, string, os, io, shutil, Mol2Writer


#-------------------------------
def solvent_header(out_file, solvent):

  if solvent =='PENT':
    out_file.write('SOLVNT=GENORG IOFR=1.3575 ALPHA=0.00 BETA=0.00 GAMMA=22.3' + '\n')
    out_file.write('& DIELEC=1.84 FACARB=0.00 FEHALO=0.00 ')
  elif solvent =='PFB':
    out_file.write('SOLVNT=GENORG IOFR=1.3777 ALPHA=0.00 BETA=0.00 GAMMA=31.7' + '\n')
    out_file.write('& DIELEC=2.03 FACARB=0.5 FEHALO=0.50 ')
  elif solvent =='CYC':
    #BWQ definitions, but it is actually hexadecane
    out_file.write('SOLVNT=GENORG IOFR=1.4345 ALPHA=0.00 BETA=0.00 GAMMA=38.93' + '\n')
    out_file.write('& DIELEC=2.06 FACARB=0.00 FEHALO=0.00 ')
  elif solvent =='RCYC':
    #BWQ really cyclohexadecane
    out_file.write('SOLVNT=GENORG IOFR=1.4266 ALPHA=0.00 BETA=0.00 GAMMA=35.5' + '\n')
    out_file.write('& DIELEC=2.02 FACARB=0.00 FEHALO=0.00 ')
  elif solvent =='TMB':
    out_file.write('SOLVNT=GENORG IOFR=1.5048 ALPHA=0.00 BETA=0.19 GAMMA=42.0' + '\n')
    out_file.write('& DIELEC=2.37 FACARB=0.67 FEHALO=0.00 ')
  elif solvent == 'DBE':
    out_file.write('SOLVNT=GENORG IOFR=1.3992 ALPHA=0.00 BETA=0.45 GAMMA=32.3' + '\n')
    out_file.write('& DIELEC=3.05 FACARB=0.00 FEHALO=0.00 ' ) 
  elif solvent =='DEE':
    out_file.write('SOLVNT=GENORG IOFR=1.3526 ALPHA=0.00 BETA=0.41 GAMMA=24.0' + '\n')
    out_file.write('& DIELEC=4.24 FACARB=0.00 FEHALO=0.00 ' ) 
  elif solvent =='EAC':
    out_file.write('SOLVNT=GENORG IOFR=1.3723 ALPHA=0.00 BETA=0.45 GAMMA=33.7' + '\n')
    out_file.write('& DIELEC=5.99 FACARB=0.00 FEHALO=0.00 ' ) 
  elif solvent =='DCA':
    out_file.write('SOLVNT=GENORG IOFR=1.4372 ALPHA=0.37 BETA=0.48 GAMMA=41.0' + '\n')
    out_file.write('& DIELEC=7.53 FACARB=0.00 FEHALO=0.00 ' ) 
  elif solvent =='DCE':
    out_file.write('SOLVNT=GENORG IOFR=1.4448 ALPHA=0.10 BETA=0.11 GAMMA=45.9' + '\n')
    out_file.write('& DIELEC=10.19 FACARB=0.00 FEHALO=0.50 ' ) 
  elif solvent =='TBU':
    out_file.write('SOLVNT=GENORG IOFR=1.3878 ALPHA=0.31 BETA=0.60 GAMMA=28.7' + '\n')
    out_file.write('& DIELEC=12.47 FACARB=0.00 FEHALO=0.00 ' ) 
  elif solvent =='ACE':
    out_file.write('SOLVNT=GENORG IOFR=1.3588 ALPHA=0.04 BETA=0.49 GAMMA=23.5' + '\n')
    out_file.write('& DIELEC=20.49 FACARB=0.00 FEHALO=0.00 ' ) 
  elif solvent =='MET':
    out_file.write('SOLVNT=GENORG IOFR=1.3288 ALPHA=0.43 BETA=0.47 GAMMA=22.1' + '\n')
    out_file.write('& DIELEC=32.63 FACARB=0.00 FEHALO=0.00 ' ) 
  elif solvent =='DMSO':
    out_file.write('SOLVNT=GENORG IOFR=1.4170 ALPHA=0.00 BETA=0.88 GAMMA=61.8' + '\n')
    out_file.write('& DIELEC=46.83 FACARB=0.00 FEHALO=0.00 ' )
  else:
    print ("Solvent Keyword is wrong")

 
#-------------------------------
def header(solvent, filename, out_file, charge):

  if solvent == 'Water':
     #out_file.write('\n')  #for amsol 4, we might need this line, but not for amsol7
     out_file.write('SOLVNT=WATER CHARGE=' +  str(charge) + ' AM1 SM5.42R 1SCF GEO-OK CART' + '\n')
  else:
    #out_file.write('\n')   #for amsol4, we might need this line, but not for amsol7
    solvent_header(out_file, solvent)
    out_file.write('CHARGE=' + str(int(charge)) + ' AM1 SM5.42R 1SCF' + '\n')
    out_file.write('& GEO-OK CART' + '\n' )

#-------------------------------
def convert2amsol(mol, solvent, dir):

	filename = mol.GetProp("_Name")
	#print (filename)

	charge = Chem.GetFormalCharge(mol)
	#print (charge)
 
	out_file_wat = open('temp.in-wat', 'w')
	out_file_hex = open('temp.in-hex', 'w')
	
	header('Water', filename, out_file_wat, charge)	
	header(solvent, filename, out_file_hex, charge)

	block = ''

	conf = mol.GetConformer(0)

	coords = conf.GetPositions()
	#print (coords)

	counter = 0
	for atom in mol.GetAtoms():
		atom_name =	atom.GetSymbol()
		line = atom_name 
		for i in range(0,3-len(atom_name)):
			line = line + ' '
		for i in coords[counter]:
			line = line + "%.5f" % i + ' 1 '
		line = line[:-1]
		line = line + '\n'

		block = block + line
		counter = counter + 1


	#print block

	out_file_wat.write(mol.GetProp("_Name") + ' ' + str(mol.GetNumAtoms()) + '\n\n')
	out_file_hex.write(mol.GetProp("_Name") + ' ' + str(mol.GetNumAtoms()) + '\n\n')
	
	out_file_wat.write(block)
	out_file_hex.write(block)

	out_file_hex.write('\n')
	out_file_wat.write('\n')
	
	out_file_wat.close()
	out_file_hex.close()
 
   

#-------------------------------
def Update_charge(mol, partial_charge_list, update_error_list):
  
  counter = 0
    
  if len(partial_charge_list) == mol.GetNumAtoms():
  
    for atom in mol.GetAtoms():
      atom.SetProp('_GasteigerCharge',str(partial_charge_list[counter]))
      counter = counter + 1 

    #print (Mol2Writer.MolToMol2Block(mol,addCharges=False))

  else:
    update_error_list.append(mol.GetProp("_Name"))
  return update_error_list, mol
  
#-------------------------------
def run_amsol():

   #os.system('amsol-mod4 temp.in-wat temp.o-wat')
   
   os.system('amsol7.1 < temp.in-wat > temp.o-wat.tmp')
   os.system('cp fort.12 temp.o-wat')
   
   #os.system('amsol-mod4 temp.in-hex temp.o-hex')

   os.system('amsol7.1 < temp.in-hex > temp.o-hex.tmp')
   os.system('cp fort.12 temp.o-hex')

     
   
#--------------------------------------
def parse_results(mol, result_tab, debug):

  tab_list = []

       
  part_pol_sol_wat = []
  part_apol_sol_wat = []

  failed_amsol_list = []

  num_atoms = mol.GetNumAtoms()
  

  try:
  #if 1 == 1:	

       result_file_wat = open('temp.o-wat', 'r')
       
       result_file_lines = result_file_wat.readlines()
       finished = False
       line_counter = 0

       while not finished:
         #print (result_file_lines[line_counter])  
         if result_file_lines[line_counter].count('kcal') == 3:
         #if result_file_lines[line_counter][62:68] == '(kcal)':
           #data starts in line after following line
           for i in range(line_counter + 2, num_atoms + line_counter + 2):
             #print (result_file_lines[i])
             part_pol_sol_wat.append(float(result_file_lines[i][23:32]))
             part_apol_sol_wat.append(float(result_file_lines[i][52:63]))
             #print (part_pol_sol_wat[-1], part_apol_sol_wat[-1])
	   #next line is total energy
           line_counter = i + 2
           pol_sol_wat = float(result_file_lines[line_counter][23:32])
           apol_sol_wat =float(result_file_lines[line_counter][52:63])
           if debug:
               print ('pol sol wat', pol_sol_wat)
               print ('apol sol wat', apol_sol_wat)
           finished = True
         line_counter = line_counter + 1	   				
       result_file_wat.close()

       #read hex file

              
       result_file_hex = open('temp.o-hex', 'r')
       result_file_lines = result_file_hex.readlines()
       
       part_pol_sol_hex = []
       part_apol_sol_hex = []
       charge_list = []
       area_list = []

       finished = False
       line_counter = 0


       while not finished:

         #print (result_file_lines[line_counter])
         #print (line_counter)


         if result_file_lines[line_counter][6:13] == 'G-P(sol':
           pol_sol_hex = float(result_file_lines[line_counter][65:74])
           if debug :
                print ('pol_sol_hex', pol_sol_hex)
         elif result_file_lines[line_counter][6:16] =='G-CDS(sol)':
           apol_sol_hex =float(result_file_lines[line_counter][65:74])
           if debug :
                print ('apol_sol_hex', apol_sol_hex)
           finished = True	

	 #get correction_factor for non-water solvents
         if result_file_lines[line_counter][4:16] == 'Contribution':
                cs_contribution = float(result_file_lines[line_counter][52:63])
                #print ('cs_contribution', cs_contribution)
                total_area = float(result_file_lines[line_counter][32:41])
                cs_contribution_area = cs_contribution/total_area

         if result_file_lines[line_counter].count('kcal') == 3:
           #data starts in line after following line
           for i in range(line_counter + 2, num_atoms + line_counter + 2):
             part_pol_sol_hex.append(float(result_file_lines[i][23:32]))
             part_apol_sol_hex.append(float(result_file_lines[i][52:63]))
             charge_list.append(float(result_file_lines[i][16:23]))
             area_list.append(float(result_file_lines[i][32:41]))
         line_counter = line_counter + 1	   				
       result_file_hex.close()
       
        


 
       
       #everything is read in
       desolv_pol = pol_sol_wat - pol_sol_hex
       desolv_apol = apol_sol_wat - apol_sol_hex

       correction_factor = cs_contribution/len(charge_list)
       if debug :
             print ('Corrrection factor :' + str(correction_factor))
       
       line = mol.GetProp("_Name") + str(num_atoms).rjust(4)  + '%5.1f'%Chem.GetFormalCharge(mol) + '%9.2f'%desolv_pol + '   000.00' +'%9.2f'%desolv_apol + '   -00.00' + '\n'

       tab_list.append(line)


       counter = 0
       #print (len(charge_list),charge_list)
       #print (len(part_pol_sol_hex))
       #print (len(part_pol_sol_wat))
       for charge in charge_list:

           #correct apol_sol for cs_contribution, relative to area of atom
      	 tab_list.append('%8.4f'%charge + '%8.2f'%(part_pol_sol_wat[counter] - part_pol_sol_hex[counter]) + '  00.30' + '%8.2f'%((part_apol_sol_wat[counter]) - (part_apol_sol_hex[counter] + cs_contribution_area*area_list[counter])) + '   -0.30' + '\n')

      	 counter = counter + 1
       failed_amsol_list,mol = Update_charge(mol,charge_list,failed_amsol_list)
       
       for lines in tab_list:
          result_tab.write(lines)
  
       
  #if 1 == 2:     
  except:
    failed_amsol_list.append(mol.GetProp("_Name"))
  if not debug:
      os.system('rm -f temp.o*')
  else:
      os.system('mv temp.o-hex temp.o-hex_' + mol.GetProp("_Name"))
      os.system('mv temp.o-wat temp.o-wat_' + mol.GetProp("_Name"))
  return failed_amsol_list,mol

#--------------------------------------
def calc_amsol(mol_list, solvent,file,debug):

	ext = '_amsol_cav_' + solvent
	dir = 'temp_files/'

	failed_amsol_list = []

	if not os.path.exists(file[:-4]+ ext):
		os.mkdir(file[:-4] + ext)
		os.chdir(file[:-4] + ext)
		os.mkdir(dir)
	    
		result_tab = open(file[:-4] + '_' + solvent + '.solv' ,'w')
		
		print ('processing file ' + file)

	 
		mol2_text = ''

		for mol in mol_list:


			convert2amsol(mol, solvent, dir)
			
			run_amsol()
			#print (mol)

	
			failed_amsol_list,mol = parse_results(mol, result_tab, debug)
			
			if debug:

				mol2_text = mol2_text + Mol2Writer.MolToMol2Block(mol,addCharges=False)

			
	  
		#clean_up  
		result_tab.close()

		if debug:
			out_file = open(file[:-4] + ext + '.mol2', 'w')
			out_file.write(mol2_text)
			out_file.close()

		if not debug:
			print ('delete files')
			shutil.rmtree(dir,1)
			os.system('rm -f *log')
			os.system('rm -f *.in')
			os.system('rm -f *.out')
			os.system('rm -f core*')
		
		os.chdir('..')
	else:
		print (file[:-4] + ext + ' already exists')

	return failed_amsol_list

#to do: write out mol2 file



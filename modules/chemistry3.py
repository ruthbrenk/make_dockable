#!/usr/bin/env python
#same usefull, reoccuring OEChem functions
#from openeye.oechem import *
import string, re, Mol2Writer, math
import my_mysql3 as mysql

from rdkit import Chem
from rdkit.Chem import AllChem
#############################################################################

def RetrieveMol2Block(fileLikeObject, delimiter="@<TRIPOS>MOLECULE"):
    """generator which retrieves one mol2 block at a time
    """
    mol2 = []
    #print (delimiter)
    for line in fileLikeObject:
        #print (line)
        if line.startswith(delimiter) and mol2:
            yield "".join(mol2)
            mol2 = []
        mol2.append(line)
    #print (mol2)
    if mol2:
        yield "".join(mol2)
#############################################################################

def Mol2Supplier (infile, sanitize=True):   #infile: name of closed file

	ifs = []
	
	for mol2 in RetrieveMol2Block(open(infile)):
		#print (mol2)
		mol = Chem.MolFromMol2Block(mol2, sanitize)
		ifs.append(mol)
	return ifs


#############################################################################
#iso: check for e/z 
#kek: kekulise 
#chiral: check for r/s
def CanSmi(mol,iso,kek,chiral):
  smiflag=OESMILESFlag_Canonical
  #next two lines to not have double bond problem
  OEFindRingAtomsAndBonds(mol)
  #OEAssignAromaticFlags(mol,OEAroModelOpenEye)
  OEAssignAromaticFlags(mol,OEAroModelMMFF)
	
  if iso:
    smiflag|=OESMILESFlag_ISOMERIC
    OEPerceiveChiral(mol)
    for bond in mol.GetBonds():
      if bond.GetOrder()!=2: continue
      HasStereoSpecified=bond.HasStereoSpecified()
      IsChiral=bond.IsChiral()
      if HasStereoSpecified and not IsChiral:
        OEThrow.Warning("correcting bond stereo inconsistency!")
        a1=bond.GetBgn()
        a2=bond.GetEnd()
        bond.SetStereo([a1,a2],OEBondStereo_CisTrans,OEBondStereo_Undefined)
  else:	
    #remove E/Z information
    for bond in mol.GetBonds():
      if bond.GetOrder()!=2: continue
      a1=bond.GetBgn()
      a2=bond.GetEnd()
      bond.SetStereo([a1,a2],OEBondStereo_CisTrans,OEBondStereo_Undefined)
  if not chiral:
    OEPerceiveChiral(mol)
    for atom in mol.GetAtoms():
      IsChiral=atom.IsChiral()
      HasStereo=atom.HasStereoSpecified()

      if IsChiral:
        nbrs=[]
        for bond in atom.GetBonds():
          nbrs.append(bond.GetNbr(atom))
        atom.SetStereo(nbrs,OEAtomStereo_Tetrahedral,OEAtomStereo_Undefined)
  else:
    smiflag|=OESMILESFlag_ISOMERIC
    OEPerceiveChiral(mol)
    for atom in mol.GetAtoms():
      IsChiral=atom.IsChiral()
      HasStereo=atom.HasStereoSpecified()
      if HasStereo and not IsChiral:
        OEThrow.Warning("correcting atom stereo inconsistency!")
        nbrs=[]
        for bond in atom.GetBonds():
          nbrs.append(bond.GetNbr(atom))
        atom.SetStereo(nbrs,OEAtomStereo_Tetrahedral,OEAtomStereo_Undefined)
      #elif  IsChiral and not HasStereo:
        #OEThrow.Warning("R/S not specified!")
       # nbrs=[]
        #for bond in atom.GetBonds():
          #nbrs.append(bond.GetNbr(atom))
        #atom.SetStereo(nbrs,OEAtomStereo_Tetrahedral,OEAtomStereo_Undefined)


  if kek:
    OEFindRingAtomsAndBonds(mol)
    OEAssignAromaticFlags(mol,OEAroModelOpenEye)
    for bond in mol.GetBonds(OEIsAromaticBond()):
      bond.SetIntType(5)
    OECanonicalOrderAtoms(mol)
    OECanonicalOrderBonds(mol)
    OEClearAromaticFlags(mol)
    OEKekulize(mol)
  #smiflag|=OESMILESFlag_Hydrogens
  smi=OECreateSmiString(mol,smiflag)
  return smi

#############################################################################
#check if molecule is planar
def check_planar(mol, atom_idx_list,debug):
	planar = False

	if len(atom_idx_list) > 3:
		conf = mol.GetConformer(0)
		coords = conf.GetPositions()
		#print len(coords), mol.NumAtoms()
		#coords ={0: (2, -1, -2), 1: (1,2,1), 2: (2,3,0), 3: (5,0,-6)} 
		#coords = {0:(-1,3,3), 1: (0,4,2), 2: (3,1,-4)}
		#print (coords)
		sub_coords = []
		
		for atom in mol.GetAtoms():
			#print (atom.GetIdx(), atom.GetIsAromatic())
			if atom.GetIdx() in atom_idx_list:
				sub_coords.append(coords[atom.GetIdx()])
		#print (sub_coords)

		#a = 0 -1
		a= []
		for i in range(0,3):
			a.append(sub_coords[0][i] - sub_coords[1][i])
		#b = 0 -2
		b= []
		for i in range(0,3):
			b.append(sub_coords[0][i] - sub_coords[2][i])
		planar = True
		#print (a,b)

		for i in range(3, len(sub_coords)):
			#c = 0 - i
			c= []
			for j in range(0,3):
				#print coords[i][j]
				c.append(sub_coords[0][j] - sub_coords[i][j])

			#check for planarity
			determinant = a[0]*b[1]*c[2] + c[0]*a[1]*b[2] + b[0]*c[1]*a[2] - (a[2]*b[1]*c[0]) - (a[0]*c[1]*b[2]) - (b[0]*a[1]*c[2])
			#print a[0], b[1], c[2] 
			#determinant = - b[0]*a[1]*c[2]

			#print determinant
			#print a
			#print b
			#print c
			if determinant > 0.1 or determinant < -0.1:
				planar = False
	else:
		#three atoms are always in one plane
		planar = True
	if debug:
		print (planar, 'planar')
	return planar

#--------------------------------------------
def get_smartslist(table):
	conn=mysql.connect2server('', 'webuser', 'purchasable')  			  
	cursor = conn.cursor ()
	rule = []
	command = 'select smarts_string, name from purchasable.' + table + ' order by name'
	print (command)
	cursor.execute(command)	
	rules = cursor.fetchall()
	for i in rules:
		smarts = i[0]
		name = i[1]
		name = name.replace(' ', '_')
		name = name.replace('-', '_')
		name = name.replace('/', '_')
		name = name.replace(',', '_')
		name = name.replace('=', '_')
		rule.append([smarts,name])
	return rule
#--------------------------------------------
def get_smartslist_file(infile):
	rule = []
	for i in infile.readlines():
		i = i.strip()
		smarts,name = i.split('\t')
		rule.append([smarts,name])
	infile.close()
	return rule
#--------------------------------------------
def get_patlist(smarts_list, atom_expr, bond_expr,max):
  patlist=[]  ### list of (pat,smarts,description) tuples
  for i in smarts_list:
      smarts=i[0]
      #print (smarts, 'get_patlist')
      pat=Chem.MolFromSmarts(smarts)
      patlist.append((pat,smarts,i[1]))
   
  return patlist
#--------------------------------------------
def Match(pat,smarts,mol,kek,verbose,exph,asym,usa):
  targetmol=OEGraphMol()
  idxll=[]  ### list of list of match symclasses for a mol
  vout='MOL: %s %s\n'%(OECreateCanSmiString(mol),mol.GetTitle())
  vout+='SMARTS: %s\n'%smarts
  if exph: OEAddExplicitHydrogens(mol)
  ## OEAssignHybridization(mol)
  OETriposAtomNames(mol)
  if kek: 
     OEClearAromaticFlags(mol)
     OEKekulize(mol)
     #print OECreateCanSmiString(mol)
  matchcount=0
  if asym: OEPerceiveSymmetry(mol)
  for match in pat.Match(mol): #used to be (mol,usa), does not work anylonager
     vout+='\t'
     if usa: vout+='USA '
     vout+='match (%d): '%(matchcount+1)
     idxl=[]  ### list of atom idxs for a match
     for matchpair in match.GetAtoms():
        vout+='%s '%matchpair.target.GetName()
        idxl.append(matchpair.target.GetSymmetryClass())

     idxll.append(idxl)
     matchcount+=1
     OESubsetMol(targetmol,match)
     vout+='\t%s'%(OECreateCanSmiString(targetmol))
     targetmol.Clear()
  if matchcount==0:
    vout+='no match\n'

  if asym and matchcount>1:
    syms=[]
    for j in range(matchcount):
      for k in range(j):
        if idxll[j]==idxll[k]:
          vout+='symmetrical match (%d~%d)'%(j+1,k+1)
          if j not in syms: syms.append(j)
    vout+='\tasymmetrical matches: %d'%(matchcount-len(syms))

  if verbose:
    sys.stderr.write('%s\n'%vout)

  return matchcount
#--------------------------------------------
def clean_smiles(smi):
    old_smi = ''
    lb = 0
    rb = 0
    #print '-------------------->'
    while (smi.find('[',lb) != -1):
         old_smi = smi
         lb = smi.find('[',lb)
         rb = smi.find(']',rb)
         if lb + 2 == rb:
                #print smi[:lb -1]
                  #print smi[lb + 1]
                  #print smi[rb + 1:]
                  smi = smi[:lb] + smi[lb + 1] + smi[rb + 1:]
         rb = rb + 1
         lb = lb + 1
    #print smi
    smi = smi.replace('[CH]','C')
    smi = smi.replace('[NH]','N')
    smi = smi.replace('[NH2]','N')
    smi = smi.replace('[cH]','c')
    smi = smi.replace('[CH2]','C')
    #smi = smi.replace('c-','c')
    smi = smi.replace('-c','c')
    #single bonds in aromatic rings must be replaced
    for counter in range(1,11):
         digit = str(counter)
         smi = smi.replace('c-'+ digit, 'c' + digit)  ##added c,n because otherwise charge can get lost, 20/02/2008
         smi = smi.replace('n-'+ digit, 'n' + digit)
         smi = smi.replace('s-'+ digit, 's' + digit)
    #smi = smi.replace('s-','s')
    smi = smi.replace('-s','s')
    #smi = smi.replace('n-','n') 
    smi = smi.replace('-n','n') 
    #smi = smi.replace('o-','o')
    smi = smi.replace('-o','o')
	
    return smi		

#-----------------------------------
#convert smiles to smarts
def convert_smarts(smarts_list):
	#exocyclic carbons can be aliphatic or aromatic
	for i in range(0,len(smarts_list)):
		smiles = smarts_list[i][0]
		smi_mol = OEGraphMol()
		OEParseSmiles(smi_mol,smiles)
		for atom in smi_mol.GetAtoms():
			if (atom.IsCarbon() and atom.GetExplicitDegree() == 1 and not atom.IsInRing()) or (not atom.IsInRing() and atom.IsAromatic()) : #terminal carbon or (non ring + aromatic)
				#mark atom with unique element => easier to change later
				atom.SetAtomicNum(20)
				atom.SetAromatic(0)
				atom.SetImplicitHCount(0)
				atom.SetFormalCharge(0)
				for bond in atom.GetBonds():
					bond.SetOrder(0)
		smarts = OECreateSmiString(smi_mol)
		smarts = smarts.replace('Ca', 'C,c')
		smarts = smarts.replace('[nH]','n')	#tautomers
		smarts = smarts.replace('[NH+]','N')	#charge
		smarts = smarts.replace('[nH+]','n')	#charge
		smarts = smarts.replace('[N+]','N')	#charge
		smarts = smarts.replace('[O+]','O')	#charge
		smarts = smarts.replace('[NH2+]','N')	#charge
		smarts = smarts.replace('[NH3+]','N')	#charge
		#smarts = smarts.replace('[CH]','C')	#why ?
		smarts = smarts.replace('-','')	#charge
		smarts = smarts.replace('CH3','C')	#why, (maybe implicit h count) ?
		smarts = smarts.replace('CH2','C')	#why, (maybe implicit h count) ?
		smarts = smarts.replace('H','')	#why, (maybe implicit h count) ?
		debug = False
		if smarts.find('+') >= 0 and debug:
			print ('charge issue', smarts)
		smarts_list[i][0] = smarts
	return smarts_list
#--------------------------------------------
def check_smi(smi,old_smi, rule, comment, title, debug):
   error = False
   if smi.find('.') != -1:
     print ('----> error <------')
     print (rule+', '+  comment+', '+  title)
     print (smi)
     error = True
   if debug and smi != old_smi:
      print ('----> reaction <------')	
      print (old_smi)
      print ('--> ' + title  +', '+ rule +', ' + comment)
      print (smi)
      print ('----> end <------')
      print
   return error
#--------------------------------------------
#### lneutralize molecules
#--------------------------------------------
def neutralize_get_rules():
	conn=mysql.connect2server('', 'webuser','purchasable')
	cursor = conn.cursor ()
	command = "SELECT smirks, type, pka, comment, id FROM purchasable.definition_neutralize_compounds order by id"
	print (command)
	cursor.execute(command)
	rows = cursor.fetchall ()

	rules = []
	for row in rows:
		rules.append([[row[0]],row[1],float(row[2]),row[3]])

	#OPENEYE has got a problem with writting [NH3+], ...
	rule_counter = 0
	for rule in rules:
		trans_counter = 0
		for transformation in rule[0]:
			educt, product = transformation.split('>>')
			#primary amines
			product = product.replace('[NH3+]', '[N+]([H])([H])[H]')
			#secundary amines
			product = product.replace('[NH2+]', '[N+]([H])([H])')
			#tertary amines
			product = product.replace('[NH+]', '[N+]([H])')
			product = product.replace('[nH+]', '[n+]([H])')
			product = product.replace('[nH]', '[n]([H])')
			reaction = educt + '>>' + product
			#print reaction
			rules[rule_counter][0][trans_counter] = reaction
			trans_counter = trans_counter + 1
		rule_counter = rule_counter + 1
	return rules
#--------------------------------------------
def neutralize_smi(smiles,rules,mol,debug):

  iso = False #keep stereoinformation
  kek = False #kekulize to get unique smiles (no, this is no longer true)
  chiral = False #keep information about chirality (could also be false, because this information is guessed anyway)

  for rule in rules:
         old_smi = smiles
         for transformation in rule[0]:
                  #make mol
                  work_mol = OEGraphMol()
                  OEParseSmiles(work_mol, smiles)
                  umr = OEUniMolecularRxn(transformation)
                  umr(work_mol)
                  smi = CanSmi(work_mol,iso,kek, chiral)
                  error = check_smi(smi, old_smi, transformation, rule[3], mol.GetTitle(), debug)
                  #smi_kek = chemistry.CanSmi(work_mol,iso,True,chiral)
                  if not error:
                           smiles = smi
                  elif error:
                           smiles = old_smi
  return smiles
#--------------------------------------------
def swap_stereoisomers(smiles):
#/.../ = \...\ => generate opposite verions

    	#generate original smiles
	smiles = mysql.decode_smiles(smiles)

	smiles = smiles.replace('/', 'X')
	smiles = smiles.replace('\\', '/')
	smiles = smiles.replace('X', '\\')

	smiles = mysql.encode_smiles(smiles)
	return smiles

#--------------------------------------------
def generateconformations(m, n,debug,rmsd_threshold = 1.0, energy_diff = 10.0):


	m = Chem.AddHs(m)
	#print ("Embed")
	ids=AllChem.EmbedMultipleConfs(m, numConfs=n, pruneRmsThresh=rmsd_threshold)  #pruneRmsThresh=0.5: minium rmsd between two molis, used 0.6 for omega
	if len(m.GetConformers()) == 0: #conformer generation failed
		return False

	id_list = list(ids)
	#res = AllChem.MMFFOptimizeMoleculeConfs(m,mmffVariant='MMFF94s',maxIters=10000) #don't get N.pl3 planar
	#print ("UFFOpt")
	res = AllChem.UFFOptimizeMoleculeConfs(m,maxIters=10000)

	#print (res)

	#get lowest energy
	energy_list = []
	for entry in res:
		#print (entry)
		if entry[0] ==0: #minimization converged
			energy_list.append(entry[1])
	energy_list.sort()
	#print (energy_list)

	removed = False

	for i in range (0, len(res)):
		entry = res[i]
		#print (energy_list[0],  energy_list[0] + energy_diff, entry[1])
		if entry[0] ==0 and (entry[1] < (energy_list[0] + energy_diff)):
			continue #keep this conformer
		else:
			#print (id_list)
			m.RemoveConformer(i) #delete the conformer
			id_list.remove(i) #delete from list
			if debug:
				print ('removed conformer', i)
				#print (entry [0], energy_list[0],  energy_list[0] + energy_diff, entry[1])
				#print (id_list)
			removed = True
	if removed: #renumber conformers
		counter = 0
		for conf in m.GetConformers():
			conf.SetId(counter)
			counter = counter + 1
	if debug:
		print (len(m.GetConformers()), "conformers generated")
		print ('rmsd cut-off:', rmsd_threshold)

	#Mol2Writer.MultiMolToMol2File([m], m.GetProp("_Name") + 'MMFFOptimizeMoleculeConfs.mol2', confId=None, addHs=True)
		

#The result is a list a containing 2-tuples: (not_converged, energy) for each conformer. If not_converged is 0, the minimization for that conformer converged

	# EmbedMultipleConfs returns a Boost-wrapped type which 
	# cannot be pickled. Convert it to a Python list, which can.
	return m
#-----------------------------
def gen_confs(smi_input_file, max_confs,debug=False,rmsd_threshold = 1.0, energy_diff = 10.0 ):
#takes an input_file

	input_smiles = Chem.SmilesMolSupplier(smi_input_file, delimiter='\t', titleLine=False)

	molis = ''

	mol_list = []


	for mol in input_smiles:
			if debug:
				print (mol.GetProp("_Name"))
			rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
			sat_rings = Chem.rdMolDescriptors.CalcNumSaturatedRings(mol)
			if mol:
				n = 1 + math.pow(5,rot_bonds)  #1 => also for aromatic rings without substituents
				if sat_rings > 0:
					n = n * (sat_rings * 6 * 3) #1 => also for aromatic rings without substiuents

				n = int(n) #math.pow returns a float
			
				if n > max_confs:
					n = max_confs

				if debug:
					print ('n', n)

				m = generateconformations(mol, n,debug,rmsd_threshold=rmsd_threshold)
				if m: #do not append failed molecules
					mol_list.append(m)
				#for id in ids:
					#writer.write(m, confId=id)
					#molis = molis + Chem.rdmolfiles.MolToMolBlock(m, confId=id)
	if debug:
		print ('finshed conf gen')

	return mol_list
#-----------------------------------
def gen_confs_single(mol, max_confs,debug=False,rmsd_threshold = 1.0, energy_diff = 10.0 ):
#takes a singe molecules as input => does not need to keep so much in memoery

	if debug:
		print (mol.GetProp("_Name"))
	rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
	sat_rings = Chem.rdMolDescriptors.CalcNumSaturatedRings(mol)
	n = 1 + math.pow(5,rot_bonds)  #1 => also for aromatic rings without substituents
	if sat_rings > 0:
		n = n * (sat_rings * 6 * 3) #1 => also for aromatic rings without substiuents
	n = int(n) #math.pow returns a float
		
	if n > max_confs:
		n = max_confs
	if debug:
		print ('n', n)
	m = generateconformations(mol, n,debug,rmsd_threshold=rmsd_threshold)

	if debug:
		print ('finshed conf gen')

	return m
#-----------------------------------		





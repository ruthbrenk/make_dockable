#!/usr/bin/env python
#same usefull, reoccuring OEChem functions
from openeye.oechem import *
import string, re
import my_mysql as mysql


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
def check_planar(mol):
	planar = False

	if mol.NumAtoms() > 3:
		coords = mol.GetCoords()
		#print len(coords), mol.NumAtoms()
		#coords ={0: (2, -1, -2), 1: (1,2,1), 2: (2,3,0), 3: (5,0,-6)} 
		#coords = {0:(-1,3,3), 1: (0,4,2), 2: (3,1,-4)}
		#print coords

		#a = 0 -1
		a= []
		for i in range(0,3):
			a.append(coords[0][i] - coords[1][i])
		#b = 0 -2
		b= []
		for i in range(0,3):
			b.append(coords[0][i] - coords[2][i])
		planar = True

		for i in range(3, len(coords)):
			#c = 0 - i
			c= []
			for j in range(0,3):
				#print coords[i][j]
				c.append(coords[0][j] - coords[i][j])

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
		#three atoms are allways in one plane
		planar = True
	#print planar
	return planar

#--------------------------------------------
def get_smartslist(table):
	conn=mysql.connect2server('', 'webuser', 'purchasable')  			  
	cursor = conn.cursor ()
	rule = []
	command = 'select smarts_string, name from purchasable.' + table + ' order by name'
	print command
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
	for i in infile.xreadlines():
		#print i
		smarts,name = string.split(string.strip(i), '\t')
		rule.append([smarts,name])
	infile.close()
	return rule
#--------------------------------------------
def get_patlist(smarts_list, atom_expr, bond_expr,max):
  patlist=[]  ### list of (pat,smarts,description) tuples
  for i in smarts_list:
      smarts=i[0]
      #print smarts
      pat=OESubSearch()
      if not pat.Init(re.sub('\s.*$','',smarts)):
        OEThrow.Fatal('Bad smarts: '+smarts)

      if atom_expr>=0 or bond_expr>=0:
        if not atom_expr>=0: atom_expr=OEExprOpts_DefaultAtoms
        if not bond_expr>=0: bond_expr=OEExprOpts_DefaultBonds
        qmol=pat.GetPattern()
        qmol.BuildExpressions(atom_expr,bond_expr)
        pat=OESubSearch(qmol)

      if max: pat.SetMaxMatches(max)
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
    while (smi.find('[',lb) <> -1):
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
			print 'charge issue', smarts
		smarts_list[i][0] = smarts
	return smarts_list
#--------------------------------------------
def check_smi(smi,old_smi, rule, comment, title, debug):
   error = False
   if smi.find('.') <> -1:
     print '----> error <------'
     print rule+', '+  comment+', '+  title
     print smi
     error = True
   if debug and smi <> old_smi:
      print'----> reaction <------'	
      print old_smi
      print '--> ' + title  +', '+ rule +', ' + comment
      print smi
      print '----> end <------'
      print
   return error
#--------------------------------------------
#### lneutralize molecules
#--------------------------------------------
def neutralize_get_rules():
	conn=mysql.connect2server('', 'webuser','purchasable')
	cursor = conn.cursor ()
	command = "SELECT smirks, type, pka, comment, id FROM purchasable.definition_neutralize_compounds order by id"
	print command
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

	smiles = smiles.replace('/', 'T')
	smiles = smiles.replace('\\', '/')
	smiles = smiles.replace('T', '\\')

	smiles = mysql.encode_smiles(smiles)
	return smiles

#--------------------------------------------

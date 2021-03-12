#!/usr/bin/env python
import sys,os
from rdkit import Chem
import rdkit.Chem.rdMolTransforms
import math

def list_dihedrals(mol):
    RotatableBond = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
    rotatablebonds = mol.GetSubstructMatches(RotatableBond)
    dihedrals = list()
    for bond in rotatablebonds:
        j,k = bond
        atomj = mol.GetAtomWithIdx(j)
        atomk = mol.GetAtomWithIdx(k)
        jneighbors = atomj.GetNeighbors()
        kneighbors = atomk.GetNeighbors()
        i = None
        for tmp in jneighbors:
            if tmp.GetIdx() != k :
                i = tmp.GetIdx()
                break
        l = None
        for tmp in kneighbors:
            if tmp.GetIdx() != j :
                l = tmp.GetIdx()
                break
        if i != None and l != None:
            namei = mol.GetAtomWithIdx(i).GetPropsAsDict()['_TriposAtomName']
            namej = mol.GetAtomWithIdx(j).GetPropsAsDict()['_TriposAtomName']
            namek = mol.GetAtomWithIdx(k).GetPropsAsDict()['_TriposAtomName']
            namel = mol.GetAtomWithIdx(l).GetPropsAsDict()['_TriposAtomName']
            allname = f"{namei}-{namej}-{namek}-{namel}"
            dihedrals.append( (allname,(i,j,k,l))  )
    return dihedrals

def compute_d(mol,dihedral):
    conformer = mol.GetConformers()[0]
    i,j,k,l = dihedral
    result = rdkit.Chem.rdMolTransforms.GetDihedralDeg(conformer, i,j,k,l)
    epsilon = 1e-3
    if result <= 0 - math.pi + epsilon :
        result += math.pi
    elif result > math.pi:
        result -= math.pi
    else:
        pass
    return (result)

def next_mol2_block(infile):
    block = ''
    for line in open(infile):
        if "@<TRIPOS>MOLECULE" in line:
            if block == '':
                pass
            else:
                yield block
                block = ''
            block = block + line
        else:
            block = block + line
    yield block

infile = sys.argv[1]

mols = [ Chem.MolFromMol2Block(x,sanitize=True) for x in next_mol2_block(infile) ]
dihedrals = list_dihedrals(mols[0])

names = [ x[0] for x in dihedrals ]
print("\t".join(names))
for mol in mols:
    values = [ "%8.3f"%compute_d(mol,x) for (i,x) in dihedrals ]
    print("\t".join(values))


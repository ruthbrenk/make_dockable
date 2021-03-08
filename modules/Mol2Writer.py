# $Id$
#
# Copyright (C) 2001-2008 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit.Chem.rdmolfiles import MolFromSmarts, MolFromMol2Block
from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges
from rdkit.Chem.rdmolops import AddHs
from rdkit import Chem


def _get_positions(mol, confId=-1):
    if mol.GetNumConformers() > 0:
        conf = mol.GetConformer(confId)
        return [conf.GetAtomPosition(i) for i in range(conf.GetNumAtoms())]
    return [(0, 0, 0) for i in range(mol.GetNumAtoms())]


class Mol2MolSupplier:
    def __init__(self, filename, *args, **kwargs):
        """Reads a multi-mol Mol2 file
          ARGUMENTS:

            - filename: the file to read  or file-like object
            - args, kwargs: arbitrary arguments to pass to internal MolFromMol2Block

          RETURNS:

            None
        """
        self.f = filename
        self._args = args
        self._kwargs = kwargs

    def __iter__(self):
        """ Iterates over molecules in file """
        block = ''
        data = ''
        n = 0
        if hasattr(self.f, 'read') and hasattr(self.f, 'close'):
            f = self.f
        else:
            f = open(self.f)
        for line in f:
            if line[:1] == '#':
                data += line
            elif line[:17] == '@<TRIPOS>MOLECULE':
                if n > 0:  #skip `zero` molecule (any preciding comments and spaces)
                    yield MolFromMol2Block(block, *self._args, **self._kwargs)
                n += 1
                block = data
                data = ''
            block += line
        # open last molecule
        if block:
            yield MolFromMol2Block(block, *self._args, **self._kwargs)
        f.close()


class Mol2Writer:
    def __init__(self, filename, *args, **kwargs):
        """Writes a multi-mol Mol2 file
          ARGUMENTS:

            - filename: the file to write or file-like object
            - args, kwargs: arbitrary arguments to pass to internal MolToMol2Block

          RETURNS:

            None
        """
        if hasattr(filename, 'write') and hasattr(self.f, 'close'):
            self.f = filename
        else:
            self.f = open(filename, 'w')
        self._args = args
        self._kwargs = kwargs

    def write(self, mol):
        """Writes a multi-mol Mol2 file
          ARGUMENTS:

            - mol: the molecule to be written

          RETURNS:

            bool
        """
        return self.f.write(MolToMol2Block(mol, *self._args, **self._kwargs))

    def close(self):
        """ Closes file for writing """
        return self.f.close()


def MolToMol2File(mol, filename, confId=-1, addHs=True):
    """Writes a Mol2 file for a molecule
      ARGUMENTS:

        - mol: the molecule
        - filename: the file to write to
        - confId: (optional) selects which conformation to output (-1 = default)
                  if set to None will return all conformers

      RETURNS:

        None
    """
    block = MolToMol2Block(mol, confId, addHs)
    open(filename, "w").writelines(block)


#added by me
def MultiMolToMol2File(mol_list, filename, confId=-1, addHs=True, append=False, addCharges=True):
    """Writes a multi Mol2 file for a molecule ensemble
      ARGUMENTS:

        - mol: the molecule ensemble
        - filename: the file to write to
        - confId: (optional) selects which conformation to output (-1 = default)
                  if set to None will return all conformers

      RETURNS:

        None
    """
    all_mols_block = ''
    for mol in mol_list:
        block = MolToMol2Block(mol, confId, addHs, addCharges)
        all_mols_block = all_mols_block + block
    if not append:
    	open(filename, "w").writelines(all_mols_block)
    else:
    	open(filename, "a").writelines(all_mols_block)

#-------------------------------


def MolToMol2Block(mol, confId=-1, addHs=True, addCharges=True):
    """Returns a Mol2 string block for a molecule
      ARGUMENTS:

        - mol: the molecule
        - confId: (optional) selects which conformation to output (-1 = default)
                  if set to None will return all conformers

      RETURNS:

        a string
    """

    #
    # References
    # - Format specs http://www.tripos.com/data/support/mol2.pdf
    # - Atom typing http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
    #

    confIds = (confId, )

    if confId == None:
        confIds = [] #RB: this is a bit dirty, because I make a list out of a tuple
        confId = mol.GetNumConformers()
        for i in range (0,confId):
               confIds.append(i)

    blocks = []

    # add explicit hydrogens (since mol2 reader requires them)
    if addHs:
        h_coords = mol.GetNumConformers() > 0 and mol.GetConformer(-1).Is3D()
        try:
            mol = AddHs(mol, addCoords=h_coords)
        except RuntimeError:
            mol = AddHs(mol, addCoords=False)

    # compute charges
    if addCharges:
        ComputeGasteigerCharges(mol)

    for confId in confIds:
        molecule = """@<TRIPOS>MOLECULE
{}
{} {} 0 0 0
SMALL
GASTEIGER\n\n""".format(
            mol.GetProp("_Name") if mol.HasProp("_Name") else "UNK",
            mol.GetNumAtoms(), mol.GetNumBonds())

        # FIXME "USER_CHARGES" could become 'Gasteiger charges'
        # FIXME "SMALL" means small molecule but could become "PROTEIN"

        pos = _get_positions(mol, confId)
        atom_lines = [
            "{:>4} {:>4} {:>13.4f} {:>9.4f} {:>9.4f} {:<5} {} {} {:>7.4f}".
            format(a.GetIdx() + 1,
                   a.GetSymbol(),
                   float(pos[a.GetIdx()][0]),
                   float(pos[a.GetIdx()][1]),
                   float(pos[a.GetIdx()][2]),
                   _sybyl_atom_type(a), 1, "UNL",
                   float(a.GetProp('_GasteigerCharge').replace(',', '.'))
                   if a.HasProp('_GasteigerCharge') else 0.0)
            for a in mol.GetAtoms()
        ]
        atom_lines = ["@<TRIPOS>ATOM"] + atom_lines
        atom_lines = "\n".join(atom_lines) + "\n"

        bond_lines = [
            "{:>5} {:>5} {:>5} {:>2}".format(
                bid + 1,
                b.GetBeginAtomIdx() + 1,
                b.GetEndAtomIdx() + 1, "ar"
                if b.GetBondTypeAsDouble() == 1.5 else "am"
                if _amide_bond(b) else str(int(b.GetBondTypeAsDouble())))
            for bid, (b) in enumerate(mol.GetBonds())
        ]
        bond_lines = ["@<TRIPOS>BOND"] + bond_lines + ["\n"]
        bond_lines = "\n".join(bond_lines)

        block = molecule + atom_lines + bond_lines
        blocks.append(block)
    return "".join(blocks)


def _sybyl_atom_type(atom):
    """ Asign sybyl atom type
    Reference #1: http://www.tripos.com/mol2/atom_types.html
    Reference #2: http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
    """
    sybyl = None
    atom_symbol = atom.GetSymbol()
    atomic_num = atom.GetAtomicNum()
    hyb = atom.GetHybridization() - 1  # -1 since 1 = sp, 2 = sp1 etc
    hyb = min(hyb, 3)
    degree = atom.GetDegree()
    aromtic = atom.GetIsAromatic()

    # define groups for atom types
    guanidine = '[NX3,NX2]([!O,!S])!@C(!@[NX3,NX2]([!O,!S]))!@[NX3,NX2]([!O,!S])'  # strict
    # guanidine = '[NX3]([!O])([!O])!:C!:[NX3]([!O])([!O])' # corina compatible
    # guanidine = '[NX3]!@C(!@[NX3])!@[NX3,NX2]'
    # guanidine = '[NX3]C([NX3])=[NX2]'
    # guanidine = '[NX3H1,NX2,NX3H2]C(=[NH1])[NH2]' # previous
    #

    if atomic_num == 6:
        if aromtic:
            sybyl = 'C.ar'
        elif degree == 3 and _atom_matches_smarts(atom, guanidine):
            sybyl = 'C.cat'
        else:
            sybyl = '%s.%i' % (atom_symbol, hyb)
    elif atomic_num == 7:
        if aromtic:
            sybyl = 'N.ar'
        elif _atom_matches_smarts(atom, 'C(=[O,S])-N'):
            sybyl = 'N.am'
        elif degree == 3 and _atom_matches_smarts(atom, '[$(N!-*),$([NX3H1]-*!-*)]'):
            sybyl = 'N.pl3'
        elif _atom_matches_smarts(atom, guanidine):  # guanidine has N.pl3
            sybyl = 'N.pl3'
        elif degree == 4 or hyb == 3 and atom.GetFormalCharge():
            sybyl = 'N.4'
        else:
            sybyl = '%s.%i' % (atom_symbol, hyb)
    elif atomic_num == 8:
        # http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
        if degree == 1 and _atom_matches_smarts(atom, '[CX3](=O)[OX1H0-]'):
            sybyl = 'O.co2'
        elif degree == 2 and not aromtic:  # Aromatic Os are sp2
            sybyl = 'O.3'
        else:
            sybyl = 'O.2'
    elif atomic_num == 16:
        # http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
        if degree == 3 and _atom_matches_smarts(atom, '[$([#16X3]=[OX1]),$([#16X3+][OX1-])]'):
            sybyl = 'S.O'
        # https://github.com/rdkit/rdkit/blob/master/Data/FragmentDescriptors.csv
        elif _atom_matches_smarts(atom, 'S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]'):
            sybyl = 'S.o2'
        else:
            sybyl = '%s.%i' % (atom_symbol, hyb)
    elif atomic_num == 15 and hyb == 3:
        sybyl = '%s.%i' % (atom_symbol, hyb)

    if not sybyl:
        sybyl = atom_symbol
    return sybyl


def _atom_matches_smarts(atom, smarts):
    idx = atom.GetIdx()
    patt = MolFromSmarts(smarts)
    for m in atom.GetOwningMol().GetSubstructMatches(patt):
        if idx in m:
            return True
    return False


def _amide_bond(bond):
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()
    if ((a1.GetAtomicNum() == 6 and
         a2.GetAtomicNum() == 7 ) or
        (a2.GetAtomicNum() == 6 and
         a1.GetAtomicNum() == 7)):
        # https://github.com/rdkit/rdkit/blob/master/Data/FragmentDescriptors.csv
        patt = MolFromSmarts('C(=O)-N')
        for m in bond.GetOwningMol().GetSubstructMatches(patt):
            if a1.GetIdx() in m and a2.GetIdx() in m:
                return True
    return False


__all__ = ["MolToMol2Block", "MolToMol2File"]


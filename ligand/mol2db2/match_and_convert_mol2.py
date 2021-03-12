#!/usr/bin/env python3.6
import sys, os
from rdkit import Chem

### Define const
babel_bin = os.environ["OBABELBASE"] + "/obabel"
dockbase = os.environ["DOCKBASE"]
try:
    matchversion = os.environ['MATCHVERSION']
except:
    matchversion = 'Highest'

python2_bin = "/usr/bin/python2"
EPS = 0.05
EPS_modify = 1e-4

infile = sys.argv[1]
outprefix = sys.argv[2]


def load_coords(ref, conf, index):
    """ Load conformer coordinates"""
    coords = list()
    conf = ref.GetConformer(id=conf)
    for i in index:
        c = conf.GetAtomPosition(i)
        x, y, z = c.x, c.y, c.z
        coords.append((x, y, z))
    print("Return coords:", coords)
    return coords


def dist_2(c1, c2):
    """ Calculate the distance**2 between two points"""
    return (c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2


def is_same(ref, i1, i2, index, eps=EPS, eps_modify=EPS_modify):
    """ A function that checks whether two conformers share the same coordinates for a part of the whole molecule.
        Here, two cutoffs are applied:
        EPS: The cutoff to define two coordinate are different.
        EPS_modify: two coordinate distance > EPS_modify are also regarded as same. This script will assign the exact same coordinates to the second conformer to the first conformer. Thus, two conformer would share the same coordinates. 
    """
    same_index = list()
    not_same = list()
    eps_2 = eps ** 2
    eps_modify_2 = eps_modify ** 2

    conf1 = ref.GetConformer(id=i1)
    conf2 = ref.GetConformer(id=i2)

    for i in index:
        tmp = conf1.GetAtomPosition(i)
        c1 = tmp.x, tmp.y, tmp.z

        tmp = conf2.GetAtomPosition(i)
        c2 = tmp.x, tmp.y, tmp.z

        if dist_2(c1, c2) > eps_2:
            not_same.append(i)
        elif eps_modify_2 < dist_2(c1, c2) <= eps_2:
            same_index.append(
                i
            )  # see them as same coordinate. But new conformation should be assign IDENTICAL coordinate.
            conf2.SetAtomPosition(i, conf1.GetAtomPosition(i))  # Assign
        else:
            same_index.append(i)

    # print("Same_index: %s" % same_index)
    if len(same_index) <= 1:
        return False, same_index
    else:
        return True, same_index


def f_AlignMolConformers(ref, atomIds, maxIters, reflect, RMSlist=None):
    """ Align conformers """
    index = atomIds
    rmsd = list()

    if (
        len(ref.GetConformers()) == 1
    ):  # if Only one conformation, then no alignation is needed.
        return [[0]]
    else:
        AlignMolConformers(
            ref, atomIds=atomIds, maxIters=500, reflect=True, RMSlist=rmsd
        )

        all_index = list()
        current = list()
        i = 0

        start = True
        while True:
            if i >= len(rmsd):
                break
            if start:
                current.append(0)
                start = False

            else:
                same_flag, same_index = is_same(ref, i - 1, i, index=index, eps=EPS)
                if same_flag:
                    current.append(i)
                    index = same_index
                else:
                    all_index.append(current)
                    current = list()
                    current.append(i)
                    index = atomIds

            i += 1

        all_index.append(current)

        return all_index

def f_AlignMolConformers_working2(ref, atomIds, maxIters, reflect, RMSlist=None):
    """ Align conformers """
    index = atomIds
    rmsd = list()

    if (
        len(ref.GetConformers()) == 1
    ):  # if Only one conformation, then no alignation is needed.
        return [[0]]
    else:
        AlignMolConformers(
            ref, atomIds=atomIds, maxIters=500, reflect=True, RMSlist=rmsd
        )

        all_families = list()
        all_indexes = list()

        # Go through all conformations
        for i in range(len(rmsd) + 1):
            flag_match = False

            # Go match all exist families ( But match most one )
            for current_family, current_index,tmpi in zip(all_families, all_indexes,range(len(all_families)) ):
                fam_ref = current_family[0]
                fam_index = current_index
                same_flag, same_index = is_same(ref, fam_ref, i, index=fam_index)
                if same_flag:  # Matched
                    flag_match = True
                    if len(same_index) == len(index):  # Matched and no change
                        all_families[tmpi].append(i)
                        pass
                    else:  # Match and change
                        all_indexes[tmpi] = same_index[:]
                        all_families[tmpi].append(i)
                    break
            if not flag_match:  # Do not match, Create new family
                all_families.append(
                    [i,]
                )
                all_indexes.append(atomIds[:])
            else:  # Do nothing
                pass

        return all_families

def f_AlignMolConformers_working3(ref, atomIds, maxIters, reflect, RMSlist=None):
    """ Align conformers """
    index = atomIds
    rmsd = list()

    if (
        len(ref.GetConformers()) == 1
    ):  # if Only one conformation, then no alignation is needed.
        return [[0]]
    else:
        AlignMolConformers(
            ref, atomIds=atomIds, maxIters=500, reflect=True, RMSlist=rmsd
        )

        all_families = list()
        all_indexes = list()

        # Go through all conformations
        for i in range(len(rmsd) + 1):
            flag_match = False

            # Go match all exist families ( But match most one )
            for current_family, current_index,tmpi in zip(all_families, all_indexes,range(len(all_families)) ):
                fam_ref = current_family[0]
                fam_index = current_index
                same_flag, same_index = is_same(ref, fam_ref, i, index=fam_index)
                if same_flag and len(same_index) == len(index):
                    flag_match = True
                    all_families[tmpi].append(i)
                    break
                # if same_flag:  # Matched
                #     flag_match = True
                #     if len(same_index) == len(index):  # Matched and no change
                #         all_families[tmpi].append(i)
                #         pass
                #     else:  # Match and change
                #         flag_match = False
                #         # all_indexes[tmpi] = same_index[:]
                #         # all_families[tmpi].append(i)
                #     break
            if not flag_match:  # Do not match, Create new family
                all_families.append(
                    [i,]
                )
                all_indexes.append(atomIds[:])
            else:  # Do nothing
                pass

        return all_families

# Core function
def f_AlignMolConformers_bak(ref, atomIds, maxIters, reflect, RMSlist=None):
    """ This function would not be used."""
    rmsd = list()

    if len(ref.GetConformers()) == 1:
        return [[0]]
    else:
        AlignMolConformers(
            ref, atomIds=atomIds, maxIters=500, reflect=True, RMSlist=rmsd
        )
        #         print("RMSD:")
        #         for tmp in rmsd:
        #             print(tmp)

        all_index = list()
        current = list()
        i = 0
        current_rmsd = 0
        start = True
        while True:
            if i >= len(rmsd):
                break
            if start:
                current_rmsd = 0
                current.append(0)
                start = False
            elif abs(current_rmsd - rmsd[i]) <= EPS:
                current.append(i)
            else:
                all_index.append(current)
                current = list()
                current.append(i)
                current_rmsd = rmsd[i]
            i += 1
        all_index.append(current)

        return all_index

# Read mol2 molecule once per time
def next_mol2_lines(infile):
    """Method to return one mol2 block once."""
    lines = list()

    for line in open(infile):
        if "@<TRIPOS>MOLECULE" in line:
            if len(lines) == 0:
                lines.append(line)
            else:
                yield lines
                lines = list()
                lines.append(line)
        else:
            lines.append(line)

    yield lines

# load all block
all_blocks = [x for x in next_mol2_lines(infile)]

# Staffs on first molecule
mol = Chem.MolFromMol2Block("".join(all_blocks[0]), removeHs=False)

# Get atom name to serial number mapping
if True:
    index_of_name = dict()
    for atom in mol.GetAtoms():
        name = atom.GetPropsAsDict()["_TriposAtomName"]
        sn = atom.GetIdx()
        index_of_name[name] = sn

# Get break down fragments
if True:
    # patt = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]')
    patt = Chem.MolFromSmarts("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]")

    # find the rotatable bonds
    bonds = mol.GetSubstructMatches(patt)

    # create an editable molecule, break the bonds, and add dummies:
    em = Chem.EditableMol(mol)
    nAts = mol.GetNumAtoms()
    for a, b in bonds:
        em.RemoveBond(a, b)
        em.AddAtom(Chem.Atom(0))
        em.AddBond(a, nAts, Chem.BondType.SINGLE)
        em.AddAtom(Chem.Atom(0))
        em.AddBond(b, nAts + 1, Chem.BondType.SINGLE)
        nAts += 2
    p = em.GetMol()
    Chem.SanitizeMol(p)

    frags = [x for x in Chem.GetMolFrags(p, asMols=True)]
    fragsindex = []
    for frag in frags:
        tmp = list()
        # print(len(frag.GetAtoms()))
        if len(frag.GetAtoms()) < 2:
            continue
        for atom in frag.GetAtoms():
            if atom.GetAtomicNum() > 1:
                name = atom.GetPropsAsDict()["_TriposAtomName"]
                tmp.append(index_of_name[name])
        if len(tmp) >= 2:
            fragsindex.append(tmp)
    print("fragsindex: %s" % fragsindex)

# Alignment molecules and output. # Second version.
if True:
    if matchversion == 'Highest':
        match_func = f_AlignMolConformers_working3
    else:
        pass
        # match_func = f_AlignMolConformers_working + matchversion 

    # Load reference molecule
    ref = Chem.MolFromMol2Block("".join(all_blocks[0]), removeHs=False)

    # load all other molecules
    all_mols = [
        Chem.MolFromMol2Block("".join(x), removeHs=False) for x in all_blocks[1:]
    ]

    # add all molecules to ref molecule's conformations
    for mol in all_mols:
        mol_conf = mol.GetConformer()
        ref.AddConformer(mol_conf, assignId=True)

    # Load pre-module
    from rdkit.Chem.rdMolAlign import AlignMolConformers

    # Aligning every frag
    i = 0
    os.system("mkdir sdf mol2 db2")
    for index in fragsindex:

        all_index = match_func(ref, atomIds=index, maxIters=500, reflect=True)
        print("Index:", index)
        print("All_index:", all_index)

        for conf_index in all_index:
            writer = Chem.SDWriter(f"sdf/output.{i}.sdf")
            for conf in conf_index:
                writer.write(ref, confId=conf)
            writer.close()
            os.system(
                f"{babel_bin} -isdf sdf/output.{i}.sdf -omol2 -O  mol2/output.{i}.mol2"
            )
            os.system(
                f"{python2_bin} {dockbase}/ligand/mol2db2/mol2db2.py -m mol2/output.{i}.mol2 -s output.solv -o db2/output.{i}.db2.gz"
            )
            i += 1

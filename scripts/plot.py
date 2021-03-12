#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys,os
import numpy as np
import math
from rdkit import Chem
from pathlib import Path
import argparse
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

parser = argparse.ArgumentParser("Analysis of dihedrals")
parser.add_argument('-i', '--input', required = True, nargs= '+', help="Input csv files")
parser.add_argument('-m','--mol2', help = 'Input mol2 molecule', default='')

args = parser.parse_args()

all_input_files = args.input
mol2 = args.mol2
# p = Path(".")
# all_input_files = list(p.glob("*.csv"))

dataset = dict()
for tmp in all_input_files:
    # name = tmp.name
    name = tmp
    tag = name[:-4]
    dataset[name] = tag

all_data = dict()
all_bonds = dict()

if mol2 :
    mol2 = Chem.MolFromMol2File(mol2)
    AllChem.Compute2DCoords(mol2)
    # Draw.MolToFile(mol2,'tmp.png')
    a = Draw.MolToMPL(mol2)
    for i in dir(a):
        print(i)
    for b in a.get_children():
        print(b)
    from matplotlib.backend_bases import RendererBase
    a.draw( renderer = RendererBase )
    a.show()
    sys.exit()

for filename,name in dataset.items():
#    print(name)
    data = []
    lines = open(filename).readlines()
    bonds = lines[0].strip().split()  # rotatable bonds list in first line
    all_bonds[name] = bonds

    for line in lines[1:]:
        data.append(line.split())
    data = np.array(data,dtype=np.float32)
    all_data[name] = data


ncate = len(all_data)
name = next(iter(all_data.keys()))

for i,bond in enumerate(all_bonds[name]):

    r = 1.0
    dr = 0.2

    outputfile = "%s.png"%(bond)
    plt.xticks( [] )
    plt.yticks( [] )
    plt.xlabel(" Dihedral distribution ")
    fig,ax = plt.subplots(figsize = (6,6))
    ax.set_xlim( [0-r-ncate*dr-0.2,r+ncate*dr+0.2] )
    ax.set_ylim( [0-r-ncate*dr-0.2,r+ncate*dr+0.2] )

    # color = colors.pop()
    colors = list(mcolors.TABLEAU_COLORS)

    for name,data in all_data.items():

        x = data[...,i]

        scale = np.random.random_sample(x.shape)
        tmpx = (scale * dr + r ) * np.cos( x*np.pi/180. )
        tmpy = (scale * dr + r ) * np.sin( x*np.pi/180. )

        ax.scatter(tmpx,tmpy,label = name,color = colors.pop() ,alpha = 0.3)

        r += dr

    plt.title( bond,fontsize ='xx-large' )
    ax.legend()
    plt.savefig(outputfile)
    plt.cla()


# for name,data in all_data.items():
#     for i,bond in enumerate(all_bonds[name]):
#         x = data[...,i]
# 
#         scale = np.random.random_sample(x.shape)
#         tmpx = (scale * dr + r ) * np.cos( x*np.pi/180. )
#         tmpy = (scale * dr + r ) * np.sin( x*np.pi/180. )
#         
#         outputfile = "%s.%s.png"%(name,bond)
#         
#         plt.figure(figsize=(6,6))
#         plt.xlim( [0-r-dr-0.2,r+dr+0.2] )
#         plt.ylim( [0-r-dr-0.2,r+dr+0.2] )
#         plt.xticks( [] )
#         plt.yticks( [] )
#         plt.xlabel(" Dihedra distribution ")
# 
#         plt.scatter(tmpx,tmpy,label = bond,color = 'g',alpha = 0.5)
#         plt.title( name,fontsize ='xx-large' )
#         plt.legend()
#         plt.savefig(outputfile)
#         plt.cla()



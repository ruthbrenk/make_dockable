#!/bin/bash

# Defination, Env
source ~/.bashrc
export PATH="$PATH:/home/soft/amber/19/bin/"
conda activate rdkit
export DOCKBASE='/pubhome/qyfu02/Software/db2_converter'
export OBABELBASE='/pubhome/qyfu02/miniconda3/envs/rdkit/bin'
export AMSOLEXE='/pubhome/qyfu02/Software/db2_converter/ligand/amsol/amsol7.1_patched'

CONF_EXE='/pubhome/qyfu02/Software/conformator/conformator'
PYTHON3='/usr/bin/python3.6'

# Read command line inputs
number=$1
max_conf=$2

# Make working dir ( also output directory, make sure this can be workable even run multiple times

mkdir $number
mv $number.smi $number

# Enter dir
cd $number

    # Generate conformations by Conformator
    $CONF_EXE -i $number.smi -o conformer.$number.mol2 -q 2 -n $max_conf --hydrogens

    # derive molecule smile and name from *.smi file
    zinc=`cat $number.smi|awk '{print $2}'`
    smile=`cat $number.smi|awk '{print $1}'`

    # derive first mol2 conformation for amsol
    $PYTHON3 $DOCKBASE/ligand/generate/derive_first_mol2.py conformer.$number.mol2 output$number.mol2

    # prapare.py ( Don't why, keep it )
    $PYTHON3 $DOCKBASE/ligand/generate/prepare.py output$number.mol2 --name=$zinc --smiles=$smile

    # AMSOL ( very complicated in it, but not much work )
    $DOCKBASE/ligand/amsol/calc_solvation.csh output$number.mol2 


    # Match and convert ( Avoid api of mol2db2.py but use match_and_convert_mol2.py instead.
    rm -rf db2
    $PYTHON3 $DOCKBASE/ligand/mol2db2/match_and_convert_mol2.py conformer.$number.mol2 output

    # clean
    rm -r mol2 sdf 

    # Merge multiple db2 fils into one
    cd db2
        zcat -f *.db2.gz > all.db2
        rm *.gz
        gzip all.db2
        mv all.db2.gz ../
    cd ..
    rmdir db2

# Back
cd ..



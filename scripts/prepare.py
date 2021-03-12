#!/usr/bin/env python
import sys,os
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='Generate db2 for smile library')
parser.add_argument('-i', type=str, required = True, help="Input file *.smi. Make sure there are compound names there")
parser.add_argument('-o', type=str, default = 'DB2_Generating', help = 'Working dir' )
parser.add_argument('-s', type=int, default = 500, help = 'Max number of compounds in each directory' ) 

args = parser.parse_args()

infile = args.i
outdir = args.o
N = args.s

p = Path('.')

work_dir = p/outdir
work_dir.mkdir(exist_ok = True)

n_split = 0
count = 0
split_dir = None
for line in open( infile ):
    if len(line.split()) < 2:
        sys.stderr.write("##### Error. There's no atom name for compound. \n")
        sys.stderr.write("##### That a required option. Otherwise, db2 file can not trace back to origin molecule\n")
        sys.stderr.write("##### Error line is:\n")
        sys.stderr.write(line)
        sys.exit(1)
    smile = line.split()[0]
    name = line.strip().split()[-1]
    newline = f"{smile}    {name}\n"
    if count % N == 0:
        n_split += 1
        split_dir = work_dir/f'split_{n_split}'
        split_dir.mkdir(exist_ok = True)
    db2_dir = split_dir/f'{count}'
    db2_dir.mkdir( exist_ok = True )
    db2_smi = db2_dir/f'{count}.smi'
    db2_smi.write_text( newline )

    count += 1
        


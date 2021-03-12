#!/usr/bin/env python
import sys,os

infile = sys.argv[1]
outfile = sys.argv[2]

def next_mol2_lines(infile):
    '''Method to return one mol2 block once.'''
    lines = list()

    for line in open(infile):
        if '@<TRIPOS>MOLECULE' in line :
            if len(lines) == 0:
                lines.append(line)
            else:
                yield lines
                lines = list()
                lines.append(line)
        else:
            lines.append(line)

    yield lines

all_blocks = [x for x in next_mol2_lines(infile)]


with open(outfile,'w') as ofp:
    for block in all_blocks:
        line2 = block[2].strip()
        newline2 = f"{line2}  0  0  0\n"

        block.pop(2)
        block.insert(2,newline2)
        block.insert(5,'\n')
        block.insert(6,'mmff94s = 1515.15\n')

        if block[-1].strip() == '':
            block.pop()

        for line in block:
            ofp.write(line)




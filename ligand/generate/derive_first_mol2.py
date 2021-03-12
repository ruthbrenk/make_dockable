#!/usr/bin/env python3.6
import sys,os

infile = sys.argv[1]
outfile = sys.argv[2]

with open(outfile,'w') as ofp:
    start = True
    for line in open(infile):
        if "@<TRIPOS>MOLECULE" in line and start:
            start = False
            ofp.write(line)
        elif "@<TRIPOS>MOLECULE" in line and not start:
            break
        else:
            ofp.write(line)




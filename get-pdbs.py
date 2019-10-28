#!/usr/bin/env python3
import Bio
from Bio.PDB import PDBList
'''Selecting structures from PDB'''
pdb1 = PDBList()
PDBlist2=[]

#this will ensure that the script can be run for any sublist without having
#to modify the actual script.
print("Which sublist of structures to analyze? (100 structures per list, sublist 1-425)")
num = input()

with open('../sublists/sublist'+num+'.txt','r') as infile:
    for line in infile:
        PDBlist2.append(line.strip())

for i in PDBlist2:
    pdb1.retrieve_pdb_file(i,pdir='../PDB',file_format='pdb')

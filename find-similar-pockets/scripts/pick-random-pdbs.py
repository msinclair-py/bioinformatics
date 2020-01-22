#!/usr/bin/env python
import random, Bio
from Bio.PDB import PDBList

pdb1 = PDBList()
PDB2list=[]

array=[]
with open('../sublists/pdb_list.txt', 'r') as infile:
	for line in infile:
		array.append(line)

samp = random.sample(range(0, len(array)), 100)
for num in samp:
	PDB2list.append(array[num].strip())
for i in PDB2list:
	pdb1.retrieve_pdb_file(i,pdir='../PDB',file_format='pdb')

#!/usr/bin/env python
import os, glob

#This block finds all the directories remaining and 
#puts them in an array
names=[]
for name in glob.glob('../PDB/*_out'):
	pdb =  name.split('/')
	names.append(pdb[2][:-4])

#This block looks for whether the pdb file has a matching
#directory or not. PDB's without pockets should not have
#directories and should be deleted by this script
for name in glob.glob('../PDB/*.ent'):
	pdb = name.split('/')
	if pdb[2][:-4] not in names:
		os.remove(name)
	else:
		continue

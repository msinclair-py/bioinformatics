#!/usr/bin/env python
import glob, os, shutil

script='/home/matt/Software/rosetta_src_2019.35.60890_bundle/tools/protein_tools/scripts/clean_pdb.py'

for i in glob.glob('../PDB/*'):
	inp = str(script)+' '+str(i)+' A'
	os.system(inp)

for i in glob.glob('*.fasta'):
	name = i.split('.')[0]+'.'+i.split('.')[1]+'.pdb'
	os.remove(i)
	shutil.move(name,'../PDB/')

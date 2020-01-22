#!/usr/bin/env python
import os, glob, subprocess

alphas = 80

fpocket = '/home/matt/Software/fpocket/bin/fpocket'

for i in glob.glob('../PDB/*.pdb'):
	args=(fpocket,'-f',i,'-i',alphas)
	str_args = [ str(x) for x in args ]
	subprocess.run(str_args)

for i in glob.glob('../PDB/*_out'):
	if len(os.listdir(i)) == 0:
		print(f'{i}: empty directory')
	else:
		print(f'{i}: has pockets')

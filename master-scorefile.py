#!/usr/bin/env python
import glob, os.path
import numpy as np

scores = np.genfromtxt('../score.txt', dtype='unicode', skip_header=1)

if os.path.exists('../winners/master_scorefile.txt'):
	mode='a'
else:
	mode='w'
locations=['1m1n', '2']

with open('../winners/master_scorefile.txt', mode) as f:
	if mode == 'w':
		f.write('PDB   Pock   Target Vol  Pock Vol    Int Vol   Int %  #hits  winners directory\n')
	for s in scores:
		for l in locations:
			if s[0] in l:
				out=f'{s[0]:<6}{s[1]:<7}{s[2]:{12}.{8}}{s[3]:{12}.{8}}{s[4]:{11}.{7}}{s[5]:{7}.{3}}{s[6]:<7}{l[0]}\n'
				f.write(out)

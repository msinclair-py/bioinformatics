#!/usr/bin/env python
import os, glob, subprocess

#need our aliases
surf='/home/matt/Software/VASPreleasePack/codeRelease/surfProcessingRelease/debug/surfProcessing'
vasp='/home/matt/Software/VASPreleasePack/codeRelease/vaspRelease/debug/vasp'


#check to see if VASP scores and VASP pockets directory exists
if not os.path.exists('../VASP/scores'):
   os.mkdir('../VASP/scores')
if not os.path.exists('../VASP/pockets'):
   os.mkdir('../VASP/pockets')

#get surf file for each structure
for pock in glob.glob('../PDB/aligned_*'):
	pocket = pock.split('/')[2]
	n = pocket.split('_')[1]
	p = pocket.split('_')[2]
	arg = (surf,'-surfProbeGen',pock,f'../VASP/pockets/{n}_{p}.SURF',3,.5)
	str_arg = [ str(x) for x in arg ]
	subprocess.run(str_arg)


#run intersect VASP on each structure, and get the original pocket volume
for sur in glob.glob('../VASP/pockets/*.SURF'):
	pocket = sur.split('/')[3]
	n = pocket.split('.')[0]  #NOTE: the format of this will be 1cpe_pock1 for example
	for init in glob.glob('../initial_structures/*.SURF'):
		info = init.split('/')[2]
		conf = info.split('.')[1]
		rot = info.split('.')[2]
		tilt = info.split('.')[3]
		trans = info.split('.')[4]
		flip = info.split('.')[5]
		arg = (vasp,'-csg',sur,init,'I',f'../VASP/intersect.{n}.{conf}.{rot}.{tilt}.{trans}.{flip}.SURF',.5)
		str_arg = [ str(x) for x in arg ]
		subprocess.run(str_arg)
	
	arg2 = (surf,'-surveyVolume',sur)
	str_arg2 = [ str(x) for x in arg2 ]
	out = open('../VASP/pockets/{n}.vol.txt','w')
	subprocess.run(str_arg2, stdout=out)

#get the score of each structure generated
for inter in glob.glob('../VASP/intersect*'):
	pocket = inter.split('/')[2]
	na = pocket.split('.')[1]
	co = pocket.split('.')[2]
	ro = pocket.split('.')[3]
	ti = pocket.split('.')[4]
	tr = pocket.split('.')[5]
	fl = pocket.split('.')[6]
	arg = (surf,'-surveyVolume',inter)
	out = open(f'../VASP/scores/{na}_{co}_{ro}_{ti}_{tr}_{fl}_iscore.txt', 'w')
	str_arg = [ str(x) for x in arg ]
	subprocess.run(str_arg, stdout=out)

#get the original volume of each pocket
for v in glob.glob('../VASP/pockets/*.SURF'):
	pocket = v.split('/')[3]
	n = pocket.split('.')[0]
	arg = (surf,'-surveyVolume',v)
	out = open(f'../VASP/pockets/{n}.vol.txt', 'w')
	str_arg = [ str(x) for x in arg ]
	subprocess.run(str_arg, stdout=out)

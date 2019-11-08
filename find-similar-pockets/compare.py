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
for pock in glob.glob('../PDB/aligned.*'):
	pocket = pock.split('/')[2]
	n = pocket.split('_')[1]
	p = pocket.split('_')[2]
	arg = (surf,'-surfProbeGen',pock,'../VASP/pockets/'+n+'_'+p+'.SURF','3','.5')
	str_arg = [ str(x) for x in arg ]
	subprocess.run(str_arg)


#run intersect VASP on each structure
for sur in glob.glob('../VASP/pockets/*.SURF'):
	pocket = sur.split('/')[3]
	n = pocket.split('.')[0]
	for init in glob.glob('../initial_structures/*.SURF'):
		rot = init.split('/')[2][:-5]
		arg = (vasp,'-csg',sur,init,'I','../VASP/'+n+'.intersect.'+rot+'.SURF','.5')
		str_arg = [ str(x) for x in arg ]
		subprocess.run(str_arg)

#get the score of each structure generated
for inter in glob.glob('../VASP/*intersect*'):
	pocket = inter.split('/')[2]
	n = pocket.split('.')[0]
	r = pocket.split('.')[3]+'.'+pocket.split('.')[4]
	arg = (surf,'-surveyVolume',inter)
	out = open(f'../VASP/scores/{n}_{r}_iscore.txt', 'w')
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

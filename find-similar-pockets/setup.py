#!/usr/bin/env python
import glob, subprocess, math
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

#this script is predicated on you already having the desired
#pocket chosen. you can input the pocket number for excision
#and further analysis. this is the same as the chain number
#in the pdb e.g. STP 14 would correspond to pocket 14
#NOTE: if you run fpocket with the -i 80 flag the target pocket
#for 1m1n is any of pockets 1-4 due to being a homotetramer
print('What is the four letter code of your target pdb?')
pdb = input()
print('What is the number of your pocket?')
print('(On the PDB it will appear as STP ##)')
pocknumber = input()
print('How many degrees would you like to rotate about z?')
deg1 = input()
print('How many degrees would you like to tilt by?')
deg2 = input()

#here are just some housekeeping variables, functions and 
#initialized arrays to be used in the following script
vmd='usr/local/bin/vmd'
surf='/home/matt/Software/VASPreleasePack/codeRelease/surfProcessingRelease/debug/surfProcessing'
rotates = int(360/float(deg1))
pocket=[]
vec_residues=[]

#function to find the residues with maximum separation. these will be needed
#to align to the z axis. returns longest distance residue pair by res id
def initialize (array):
	D=pdist(array)
	D=squareform(D)
	temp = np.where(D == D.max())
	return temp[0]


#orient the protein along the z axis, add 1 to the output of temp because
#while both temp and residues array index from 0, there is a header in pdb
def orient (array, residues):
	fpath = f'../initial_structures/cent_{pdb}.pocket.pdb'
	res1 = residues[0]+1
	res2 = residues[1]+1
	with open(fpath, 'r') as infile:
		for i, line in enumerate(infile):
			if i == res1:
				coord1 = line.split()[6:9]
				c1 = [float(i) for i in coord1]
			elif i == res2:
				coord2 = line.split()[6:9]
				c2 = [float(i) for i in coord2]
	
	#get the normalized vector of the longest dist pair
	d = math.sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2+(c1[2]-c2[2])**2)
	norm = [(c2[0]-c1[0])/d, (c2[1]-c1[1])/d, (c2[2]-c1[2])/d]

	#pass the normalized vector to vmd to do alignment
	na = str(fpath.split('/')[2])
	a1,a2,a3 = str(norm[0]),str(norm[1]),str(norm[2])
	arguments = ['vmd','-dispdev','text','-e','align.tcl','-args',na,a1,a2,a3,1]
	string_args = [ str(x) for x in arguments ]
	subprocess.run(string_args)


###body of code starts here; get the pocket from fpocket output
with open('../initial_structures/'+pdb+'_out/'+pdb+'_out.pdb',
	'r') as infile:
	for line in infile:
		if "STP" in line:
			if line.split()[5] == pocknumber:
				pocket.append(line)
			else:
				continue
		else:
			continue

#create new pocket pdb structure
with open('../initial_structures/'+pdb+'.pocket.pdb', 'w') as out:
	[out.write(line) for line in pocket]

#pass details to vmd to center the pocket before alignment
args = ['vmd','-dispdev','text','-e','center.tcl','-args',1,pdb]
str_args = [ str(x) for x in args ]
subprocess.run(str_args)
print('after center.tcl')
#acquire coords of pocket in numpy matrix
coords = np.genfromtxt(f'../initial_structures/cent_{pdb}.pocket.pdb',
	dtype=float, skip_header=1, skip_footer=1, usecols=(6,7,8))
print('got coords')
#get longest axis and align
vec_residues=initialize(coords)
print('before orient')
orient(coords, vec_residues)
print('after orient')

#acquire aligned coords
aligned_coords = np.genfromtxt(f'../initial_structures/aligned.cent_{pdb}.pocket.pdb',
	dtype=float, skip_header=1, skip_footer=1, usecols=(6,7,8))

#do rotations in vmd because its very fast and easy
args = ['vmd','-dispdev','text','-e','rotater.tcl','-args',pdb,deg1,deg2]
str_args = [str(x) for x in args ]
subprocess.run(str_args)

#perform SURF on all the various nitrogenase pockets
for i in glob.glob('../initial_structures/*tilt*'):
	stuff=i.split('/')[2]
	rotation = stuff.split('.')[1]
	tstat = stuff.split('.')[2]
	args = (surf,'-surfProbeGen',i,'../initial_structures/'+pdb+'.'+rotation+'.'+tstat+'.SURF',3,.5)
	str_args = [ str(x) for x in args ]
	subprocess.run(str_args)

#get the volume of the intitial pocket for reference in scorefile
vargs = (surf,'-surveyVolume','../initial_structures/'+pdb+'.0.notilt.SURF')
str_vargs = [ str(x) for x in vargs ]
out = open('../initial_structures/vol.txt', 'w')
subprocess.run(str_vargs, stdout=out)

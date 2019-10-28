#!/usr/bin/env python
import glob, subprocess, math
import numpy as np
from scipy.spatial.distance import pdist, squareform

#get the residue pair corresponding to the longest axis through
#the pocket and return it
def initialize (array, name):
	D=pdist(array)
	D=squareform(D)
	temp = np.where(D == D.max())
	return temp[0]

#acquire names of pockets and coords of pocket in separate numpy
#matrices, name_array and coords respectively
name_array=[]
for name in glob.glob('../PDB/cent_*.pdb'):
	name_array.append(name.split('/')[2])

coords=[]
pair=[]
for entry in name_array:
	coords = np.genfromtxt('../PDB/'+entry, dtype=float, skip_header=1,
		skip_footer=1, usecols=(6,7,8))

#run all the functions
	pair.append(initialize(coords, entry))

#the ith entry in both the pair and name_array lists correspond
#to each other. res1 and res2 are stored as the pair ID + 1 due
#to the header present in the input pdb file. both the pdb and
#the pair list will index from 0 so this doesn't need to be taken
#into account
for i in range(len(pair)):
	res1 = pair[i][0]+1
	res2 = pair[i][1]+1
	with open('../PDB/'+name_array[i], 'r') as infile:
		for j, line in enumerate(infile):
			if j == res1:
				coord1 = line.split()[6:9]
				c1 = [float(k) for k in coord1]
			elif j == res2:
				coord2 = line.split()[6:9]
				c2 = [float(k) for k in coord2]
	
	#to get the vector normal calculate the net distance and average by the
	#scalar distance of the points
	d = math.sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2+(c1[2]-c2[2])**2)
	norm = [(c2[0]-c1[0])/d, (c2[1]-c1[1])/d, (c2[2]-c1[2])/d]

	#call the vmd script to do the aligning by using a fancy bash argument
	#assignment variable and the subprocess utility. NOTE: this utility uses
	#subprocess.run rather than subprocess.call as of ~python 3.5
	na = '../PDB/'+str(name_array[i])
	a1,a2,a3 = str(norm[0]),str(norm[1]),str(norm[2])
	args = ['vmd','-dispdev','text','-e','align.tcl','-args',na,a1,a2,a3,2]
	str_args = [ str(x) for x in args ]
	subprocess.run(str_args)

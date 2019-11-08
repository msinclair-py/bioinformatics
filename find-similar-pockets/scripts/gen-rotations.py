#!/usr/bin/env python
import math, subprocess, glob
import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

#need this location since bash aliases don't work in python
vmd='usr/local/bin/vmd'

#obtain the name of the protein
for struc in glob.glob('../initial_structures/*'):
	line = struc.split('/')[2]
	if 'cent' == line[:4]:
		a = line.split('.')[0]
		name = a[-4:]

#pdb making function
def make_pdb (array, deg, tilt):
	with open(f'../initial_structures/{name}_{deg}{tilt}', 'w') as f:
		with open(f'../initial_structures/cent_{name}.pocket.pdb.aligned', 'r') as infile:
			t=-1
			for line in infile:
				if "CRYST1" in line:
					f.write(line)
				elif "END" in line:
					f.write(line)
				else:
					x="%.3f" % array[t][0]
					y="%.3f" % array[t][1]
					z="%.3f" % array[t][2]
					f.write(line[:32]+str(x)+'  '+str(y)+'  '+str(z)+line[55:])
				t+=1

#find residues to align along z axis. this essentially calculates the distance
#between every combination of points and returns the longest distance pair
def initialize (array):
	D=pdist(array)
	D=squareform(D)
	temp = np.where(D == D.max())
	return temp[0]

#orient the protein along the z axis, add 1 to the output of temp
#because while both temp and residues array index from 0, there is 
#a header in the pdb file
def orient (array, residues):
	fpath = f'../initial_structures/cent_{name}.pocket.pdb'
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

#need distance to normalize the longest vector
	d = math.sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2+(c1[2]-c2[2])**2)
	norm = [(c2[0]-c1[0])/d, (c2[1]-c1[1])/d, (c2[2]-c1[2])/d]

#pass the normalized vector to vmd align.tcl script to do aligning	
	na = str(fpath.split('/')[2])
	a1,a2,a3 = str(norm[0]),str(norm[1]),str(norm[2])
	args = ['vmd','-dispdev','text','-e','align.tcl','-args',na,a1,a2,a3,1]
	str_args = [ str(x) for x in args ]
	subprocess.run(str_args)


#rotation function
def rotate (array, tilt):
	if tilt == 1:
		#this will tilt the pocket positively along x
		tilt = "tilted"
		tilt_pos = R.from_euler('x', 30, degrees=True)
		for i in range(len(array)):
			array[i] = tilt_pos.apply(array[i])
	else:
		tilt = "notilt"
		
	#z rotation scipy method
	z_rot = R.from_euler('z', 15, degrees=True)	
	#rotate about z axis, all tilts undergo this
	for i in range(24):
		for j in range(len(array)):
			array[i] = z_rot.apply(array[i])
		deg = i*15
		make_pdb(array, deg, tilt)


#acquire coords of pocket in numpy matrix
coords = np.genfromtxt(f'../initial_structures/cent_{name}.pocket.pdb', dtype=float, skip_header=1,
	skip_footer=1, usecols=(6,7,8))

#run all the functions
vec_residues=[]
vec_residues=initialize(coords)
orient(coords, vec_residues)

#acquire the now aligned coords
aligned_coords = np.genfromtxt(f'../initial_structures/cent_{name}.pocket.pdb.aligned', dtype=float,
	skip_header=1, skip_footer=1, usecols=(6,7,8))

#call the rotation functions
rotate(aligned_coords, 0)
rotate(aligned_coords, 1)

#!/usr/bin/env python
import glob, subprocess, math
import numpy as np
from scipy.spatial.transform import Rotation as R
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
rotates = int(360/float(deg1))
pocket=[]
vec_residues=[]

#pdb making function for rotated pockets
def make_pdb (array, deg, tilt):
	with open(f'../initial_structures/{pdb}_{deg}_{tilt}', 'w') as f:
		with open(f'../initial_structures/cent_{pdb}.pocket.pdb.aligned', 'r') as infile:
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


#rotation function
def rotate (array, tilt):
	#z rotation and tilt scipy method
	z_rot = R.from_euler('z', deg1, degrees=True)
	tilt_pos = R.from_euler('x', deg2, degrees=True)

	#tilt or don't
	if tilt == 1:
		tilt = "tilted"
		for i in range(len(array)):
			array[i] = tilt_pos.apply(array[i])
	else:
		tilt = "notilt"

	#rotate about z axis
	for i in range(rotates):
		for j in range(len(array)):
			array[i] = z_rot.apply(array[i])
		deg = i*15
		make_pdb(array, deg, tilt)


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

#acquire coords of pocket in numpy matrix
coords = np.genfromtxt(f'../initial_structures/cent_{pdb}.pocket.pdb',
	dtype=float, skip_header=1, skip_footer=1, usecols=(6,7,8))

#get longest axis and align
vec_residues=initialize(coords)
orient(coords, vec_residues)

#acquire aligned coords
aligned_coords = np.genfromtxt(f'../initial_structures/cent_{pdb}.pocket.pdb.aligned',
	dtype=float, skip_header=1, skip_footer=1, usecols=(6,7,8))

#run rotations
rotate(aligned_coords, 0)
rotate(aligned_coords, 1)

#perform SURF on all the various nitrogenase pockets
subprocess.run('./pysurf')

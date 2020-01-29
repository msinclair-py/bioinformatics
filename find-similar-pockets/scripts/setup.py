#!/usr/bin/env python
import glob, subprocess, math
import numpy as np
import scipy
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
pdb = '1m1n' #input()
print('What is the number of your pocket?')
print('(On the PDB it will appear as STP ##)')
pocknumber = '1' #input()
print('How many degrees would you like to rotate about z?')
deg1 = '30' #input()
print('How many degrees would you like to tilt by?')
deg2 = '20' #input()
print('How many angstroms would you like to translate by?')
tr = 2 #input()

#here are just some housekeeping variables, functions and 
#initialized arrays to be used in the following script
surf='/home/matt/Software/VASPreleasePack/codeRelease/surfProcessingRelease/debug/surfProcessing'
rotates = int(360/float(deg1))
conformations = [-30,0,30,150,180,210]
pocket=[]
vec_residues=[]


#pdb making function for rotated pockets
def make_pdb (array,conf,deg,tilt,trans,flipped):
	with open(f'../initial_structures/{pdb}_{conf}_{deg}_{tilt}_{trans}_{flipped}', 'w') as f:
		with open(f'../initial_structures/{pdb}.pocket.pdb', 'r') as infile:
			t=0
			for line in infile:
				x=f'{array[t][0]:.3f}'
				y=f'{array[t][1]:.3f}'
				z=f'{array[t][2]:.3f}'
				f.write(f'{line[:27]}{x:>11}{y:>8}{z:>8} {line[55:]}')
				t+=1


#function to find the residues with maximum separation. these will be needed
#to align to the z axis. returns longest distance residue pair by res id
def initialize (array):
	D=pdist(array)
	D=squareform(D)
	temp = np.where(D == D.max())
	return temp[0]


#centers the coordinate array
def center(array):
	centered_array=np.zeros((array.shape[0],array.shape[1]))
	com=array.mean(0)
	for i in range(array.shape[0]):
		centered_array[i]=array[i]-com
	
	return centered_array


#aligns the longest axis to the z axis
def align(array,residues):
	aligned = np.zeros((array.shape[0],array.shape[1]))
	res1 = residues[0]
	res2 = residues[1]
	c1 = array[res1]
	c2 = array[res2]
	b = np.array([0,0,1])

	d = math.sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2+(c1[2]-c2[2])**2)
	norm = np.array([(c2[0]-c1[0])/d, (c2[1]-c1[1])/d, (c2[2]-c1[2])/d])

	#compute rotation angles theta, psi and phi by using trig on 2d projections
	xyprojection = np.array([norm[0], norm[1], 0])
	xzprojection = np.array([norm[0], 0, norm[2]])
	yzprojection = np.array([0, norm[1], norm[2]])
	theta = np.arccos(np.dot(xyprojection,[0,1,0]))
	psi = np.arccos(np.dot(xzprojection,[1,0,0]))
	phi = np.arccos(np.dot(yzprojection,[0,0,1]))
	rot_theta = R.from_euler('z',theta,degrees=False)
	rot_psi = R.from_euler('y',psi,degrees=False)
	rot_phi = R.from_euler('x',phi,degrees=False)

	al1 = rot_theta.apply(array)
	al2 = rot_psi.apply(al1)
	aligned = rot_phi.apply(al2)

	return aligned


#conformation generation/flow control function
def conformation(array):
	for i in range(len(conformations)):
		deg = conformations[i]
		conform = R.from_euler('z',deg,degrees=True)
		conformation = conform.apply(array)
		make_pdb(conformation,f'conf{i}',0,0,0,0)
		tilt(conformation,f'conf{i}',0,0,0)
		flip(conformation,f'conf{i}',0)
		translate(conformation,f'conf{i}')


#rotation function
def rotate(array,conf,tilted,trans,flipped):
	rotated_array=np.array([])

	for i in range(rotates):
		deg = i*int(deg1)
		z_rot = R.from_euler('z',deg,degrees=True)
		rotated_array = z_rot.apply(array)
		make_pdb(rotated_array,conf,deg,tilted,trans,flipped)

	if tilted == 1:
		tilt(array,conf,tilted,trans,flipped)


#tilt function
def tilt(array,conf,tlt,trans,flipped):
	tlt+=1
	tilted_array=np.array([])
	tilter = R.from_euler('x', deg2, degrees=True)
	tilted_array = tilter.apply(array)
	rotate(tilted_array,conf,tlt,trans,flipped)


#translation function
def translate(arr,conf):
	xpos=np.zeros((arr.shape[0],arr.shape[1]))
	xneg=np.zeros((arr.shape[0],arr.shape[1]))
	ypos=np.zeros((arr.shape[0],arr.shape[1]))
	yneg=np.zeros((arr.shape[0],arr.shape[1]))
	zpos=np.zeros((arr.shape[0],arr.shape[1]))
	zneg=np.zeros((arr.shape[0],arr.shape[1]))
	x=np.array([tr,0,0])
	y=np.array([0,tr,0])
	z=np.array([0,0,tr])

	for i in range(len(arr)):
		xpos[i]=arr[i]+x
		xneg[i]=arr[i]-x
		ypos[i]=arr[i]+y
		yneg[i]=arr[i]-y
		zpos[i]=arr[i]+z
		zneg[i]=arr[i]-z

	#iterate through list of translations, performing all previous
	#conformational changes
	translations = [xpos,xneg,ypos,yneg,zpos,zneg]
	for i in range(1,len(translations)):
		make_pdb(translations[i],conf,0,0,i,0)
		tilt(translations[i],conf,0,i,0)
		flip(translations[i],conf,i)


#180 degree flip function
def flip(arr,conf,trans):
	flp=np.array([])
	f = R.from_euler('x', 180, degrees=True)
	flp = f.apply(arr)
	make_pdb(flp,conf,0,0,trans,1)
	tilt(flp,conf,0,trans,1)


###body of code starts here; get the pocket from fpocket output
with open(f'../initial_structures/{pdb}_out/{pdb}_out.pdb',
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
with open(f'../initial_structures/{pdb}.pocket.pdb', 'w') as out:
	[out.write(line) for line in pocket]

initcoords = np.genfromtxt(f'../initial_structures/{pdb}.pocket.pdb',
	dtype=float, usecols=(6,7,8))

cent = center(initcoords)
vec_residues = initialize(cent)
aligned_coords = align(cent, vec_residues)

#run rotations, the conformation function orchestrates the order of operations
conformation(aligned_coords)

#pass each aligned structure through SURF to generate the pocket surface
for i in glob.glob('../initial_structures/*conf*'):
	stuff=i.split('/')[2]
	conf=stuff.split('_')[1]
	rotation=stuff.split('_')[2]
	tilt=stuff.split('_')[3]
	tran=stuff.split('_')[4]
	flip=stuff.split('_')[5]
	args=(surf,'-surfProbeGen',i,f'../initial_structures/{pdb}.{conf}.{rotation}.{tilt}.{tran}.{flip}.SURF',3,.5)
	str_args = [ str(x) for x in args ]
	subprocess.run(str_args)

#get the volume of the initial pocket for reference in the scorefile
vargs = (surf,'-surveyVolume',f'../initial_structures/{pdb}.conf0.0.0.0.0.SURF')
str_vargs = [ str(x) for x in vargs ]
out = open('../initial_structures/vol.txt','w')
subprocess.run(str_vargs, stdout=out)

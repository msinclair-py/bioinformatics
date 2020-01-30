#!/usr/bin/env python
import glob, subprocess, math
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.spatial.transform import Rotation as R

#write out new pdbs after alignment
def make_pdb(array,name):
	with open(f'../PDB/aligned_{name[:-4]}','w') as f:
		with open(f'../PDB/{name}','r') as infile:
			i=0
			for line in infile:
				x = f'{array[i][0]:.3f}'
				y = f'{array[i][1]:.3f}'
				z = f'{array[i][2]:.3f}'
				f.write(f'{line[:27]}{x:>11}{y:>8}{z:>8} {line[55:]}')
				i+=1

#get the residue pair corresponding to the longest axis through
#the pocket and return it
def initialize (array):#, name):
	D=pdist(array)
	D=squareform(D)
	temp = np.where(D == D.max())
	return temp[0]

def principal (array):  #want axis 1
	inertia = array.T @ array
	e_values, e_vectors = np.linalg.eig(inertia)
	order = np.argsort(e_values)
	eval3,eval2,eval1 = e_values[order]
	axis3,axis2,axis1 = e_vectors[:,order].T
	return axis1

#find center of mass and translate coordinate array to be centered
#about the origin
def center(array):
	centered_array=np.zeros((array.shape[0],array.shape[1]))
	com=array.mean(0)
	for i in range(array.shape[0]):
		centered_array[i]=array[i]-com
	return centered_array


def align(array,a):
	aligned=np.zeros((array.shape[0],array.shape[1]))
	b = [0,0,1]
	v = np.cross(a,b)
	c = np.dot(a,b)
	I = np.eye(3,3)
	vx = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
	R = I + vx + (vx @ vx)/(1+c)
	aligned = R @ array.T
	return aligned.T

#acquire names of pockets and coords of pocket in separate numpy
#matrices, name_array and coords respectively
name_array=[]
for name in glob.glob('../PDB/*_pock*'):
	name_array.append(name.split('/')[2])

for entry in name_array:
	coords = np.genfromtxt(f'../PDB/{entry}',dtype=float,usecols=(6,7,8))

#run all the functions
	centered=center(coords)
#	pair=initialize(coords)
	vector=principal(coords)
	aligned=align(centered,vector)#pair)
	make_pdb(aligned,entry)

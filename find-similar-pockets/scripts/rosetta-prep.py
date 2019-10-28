#!/usr/bin/env python
import numpy as np

#filt is the percentage filter cutoff and hilt is the hit filter
filt = .7
hilt = 0

lst = np.genfromtxt('../score.txt', dtype='unicode', skip_header=1,
	usecols=(0,1,5,6))

#using relevant score information including pdb id, pocket number, intersect 
#percentage and number of orientation hits for the filter, get the unique
#residue numbers into a pose file for rosetta
for i in range(len(lst)):
	array=[]
	pose_array=[]
	if float(lst[i][2]) > float(filt):
		if float(lst[i][3]) > float(hilt):
			a=lst[i][0]
			b=lst[i][1][-1]
			fil = f'../PDB/pdb{a}_out/pockets/pocket{b}_atm.pdb'
			pose_array = np.genfromtxt(fil, dtype='unicode', skip_header=20,
				skip_footer=2, usecols=5)
			
			#go through pose_array and take only the unique entries
			for i in range(len(pose_array)):
				if pose_array[i] not in array:
					if "." in pose_array[i]:
						continue
					else:
						array.append(pose_array[i])

			#convert array to int type since it is currently in strings, it can't be
			#sorted properly if it is full of strings
			array = [int(x) for x in array]

			#output the sorted array into space delimited pose file
			with open(f'../rosetta/{a}_pock{b}.pos', 'w') as outfile:
				for line in sorted(array):
					if line == sorted(array)[-1]:
						outfile.write(str(line))
					else:
						outfile.write(str(line)+' ')

#!/usr/bin/env python
import glob

###cutoff for maximum number of alpha spheres in pocket
cutoff=150

###function to evaluate number of unique pockets found
def unique(list1):
	unique_list=[]
	for x in list1:
		if x not in unique_list:
			unique_list.append(x)
	return unique_list

###get the names of the target pdbs
name_array=[]
for name in glob.glob('../PDB/*_out/*_out.pdb'):
	name_array.append(name.split('/')[3])

###pull out all the pocket coords and the pocket id numbers
for entry in name_array:
	pock_array=[]
	pocket_id=[]
	uniq=[]
	with open('../PDB/'+entry[:-4]+'/'+entry, 'r') as infile:
		for line in infile:
			if not "STP" in line:
				continue
			elif "STP" in line:
				pock_array.append(line)
				pocket_id.append(line.split()[5])

	###now determine the number of output files necessary for each
	###pdb by finding unique pocket ids and output that many new pocket
	###pdb files
	uniq = unique(pocket_id)
	for i in range(0, len(uniq)):
		j=i+1
		a = pock_array.count(j)
		if a < cutoff:
			with open(f'../PDB/{entry.split("_")[0]}_pock{j}.pdb', 'w') as outfile:
				for line in pock_array:
					if line.split()[5] == uniq[i]:
						outfile.write(line)
					else:
						continue
		else:
			continue

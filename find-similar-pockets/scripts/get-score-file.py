#!/usr/bin/env python
import glob, os
import numpy as np

#filter for percent of intersect volume from vasp
#modify this here as it gets called in the script below
filt=.7

#get list of names for corresponding scores
name_array=[]
raw_name=[]
pdbs=[]
inter=[]
pockets=[]
for name in glob.glob('../VASP/scores/*iscore.txt'):
	item = name.split('/')[3]
	raw_name.append(item)
	if item not in pdbs:
		pdbs.append(item.split('_')[0])

for name in raw_name:
	name_array.append(name.split('_'))

print('-----Getting scores-----')
#obtain all scores as well as names into a new array
#sorting for difference surf and intersect surf
for i, line in enumerate(name_array):
	name = f'../VASP/scores/{raw_name[i]}'
	a=name_array[i]
	with open(name, 'r') as ifile:
		for line in ifile:
			pass
		iline = line.split()[3]
		inter.append([a[0],a[1],a[2],a[3],a[4],a[5],a[6],iline])
		pockets.append([a[0],a[1]])

#######NOTE: pockets positioning will NOT match inter array after sorting; should eliminate
#######why do i use two arrays with identical information???? lets make it one array
#######and then do all iterative operations at once

#determine the highest intersect/lowest diff score
print('-----Sorting scores-----')
inter_sort=[]
inter_final=[]
inter_sort = sorted(inter,key=lambda x: x[-1],reverse=True)

#get the volume of the target pocket, and then for each pocket
with open('../initial_structures/vol.txt', 'r') as f:
	for line in f:
		pass
	vline = line.split()[-1]

vols={}
for p in inter_sort:#pockets:
	with open(f'../VASP/pockets/{p[0]}_{p[1]}.vol.txt', 'r') as f:
		for line in f:
			pass
		key = f'{p[0]} {p[1]}'
		vol = line.split()[3]
		vols.update({key : vol})


#take the first entry of each name since the list is ordered
#from highest to lowest. also generate the hits dictionary to
#later pull out the number of hits satisfying the filter
tracker=[]
hits={}

for j in range(len(inter_sort)):
	count=0
	if [inter_sort[j][0],inter_sort[j][1]] not in tracker:
		tracker.append([inter_sort[j][0],inter_sort[j][1]])
		inter_final.append(inter_sort[j])
	
	if float(inter_sort[j][4])/float(vline) > float(filt):
		count += 1
	
	key = f'{inter_sort[j][0]} {inter_sort[j][1]}'
	hits.update({key : count})

#	for line in inter_sort:
#		if pockets[j][0] == line[0]:
#			if pockets[j][1] == line[1]:
#				if [line[0], line[1]] not in tracker:
#					tracker.append([line[0], line[1]])
#					inter_final.append(line)
#					if float(line[4])/float(vline) > float(filt):
#						count+=1
#				else:
#					if float(line[4])/float(vline) > float(filt):
#						count+=1
#	key = f'{pockets[j][0]} {pockets[j][1]}'
#	hits.update({key : count})


###what does this loop do?
#for line in inter_sort:
#	if pockets[2][0] == line[0]:
#		if pockets[2][1] == line[1]:
#			continue


#put vline in score file and # of hits that satisfy filter
#combine into output file
print('-----Outputting scores-----')
with open('../score.txt','w') as outfile:
	outfile.write('PDB   Pock   Target Vol  Pock Vol     Int Vol   Int %  # hits\n')
	for i in range(len(inter_final)):
		a=inter_final[i][0]
		p=inter_final[i][1]
		b=hits.get(f'{a} {p}')
		v=vols.get(f'{a} {p}')
		c=float(inter_final[i][4])/float(vline)
		d=inter_final[i][4]
		out=f"{a:<6}{p:<7}{vline:{12}.{8}}{v:{12}.{8}} {d:{9}.{7}} {c:{6}.{3}}    {b}\n"
		outfile.write(out)

#put new scorefile into master scorefile if it exists, otherwise create it
s = np.genfromtxt('../score.txt', dtype='unicode', skip_header=1)
if os.path.exists('../winners/master_scorefile.txt'):
	mode='a'
else:
	mode='w'

with open('../winners/master_scorefile.txt', mode) as f:
	if mode == 'w':
		f.write('PDB   Pock   Target Vol  Pock Vol    Int Vol   Int %  #hits\n')
	for i in range(len(s)):
		pdb=f'{s[i][0]:<6}'
		pock=f'{s[i][1]:<7}'
		v1=f'{s[i][2]:{12}.{8}}'
		v2=f'{s[i][3]:{12}.{8}}'
		d=f'{s[i][4]:{11}.{7}}'
		dp=f'{s[i][5]:{7}.{3}}'
		h=f'{s[i][6]:<6}'
		f.write(pdb+pock+v1+v2+d+dp+h+'\n')

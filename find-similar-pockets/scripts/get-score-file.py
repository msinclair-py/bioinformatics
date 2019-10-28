#!/usr/bin/env python
import glob

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
	raw_name.append(name.split('/')[3])
	item = name.split('/')[3][:4]
	if item not in pdbs:
		pdbs.append(item)
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
		inter.append([a[0],a[2],a[3],a[4],a[1],a[5],iline])
		pockets.append([a[0],a[1]])

#determine the highest intersect/lowest diff score
print('-----Sorting scores-----')
inter_sort=[]
inter_final=[]
inter_sort = sorted(inter,key=lambda x: x[-1],reverse=True)

#get the volume of the target pocket, and then for each pocket
with open('../initial_structures/vol.txt', 'r') as f:
	for line in f:
		pass
	vline = line.split()[3]

vols={}
for p in pockets:
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

for j in range(len(pockets)):
	count=0
	for line in inter_sort:
		if pockets[j][0] == line[0]:
			if pockets[j][1] == line[4]:
				if [line[0], line[4]] not in tracker:
					tracker.append([line[0], line[4]])
					inter_final.append(line)
					if float(line[6])/float(vline) > float(filt):
						count+=1
				else:
					if float(line[6])/float(vline) > float(filt):
						count+=1
	key = f'{pockets[j][0]} {pockets[j][1]}'
	hits.update({key : count})

for line in inter_sort:
	if pockets[2][0] == line[0]:
		if pockets[2][1] == line[4]:
			continue

#put vline in score file and # of hits that satisfy filter
#combine into output file
print('-----Outputting scores-----')
with open('../score.txt','w') as outfile:
	outfile.write('PDB   Pock   Target Vol  Pock Vol     Int Vol   Int %  # hits\n')
	for i in range(len(inter_final)):
		a=inter_final[i][0]
		p=inter_final[i][4]
		b=hits.get(a+' '+p)
		v=vols.get(a+' '+p)
		c=float(inter_final[i][6])/float(vline)
		d=inter_final[i][6]
		out=f"{a:<6}{p:<7}{vline:{12}.{8}}{v:{12}.{8}} {d:{9}.{7}} {c:{6}.{3}}    {b}\n"
		outfile.write(out)


#!/usr/bin/env python
import glob, shutil, os
import numpy as np

#make sure all the winners directories exist
if not os.path.exists('../winners/winners1'):
	os.mkdir('../winners/winners1')
if not os.path.exists('../winners/scorefiles'):
	os.mkdir('../winners/scorefiles')

#get the name/pocket for each incoming structure
uniques=[]
names=[]
for name in glob.glob('../PDB/*.aligned'):
	init=name.split('/')[2]
	if init.split('_')[1] not in uniques:
		uniques.append(init.split('_')[1])
	names.append([init.split('_')[1],init.split('_')[2][:5]])

#find out what the score text file will be (score-#.txt)
scorefiles=[]
for scores in glob.glob('../winners/scorefiles/*'):
	scorefiles.append(scores.split('-')[1][:-4])

scorefiles = sorted([int(x) for x in scorefiles])
if os.path.exists('../winners/scorefiles/score-1.txt'):
	new_score = scorefiles[-1] + 1
else:
	new_score = 1

#find out which folders the pdbs are going to (winners#)
winners=[]
for winner in glob.glob('../winners/winners*'):
	winners.append(winner.split('/')[2][7:])

winners = sorted([int(x) for x in winners])
latest_winners = winners[-1]
new_winners = winners[-1] + 1


#check number of files in most recent winners directory
count=0
for pdb in glob.glob(f'../winners/winners{latest_winners}/*.ent'):
	count+=1
num=100-int(count)
#num is the maximum number of winning pdbs we can put in the next directory

#figure out where everything is going, store locations in an array for master
#scorefile

locations=[]
if num > len(uniques):
	for name in uniques:
		locations.append([name, latest_winners])
		for stuff in glob.glob(f'../PDB/*{name}*'):
			shutil.move(stuff, f'../winners/winners{latest_winners}')
else:
	os.mkdir(f'../winners/winners{new_winners}')
	for name in uniques[num:]:
		locations.append([name, latest_winners])
		for stuff in glob.glob(f'../PDB/*{name}*'):
			shutil.move(stuff, f'../winners/winners{latest_winners}')
	for name in uniques[:num]:
		locations.append([name, new_winners])
		for stuff in glob.glob(f'../PDB/*{name}*'):
			shutil.move(stuff, f'../winners/winners{new_winners}')

#append new scores and locations to master scorefile
scores = np.genfromtxt('../score.txt', dtype='unicode', skip_header=1)

#if master scorefile exists, open as append, else open as write
#write will overwrite existing files while append adds to existing file
if os.path.exists('../winners/master_scorefile.txt'):
	mode='a'
else:
	mode='w'

#if we just made the master scorefile then add header and then append each score
with open('../winners/master_scorefile.txt', mode) as f:
	if mode == 'w':
		f.write('PDB   Pock   Target Vol  Pock Vol     Int Vol   Int %  #hits  winners directory\n')
	for s in scores:
		for i in range(len(locations)):
			if locations[i][0] in s:
				l=locations[i][1]
				out=f'{s[0]:<6}{s[1]:<7}{s[2]:{12}.{8}}{s[3]:{12}.{8}}{s[4]:{11}.{7}}{s[5]:{7}.{3}}{s[6]:<7}{l}\n'
				f.write(out)

#now that master scorefile is updated move the score.txt file to the scorefiles archive
shutil.move('../score.txt', f'../winners/scorefiles/score-{new_score}.txt')

#finally we need to clear out the VASP directory or all the structures in it will get run again with
#the next analysis
shutil.rmtree('../VASP/scores')
shutil.rmtree('../VASP/pockets')
for surf in glob.glob('../VASP/*'):
	os.remove(surf)

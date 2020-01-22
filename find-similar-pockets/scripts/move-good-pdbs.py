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
for name in glob.glob('../PDB/aligned.*'):
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
for pdb in glob.glob(f'../winners/winners{latest_winners}/aligned.*'):
	count+=1
num=100-int(count)
#num is the maximum number of winning pdbs we can put in the next directory

if num > len(uniques):
	for name in uniques:
		print(name)
		for stuff in glob.glob(f'../PDB/*{name}*'):
			shutil.move(stuff, f'../winners/winners{latest_winners}')
else:
	os.mkdir(f'../winners/winners{new_winners}')
	for name in uniques[num:]:
		for stuff in glob.glob(f'../PDB/*{name}*'):
			shutil.move(stuff, f'../winners/winners{latest_winners}')
	for name in uniques[:num]:
		for stuff in glob.glob(f'../PDB/*{name}*'):
			shutil.move(stuff, f'../winners/winners{new_winners}')

#now that master scorefile is updated move the score.txt file to the scorefiles archive
shutil.move('../score.txt', f'../winners/scorefiles/score-{new_score}.txt')


#finally we need to clear out the VASP directory or all the structures in it will get run again with
#the next analysis
shutil.rmtree('../VASP/scores')
shutil.rmtree('../VASP/pockets')
for sf in glob.glob('../VASP/*'):
	os.remove(sf)

#!/usr/bin/env python
from itertools import islice, zip_longest

#n is the number of entries per list and is the modifiable of this
#script
n = 100

#this function determines the number of lines in the input file
def file_len(fname):
	for i, l in enumerate(fname):
		pass
	return i + 1

#this function aggregates elements from the iterable.
#documentation on the zip_longest function from itertools
#explains this better than I am capable of
def grouper(iterable, n, fillvalue=""):
	args = [iter(iterable)] * n
	return zip_longest(fillvalue=fillvalue, *args)

#the blank array, l, is instantiated so that we can fill it with 
#input lines from the main pdb list. strip is used to ensure no
#whitespace is accidentally included which can mess up the final step
l = []
with open('../sublists/pdb_list.txt', 'r') as infile:
	for line in infile:
		l.append(line.strip())

#finally we iterate over the grouper function results to a series
#of output sublists
it = iter(l)
for num, chunk in enumerate(iter(lambda: list(islice(it, n)), []), 1):
	with open('../sublists/sublist{}.txt'.format(num), 'w') as fout:
		fout.write("\n".join(chunk))

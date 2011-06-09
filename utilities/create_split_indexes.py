from numpy import *
import sys
from random import shuffle
from os import system


def read_input(genefn):
    classesinds = [[],[]]

    for ind,line in enumerate(open(genefn,'r')):
	if line[0] == '#':
	    continue

	parts = line.rstrip().rsplit()
	if parts[0] == '1':
	    classesinds[0].append(ind)
	else:
	    classesinds[1].append(ind)


    return classesinds




def create_splits(classesinds,trainprop,numrepeats):

    cuts = []
    Npos = len(classesinds[0])
    Nneg = len(classesinds[1])
    posinds = range(Npos)
    neginds = range(Nneg)


    for repeatind in range(numrepeats):

	shuffle(posinds)
	shuffle(neginds)

	print 'Splitting data with repeat #%d' % repeatind

	onecut = []

	for ind in posinds[:int(trainprop*Npos)]:
	    onecut.append(classesinds[0][ind])

	for ind in neginds[:int(trainprop*Nneg)]:
	    onecut.append(classesinds[1][ind])

	cuts.append(onecut)

    return cuts



def create_cv_splits(classesinds,numrepeats):

    cuts = []
    Npos = len(classesinds[0])
    Nneg = len(classesinds[1])
    posinds = range(Npos)
    neginds = range(Nneg)

    shuffle(posinds)
    shuffle(neginds)

    numpos_percut = Npos/numrepeats
    numneg_percut = Nneg/numrepeats

    for repeatind in range(numrepeats):


	print 'Splitting data with repeat #%d' % repeatind

	onecut = []

	for ind in range(Npos):
	    if ind not in posinds[repeatind*numpos_percut:((repeatind+1)*numpos_percut if repeatind < (numrepeats-1) else -1)]:
		onecut.append(classesinds[0][ind])

	for ind in range(Nneg):
	    if ind not in neginds[repeatind*numneg_percut:((repeatind+1)*numneg_percut if repeatind < (numrepeats-1) else -1)]:
		onecut.append(classesinds[1][ind])

	cuts.append(onecut)

    return cuts



def store_cuts_indexes(cuts,splitfns):

    cutfile = open(splitfns,'w')

    for cutind,cut in enumerate(cuts):
	for ind in cut:
	    print >> cutfile, '%d ' % ind,
	print >> cutfile


    cutfile.close()






if __name__=='__main__':
    args = sys.argv[1:]


    if len(args) == 0:
	print '=====INPUTS=====\ninputfn = args[0]\nsplitfns = args[1]\ntrainproportion = float(args[2])\nnumrepeats = int(args[3])\ndocv = args[4] == "1"\ndo_part_and_cv = args[5] == "1"\n================='
	exit()


    inputfn = args[0]
    splitfns = args[1]
    trainproportion = float(args[2])
    numrepeats = int(args[3])
    docv = args[4] == '1'
    do_part_and_cv = args[5] == '1'

    print 'Reading data'
    classinds = read_input(inputfn)
    print 'Splitting data'
    if docv:
	cuts = create_cv_splits(classinds,numrepeats)
    elif do_part_and_cv:
	cuts = create_cv_and_part_splits(classinds,numrepeats)
    else:
	cuts = create_splits(classinds,trainproportion,numrepeats)

    print 'Storing Cuts'
    store_cuts_indexes(cuts,splitfns)







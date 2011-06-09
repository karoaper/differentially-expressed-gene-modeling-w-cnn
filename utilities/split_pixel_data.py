from numpy import *
import sys
from random import shuffle
from os import system

def read_input(inputfn):
    input = open(inputfn,'r').readlines()

    header = ''
    data = [[],[]]

    for line in input:
	if line[0] == '#':
	    header += line
	else:
	    parts = line.rstrip().rsplit()
	    if parts[0] == '1':
		data[0].append(line)
	    else:
		data[1].append(line)


    shuffle(data[0])
    shuffle(data[1])
    return [header,data]


def split_data(data,trainprop,numrepeats):

    cuts = []
    Npos = len(data[0])
    Nneg = len(data[1])
    for repeatind in range(numrepeats):
	print 'Splitting data with repeat #%d' % repeatind

	shuffle(data[0])
	shuffle(data[1])

	onecut = [[],[]]

	onecut[0].extend(data[0][:int(trainprop*Npos)])
	onecut[0].extend(data[1][:int(trainprop*Nneg)])

	onecut[1].extend(data[0][int(trainprop*Npos):])
	onecut[1].extend(data[1][int(trainprop*Nneg):])

	cuts.append(onecut)

    return cuts


def store_cuts_update_confs(conffn,inputfn,cuts,header):

    conffn_parts = conffn.rsplit('.')
    inputfn_parts = inputfn.rsplit('.')
    for cutind,cut in enumerate(cuts):
#	system('cp %s %s_cut%d.conf' % (conffn,''.join(conffn_parts[:-1]),cutind))
	Ntrain = len(cut[0])
	Ntest = len(cut[1])


	trainfn = '%s_cut%d_train.dat' % (''.join(inputfn_parts[:-1]),cutind)
	testfn = '%s_cut%d_test.dat' % (''.join(inputfn_parts[:-1]),cutind)
	outputfile_train = open(trainfn,'w')
	outputfile_test = open(testfn,'w')
        trainfn = trainfn.rsplit('/')[-1]
        testfn = testfn.rsplit('/')[-1]

#	print conffn
	
#	print 'sed \'s/testing_size[ ]* =[ ]*[0-9]*/testing_size = %d/\' %s | sed \'s/training_size[ ]* =[ ]*[0-9]*/training_size = %d/\' | sed \'s/train_name[ ]* =[ ]*.*/train_name = %s/\' | sed \'s/test_name[ ]* =[ ]*.*/test_name = %s/\' >  %s_cut%d.conf' % (Ntest,conffn,Ntrain,trainfn,testfn,''.join(conffn_parts[:-1]),cutind)
	system('sed \'s/testing_size[ ]* =[ ]*[0-9]*/testing_size = %d/\' %s | sed \'s/training_size[ ]* =[ ]*[0-9]*/training_size = %d/\' | sed \'s/train_name[ ]* =[ ]*.*/train_name = %s/\' | sed \'s/test_name[ ]* =[ ]*.*/test_name = %s/\' >  %s_cut%d.conf' % (Ntest,conffn,Ntrain,trainfn,testfn,''.join(conffn_parts[:-1]),cutind))

	print >> outputfile_train, header,
	print >> outputfile_test, header,

	for datum in cut[0]:
	    print >> outputfile_train, datum,
	for datum in cut[1]:
	    print >> outputfile_test, datum,

	outputfile_train.close()
	outputfile_test.close()






if __name__=='__main__':
    args = sys.argv[1:]

    inputfn = args[0]
    conffn = args[1]
    trainproportion = float(args[2])
    numrepeats = int(args[3])

    print 'Reading data'
    [header,data] = read_input(inputfn)
    print 'Splitting data'
    cuts = split_data(data,trainproportion,numrepeats)
    print 'Storing Cuts'
    store_cuts_update_confs(conffn,inputfn,cuts,header)


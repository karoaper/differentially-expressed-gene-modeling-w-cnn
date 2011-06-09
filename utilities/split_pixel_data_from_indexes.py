from numpy import *
import sys
from random import shuffle
from os import system

def read_input(inputfn):
    input = open(inputfn,'r').readlines()

    header = ''
    data = []

    for line in input:
	if line[0] == '#':
	    header += line
	else:
	    data.append(line)

    return [header,data]


def split_data(data,cutsinds):

    cuts = []
    for repeatind,cutinds in enumerate(cutsinds):
	print 'Splitting data with repeat #%d' % repeatind

	onecut = [[],[]]


	onecut[0].extend([data[i] for i in range(len(data)) if i in cutinds])
	onecut[1].extend([data[i] for i in range(len(data)) if not i in cutinds])

	cuts.append(onecut)

    return cuts


def store_cuts_update_confs(conffn,inputfn,cuts,header,outputfolder):

    inputfn_parts = inputfn.rsplit('/')[-1].rsplit('.')

#   inputfn_nofolder = inputfn.rsplit('/')[-1].


#   inputfn_folder = '\/'.join(inputfn.rsplit('/')[:-1])
    conffn = conffn.replace(' ','\ ')
    conffn_parts = conffn.rsplit('.')
    for cutind,cut in enumerate(cuts):
#	system('cp %s %s_cut%d.conf' % (conffn,''.join(conffn_parts[:-1]),cutind))
	Ntrain = len(cut[0])
	Ntest = len(cut[1])


#	trainfn = '%s_cut%d_train.dat' % (''.join(inputfn_parts[:-1]),cutind)
#	testfn = '%s_cut%d_test.dat' % (''.join(inputfn_parts[:-1]),cutind)

	trainfn = '%s/%s_cut%d_train.dat' % (outputfolder,''.join(inputfn_parts[:-1]),cutind)
	testfn = '%s/%s_cut%d_test.dat' % (outputfolder,''.join(inputfn_parts[:-1]),cutind)


	outputfile_train = open(trainfn,'w')
	outputfile_test = open(testfn,'w')
	trainfn = trainfn.rsplit('/')[-1]
	testfn = testfn.rsplit('/')[-1]

#	print conffn

#	print 'sed \'s/testing_size[ ]* =[ ]*[0-9]*/testing_size = %d/\' %s | sed \'s/training_size[ ]* =[ ]*[0-9]*/training_size = %d/\' | sed \'s/train_name[ ]* =[ ]*.*/train_name = %s/\' | sed \'s/test_name[ ]* =[ ]*.*/test_name = %s/\' | sed \'s/root[ ]*=[ ]*.*/root = %s/\' >  %s_cut%d.conf' % (Ntest,conffn,Ntrain,trainfn,testfn,inputfn_folder,''.join(conffn_parts[:-1]),cutind)
#	print('sed \'s/testing_size[ ]* =[ ]*[0-9]*/testing_size = %d/\' %s | sed \'s/training_size[ ]* =[ ]*[0-9]*/training_size = %d/\' | sed \'s/train_name[ ]* =[ ]*.*/train_name = %s/\' | sed \'s/test_name[ ]* =[ ]*.*/test_name = %s/\' | sed \'s/root[ ]*=[ ]*.*/root = %s/\' >  %s_cut%d.conf' % (Ntest,conffn,Ntrain,trainfn,testfn,outputfolder,''.join(conffn_parts[:-1]),cutind))

	system('sed \'s/testing_size[ ]* =[ ]*[0-9]*/testing_size = %d/\' %s | sed \'s/training_size[ ]* =[ ]*[0-9]*/training_size = %d/\' | sed \'s/train_name[ ]* =[ ]*.*/train_name = %s/\' | sed \'s/test_name[ ]* =[ ]*.*/test_name = %s/\' | sed \'s/root[ ]*=[ ]*.*/root = %s/\' >  %s_cut%d.conf' % (Ntest,conffn,Ntrain,trainfn,testfn,outputfolder.replace('/','\/'),''.join(conffn_parts[:-1]),cutind))




	print >> outputfile_train, header,
	print >> outputfile_test, header,

	for datum in cut[0]:
	    print >> outputfile_train, datum,
	for datum in cut[1]:
	    print >> outputfile_test, datum,

	outputfile_train.close()
	outputfile_test.close()



def read_cuts(cutindsfn):
    cuts = [[int(x) for x in line.rstrip().rsplit()]  for line in open(cutindsfn,'r')]
    return cuts


if __name__=='__main__':

    args = sys.argv[1:]

    if len(args) == 0:
	print '=====INPUTS=====\ninputfn = args[0]\nconffn = args[1]\ncutindsfn = args[2]\noutputfolder = args[3]\n================='
	exit()


    inputfn = args[0]
    conffn = args[1]
    cutindsfn = args[2]
    outputfolder = args[3]

    print 'Reading data'
    [header,data] = read_input(inputfn)
    print 'Reading cut indexes'
    cuts = read_cuts(cutindsfn)
    print 'Splitting data'
    cuts = split_data(data,cuts)
    print 'Storing Cuts'
    store_cuts_update_confs(conffn,inputfn,cuts,header,outputfolder)




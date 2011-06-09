from numpy import *
from re import *
import sys
from random import shuffle


def read_genes(genefn,winsize):
    genes = [[x.rstrip().rsplit(' ')[0],x.rstrip().rsplit(' ')[1],int(x.rstrip().rsplit(' ')[2]),x.rstrip().rsplit(' ')[3],int(x.rstrip().rsplit(' ')[4])] for x in open(genefn,'r').readlines()]


    chroms = unique([x[1] for x in genes])
    chroms_up = unique([x[1] for x in genes if x[-1] == 1])
    chroms_down = unique([x[1] for x in genes if x[-1] == -1])


    genedict = {}

    for chr in chroms:
	genedict[chr] = []


    for gene in genes:
	genedict[gene[1]].append(gene)
	genedict[gene[1]][-1].extend([(gene[2]-winsize/2) if gene[3] == '+' else (gene[2] + winsize/2),(gene[2] + winsize/2) if gene[3] == '+' else (gene[2] - winsize/2)])




    return genedict,chroms_up,chroms_down



def read_peaks(buffersize,peaksfns):

    allpeaks = {}

    for peakfn in peaksfns:
	allpeaks[peakfn] = {}
	peaksfile = open(peakfn,'r')
	for peak in peaksfile:
	    if peak[0:5] == 'track':
		continue

	    parts = peak.rstrip().rsplit('\t')

	    if allpeaks[peakfn].has_key(parts[0]):
		allpeaks[peakfn][parts[0]].append([int(parts[1])-buffersize,int(parts[2])+buffersize,float(parts[4])])
	    else:
		allpeaks[peakfn][parts[0]] = [[int(parts[1])-buffersize,int(parts[2])+buffersize,float(parts[4])]]


	for key in allpeaks[peakfn].keys():
	    allpeaks[peakfn][key].sort()

    return allpeaks



def remove_overlaps(genes):
    numpre = 0
    numpost = 0
#   chroms = gunique([x[1] for x in genes])


    for chr in genes:
	numpre += len(genes[chr])
	ind = 0
	while ind < len(genes[chr]):
	    gene = genes[chr][ind]

	    isoverlapped = False
	    orient = gene[3]
	    for othergene in [x for x in genes[chr] if x[3] == orient]:
		if othergene[0] == gene[0]:
		    continue

		if (orient == '+' and (gene[6] > othergene[5] and gene[6] <= othergene[6] or gene[5] >= othergene[5] and gene[5] < othergene[6])) or (orient == '-' and (gene[6] < othergene[5] and gene[6] >= othergene[6] or gene[5] <= othergene[5] and gene[5] > othergene[6])):
		    isoverlapped = True
		    break

	    if isoverlapped:
		genes[chr].pop(ind)
	    else:
		ind += 1

	numpost += len(genes[chr])
    print 'Number of genes prior to/after culling overlapping genes: %d/%d' % (numpre,numpost)
    return genes






def remove_redundants(genes):
    numpre = 0
    numpost = 0

    numpos = 0
    numneg = 0
    for genekey in genes.keys():
	for gene in genes[genekey]:
	    if gene[4] == 1:
		numpos += 1
	    else:
		numneg += 1

    print 'Number of genes prior to culling redundant genes: %d(%d/%d)' % (numpos+numneg,numpos,numneg)


    for chr in genes:
	numpre += len(genes[chr])
	ind = 0
	while ind < len(genes[chr]):
	    gene = genes[chr][ind]

	    isredundant = False
	    orient = gene[3]
	    for othergene in [x for x in genes[chr] if x[3] == orient]:
		if othergene[0] == gene[0]:
		    continue

		if gene[2] == othergene[2]:
		    isredundant = True
		    break

	    if isredundant:
		genes[chr].pop(ind)
	    else:
		ind += 1

	numpost += len(genes[chr])


    numpos = 0
    numneg = 0
    for genekey in genes.keys():
	for gene in genes[genekey]:
	    if gene[4] == 1:
		numpos += 1
	    else:
		numneg += 1


    print 'Number of genes after culling redundant genes: %d(%d/%d)\n' % (numpost,numpos,numneg)

    return genes





def filter_by_peaks(genes,allpeaks):
    numpre = 0
    numpost = 0
#   chroms = unique([x[1] for x in genes])

    numpos = 0
    numneg = 0
    for genekey in genes.keys():
	for gene in genes[genekey]:
	    if gene[4] == 1:
		numpos += 1
	    else:
		numneg += 1

    print 'Number of genes prior to culling blank genes: %d(%d/%d)' % (numpos+numneg,numpos,numneg)


    for chr in genes:
	numpre += len(genes[chr])
	ind = 0
	while ind < len(genes[chr]):
	    gene = genes[chr][ind]

	    is_empty = True
	    orient = gene[3]
	    for peakfn in allpeaks:
		for peak in allpeaks[peakfn][chr]:
		    if (orient == '+' and (peak[0] >= gene[5] and peak[0] < gene[6] or peak[1] > gene[5] and peak[1] <= gene[6])) or (orient == '-' and (peak[1] <= gene[5] and peak[1] > gene[6] or peak[0] >= gene[6] and peak[0] < gene[5])):
			is_empty = False
			break

	    if is_empty:
		genes[chr].pop(ind)
	    else:
		ind += 1

	numpost += len(genes[chr])

    numpos = 0
    numneg = 0
    for genekey in genes.keys():
	for gene in genes[genekey]:
	    if gene[4] == 1:
		numpos += 1
	    else:
		numneg += 1


    print 'Number of genes after culling blank genes: %d(%d/%d)\n' % (numpost,numpos,numneg)
    return genes


def separate_by_upordown(genes,chroms_up,chroms_down):


    print genes.keys()


    print chroms_up,chroms_down

    genedict = {1:{},-1:{}}

    chrsize_up = {}
    chrsize_down = {}

    sizeup = 0
    sizedown = 0





    for chr in chroms_up:
	genedict[1][chr] = []
    for chr in chroms_down:
	genedict[-1][chr] = []


    for genekey in genes.keys():
	for gene in genes[genekey]:
	    genedict[gene[4]][genekey].append(gene)




    return genedict



def sample_genes(genes,N,classratio_up,chroms_up,chroms_down):


    numpos = 0
    numneg = 0
    for genekey in genes.keys():
	for gene in genes[genekey]:
	    if gene[4] == 1:
		numpos += 1
	    else:
		numneg += 1

    print "Number of positive genes: %d\nNumber of negative genes: %d" % (numpos,numneg)

    if classratio_up == 'fromset':
	classratio_up = numpos*1.0/numpost
	print "Computed classratio_up is %f" % (classratio_up)
    else:
	classratio_up = float(classratio_up)


#   chroms_up = unique([x[1] for x in genes if x[-1] == 1])
#   chroms_down = unique([x[1] for x in genes if x[-1] == -1])


    for chr in chroms_up:
	sizeup += len(genes[1][chr])
	chrsize_up[chr] = len(genes[1][chr])
    for chr in chroms_down:
	sizedown += len(genes[-1][chr])
	chrsize_down[chr] = len(genes[-1][chr])

    sampledgenes = {1:{},-1:{}}

    if outputrest:
	sampledgenes_rest = {1:{},-1:{}}


    Nup = N*classratio_up
    Ndown = N-Nup

    Nup_real= 0
    Ndown_real = 0

    for chr in chroms_up:
	numperchr = Nup * chrsize_up[chr]/sizeup
	allNs = range(chrsize_up[chr])
	shuffle(allNs)

	Nup_real += len(allNs[:int(ceil(numperchr))])
	sampledgenes[1][chr] = [genes[1][chr][i] for i in allNs[:int(ceil(numperchr))]]

	if outputrest:
	    sampledgenes_rest[1][chr] = [genes[1][chr][i] for i in allNs[int(ceil(numperchr)):]]


    for chr in chroms_down:
	numperchr = Ndown * chrsize_down[chr]/sizedown
	allNs = range(chrsize_down[chr])
	shuffle(allNs)

	Ndown_real += len(allNs[:int(ceil(numperchr))])
	sampledgenes[-1][chr] = [genes[-1][chr][i] for i in allNs[:int(ceil(numperchr))]]

	if outputrest:
	    sampledgenes_rest[-1][chr] = [genes[-1][chr][i] for i in allNs[int(ceil(numperchr)):]]


    return sampledgenes




def store_genes(genes,outfn,outputrest,chroms_up,chroms_down):
    outputfile = open(outfn,'w')

#   chroms_up = unique([x[1] for x in genes if x[-1] == 1])
#   chroms_down = unique([x[1] for x in genes if x[-1] == -1])



    Nup = 0
    Ndown  = 0
    for chr in chroms_up:
	Nup += len(genes[1][chr])
	for gene in genes[1][chr]:
	    print >> outputfile, '%s %s %d %s %d' % (gene[0],gene[1],gene[2],gene[3],gene[4])

    for chr in chroms_down:
	Ndown += len(genes[-1][chr])
	for gene in genes[-1][chr]:
	    print >> outputfile, '%s %s %d %s %d' % (gene[0],gene[1],gene[2],gene[3],gene[4])


    print '%d Up Genes/%d Down Genes Stored' % (Nup,Ndown)


    if outputrest:
	outputfile = open('%s_rest.dat' % outfn,'w')

	Nup = 0
	Ndown  = 0
	for chr in chroms_up:
	    Nup += len(sampledgenes_rest[1][chr])
	    for gene in sampledgenes_rest[1][chr]:
		print >> outputfile, '%s %s %d %s %d' % (gene[0],gene[1],gene[2],gene[3],gene[4])

	for chr in chroms_down:
	    Ndown += len(sampledgenes_rest[-1][chr])
	    for gene in sampledgenes_rest[-1][chr]:
		print >> outputfile, '%s %s %d %s %d' % (gene[0],gene[1],gene[2],gene[3],gene[4])


	print 'The remaining %d Up Genes/%d Down Genes Stored' % (Nup,Ndown)



def main(args):


    if len(args) == 0:
	print '=====INPUTS=====\ngenefn = args[0]\noutfn = args[1]\nclassratio_up = args[2]\nN = args[3]\nwinsize = int(args[4])\nnooverlaps = args[5] == "1"\noutputrest = args[6] == "1"\nfilterbypeaks = args[7] == "1"\n================='
	exit()



    genefn = args[0]
    outfn = args[1]
    classratio_up = args[2]
    N = args[3]
    winsize = int(args[4])
    nooverlaps = args[5] == '1'
    noredundants = args[6] == '1'
    outputrest = args[7] == '1'
    filterbypeaks = args[8] == '1'







    genes,chroms_up,chroms_down = read_genes(genefn,winsize)




    if filterbypeaks:
	buffersize = int(args[9])
	numpeakfns = int(args[10])
	peaksfns = args[11:11+numpeakfns]
	allpeaks = read_peaks(buffersize,peaksfns)





    if nooverlaps:
	genes = remove_overlaps(genes)

    elif noredundants:
	genes = remove_redundants(genes)


    if filterbypeaks:
	genes = filter_by_peaks(genes,allpeaks)







    numpos = 0
    numneg = 0
    for genekey in genes.keys():
	for gene in genes[genekey]:
	    if gene[4] == 1:
		numpos += 1
	    else:
		numneg += 1




    genes = separate_by_upordown(genes,chroms_up,chroms_down)

    print "Number of positive genes: %d\nNumber of negative genes: %d" % (numpos,numneg)



    if N != 'all':
	N = int(N)
	genes = sample_genes(genes,N,chroms_up,chroms_down)



    store_genes(genes,outfn,outputrest,chroms_up,chroms_down)





if __name__ == '__main__':



    main(sys.argv[1:])










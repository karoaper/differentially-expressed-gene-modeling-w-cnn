from numpy import *
import os
import sys
import matplotlib
matplotlib.use("Agg")
from pylab import *

from matplotlib.backends.backend_pdf import PdfPages



def ma(data,masamples):
    w=ones(masamples,'d')

    y=convolve(w/w.sum(),data,mode='same')

    newl = size(data)-masamples+1
    newl2 = ceil((len(y)-newl)/2.0)
    return y[newl2:newl2+newl]




def main(args):

    if len(args) == 0:
	print '=====INPUTS=====\noutputfn = args[0]\nmasamples = int(args[1])\nfiles = args[2:]\n================='
	exit()

    outputfn = args[0]
    masamples = int(args[1])
    topresultsnum = int(args[2])
    files = args[3:]

    all_maxtest = []
    all_maxtest_smooth = []

    pdf = PdfPages(outputfn)

    average_maxtrain_accs_inbins = []
    average_maxtrain_accs_inbins_smooth = []
    alltopindeces = []
    alltopindeces_smooth = []

    for fnind,fn in enumerate(files):
	trainvals = []
	testvals = []
	trainlines = os.popen('egrep " correct=" %s' % fn.replace(' ','\ ')).readlines()
	for line in trainlines:
	    trainvals.append(float(line.rstrip('% \n').rsplit('=')[-1]))
	trainlines = os.popen('egrep "test_correct=" %s' % fn.replace(' ','\ ')).readlines()
	for line in trainlines:
	    testvals.append(float(line.rstrip('% \n').rsplit('=')[-1]))

	trainvals = array(trainvals)
	testvals = array(testvals)

#	print trainvals

	trainvals_smooth = ma(trainvals,masamples)
	testvals_smooth = ma(testvals,masamples)


	n = len(trainvals)
	n1 = len(trainvals_smooth)


	maxtrain = [argmax(trainvals),trainvals[argmax(trainvals)]]
	maxtest = [argmax(testvals),testvals[argmax(testvals)]]

	all_maxtest.append(testvals[argmax(testvals)])

	maxtrain_smooth = [argmax(trainvals_smooth),trainvals_smooth[argmax(trainvals_smooth)]]
	maxtest_smooth = [argmax(testvals_smooth),testvals_smooth[argmax(testvals_smooth)]]

	all_maxtest_smooth.append(testvals_smooth[argmax(testvals_smooth)])
	print '%g,%g' % (all_maxtest[-1],all_maxtest_smooth[-1])
	rangeofindeces = sorted(argsort(testvals)[-topresultsnum:])
	rangeofindeces_smooth = sorted(argsort(testvals_smooth)[-topresultsnum:])


	alltopindeces.append(rangeofindeces)
	alltopindeces_smooth.append(rangeofindeces_smooth)

	f1 = figure(2)
	title('Optimal trainvals')
	for val in alltopindeces[-1]:
#	    print val,len(trainvals),len(testvals)
	    plot(fnind,trainvals[val],'ko')


	f1 = figure(3)
	title('Optimal smoothed trainvals')

	for val in alltopindeces_smooth[-1]:
#	    print val,len(trainvals_smooth),len(testvals_smooth)
	    plot(fnind,trainvals_smooth[val],'ko')


	f1 = figure(4)
	title('Optimal testvals')

	for val in alltopindeces[-1]:
#	    print val,len(trainvals),len(testvals)
	    plot(fnind,testvals[val],'ko')


	f1 = figure(5)
	title('Optimal smoothed testvals')


	for val in alltopindeces_smooth[-1]:
#	    print val,len(trainvals_smooth),len(testvals_smooth)
	    plot(fnind,testvals_smooth[val],'ko')






#	print rangeofindeces

	rangeofindeces_partitions = [[rangeofindeces[0],rangeofindeces[0]]] + [[rangeofindeces[i],rangeofindeces[i+1]] for i in range(len(rangeofindeces)-1) if rangeofindeces[i+1]-rangeofindeces[i]>5] + [[rangeofindeces[-1],rangeofindeces[-1]]]
	rangeofindeces_partitions_smooth = [[rangeofindeces_smooth[0],rangeofindeces_smooth[0]]] + [[rangeofindeces_smooth[i],rangeofindeces_smooth[i+1]] for i in range(len(rangeofindeces_smooth)-1) if rangeofindeces_smooth[i+1]-rangeofindeces_smooth[i]>5] + [[rangeofindeces_smooth[-1],rangeofindeces_smooth[-1]]]

#	print rangeofindeces_partitions
	bins = [[rangeofindeces_partitions[i][1],rangeofindeces_partitions[i+1][0]] for i in range(len(rangeofindeces_partitions)-1)]
#	print bins

	bins_smooth = [[rangeofindeces_partitions_smooth[i][1],rangeofindeces_partitions_smooth[i+1][0]] for i in range(len(rangeofindeces_partitions_smooth)-1)]
#	print bins_smooth


	maxtrain_accs_inbins = [mean(trainvals[a:b+1]) for [a,b] in bins]
	print maxtrain_accs_inbins

	maxtrain_accs_inbins_smooth = [mean(trainvals_smooth[a:b+1]) for [a,b] in bins_smooth]
	print maxtrain_accs_inbins_smooth



	average_maxtrain_accs_inbins.append(mean(maxtrain_accs_inbins))
	average_maxtrain_accs_inbins_smooth.append(mean(maxtrain_accs_inbins_smooth))

	f1 = figure(1)

	plot(trainvals,'y',label='train')
	plot(testvals,'m',label='test')


	plot(range(n-n1,n),trainvals_smooth,'g',label='train_smooth',linewidth=2)
	plot(range(n-n1,n),testvals_smooth,'r',label='test_smooth',linewidth=2)

	f1.set_figwidth(10)
	yticks(arange(floor(min(min(trainvals),min(testvals))),ceil(max(max(trainvals),max(testvals)))+5,5),arange(floor(min(min(trainvals),min(testvals))),ceil(max(max(trainvals),max(testvals)))+5,5))
	ylim(floor(min(min(trainvals),min(testvals)))-5,100)
	grid()


	for bin in bins_smooth:
	    fill_between(range(bin[0]+n-n1,bin[1]+1+n-n1),ceil(max(max(trainvals_smooth),max(testvals_smooth))),0,alpha=0.3,color='y')
	    text(bin[0]+n-n1,ceil(max(max(trainvals_smooth),max(testvals_smooth)))+3,'%.2f/%.2f\n%.2f/%.2f' % (mean(trainvals_smooth[bin[0]:bin[1]+1]),max(trainvals_smooth[bin[0]:bin[1]+1]),mean(testvals_smooth[bin[0]:bin[1]+1]),max(testvals_smooth[bin[0]:bin[1]+1])),fontsize=5)


	title ('Max train = %f/%d, max train_smooth = %f/%d\nmax test = %f/%d, max test_smooth = %f/%d' % (maxtrain[1],maxtrain[0],maxtrain_smooth[1],maxtrain_smooth[0],maxtest[1],maxtest[0],maxtest_smooth[1],maxtest_smooth[0]))
	pdf.savefig(dpi=500,bbox_inches='tight')
	close()


    f1 = figure(2)
#   for ind in xrange(len(alltopindeces)):
#	for val in alltopindeces[ind]:
#	    print val,len(trainvals),len(testvals)

#	    plot(ind,trainvals[val],'ko')
    pdf.savefig(dpi=500,bbox_inches='tight')
    close()


    f1 = figure(3)
#   for ind in xrange(len(alltopindeces_smooth)):
#	for val in alltopindeces_smooth[ind]:
#	    print val,len(trainvals_smooth),len(testvals_smooth)
#	    plot(ind,trainvals_smooth[val],'ko')
    pdf.savefig(dpi=500,bbox_inches='tight')
    close()

    f1 = figure(4)
    pdf.savefig(dpi=500,bbox_inches='tight')
    close()

    f1 = figure(5)
    pdf.savefig(dpi=500,bbox_inches='tight')
    close()





    pdf.close()

    print 'Avg/Median max test = %f/%f\nAvg/Median max test smoothed = %f/%f\nAvg/Median/Stdev opt train acc = %f/%f/%f\nAvg/Median/Stdev opt train acc smooth = %f/%f/%f' % (mean(all_maxtest),median(all_maxtest),mean(all_maxtest_smooth),median(all_maxtest_smooth),mean(average_maxtrain_accs_inbins),median(average_maxtrain_accs_inbins),std(average_maxtrain_accs_inbins),mean(average_maxtrain_accs_inbins_smooth),median(average_maxtrain_accs_inbins_smooth),std(average_maxtrain_accs_inbins_smooth))


if __name__ == '__main__':
    main(sys.argv[1:])




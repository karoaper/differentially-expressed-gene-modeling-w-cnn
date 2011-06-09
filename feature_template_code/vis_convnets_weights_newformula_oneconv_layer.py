import matplotlib
matplotlib.use("Agg")
from pylab import *
from numpy import *
import sys
import rpy2.robjects as r
from copy import copy
import os
from matplotlib.backends.backend_pdf import PdfPages
import mpl_toolkits.mplot3d.axes3d as p3





def plot_bars(features,feature_names,fn,maxval,show_center=True,normalize_scales=True):
    N = len(features)
    f1 = figure(1)
    feature_names = [str(x) for x in feature_names]

    for ind,onefeature in enumerate(features):

	ax = f1.add_subplot(N,1,ind+1)
	N4 = len(onefeature)

	fill(arange(N4+2),[0]+onefeature+[0],label=feature_names[ind])

	if normalize_scales:
	    ylim(-maxval,maxval)

	ymin, ymax = ylim()
	if show_center:
	    plot([N4/2,N4/2],[ymax,ymin],'k',linewidth=4)

	xticks(list(arange(0,N4,2)),[str((x-N4/2)*0.4) + 'Kb' for x in xrange(0,N4,2)])

	xlim(0,N4)
	for tick in ax.xaxis.get_major_ticks():
	    tick.label1.set_fontsize(17)
	ylabel(feature_names[ind].replace('_',' '),fontsize=30)

	grid()

    f1.set_figwidth(80)
    f1.set_figheight(60)
    savefig(fn)
    close()



def mycorrcoef(data,mattype='cov'):

    if mattype=='corr':
	return corrcoef(data),-1,1,'RdBu_r'
    elif mattype=='cov':
	return cov(data),-abs(cov(data)).max(),abs(cov(data)).max(),'RdBu_r'
    elif mattype=='angle':
	veclengths = [sqrt(dot(a,a)) for a in data]
	newdata = [array(data[i])/veclengths[i] for i in xrange(len(data))]

	x = pi/4.0-arccos(array([[min(1,dot(a,b)) for a in newdata] for b in newdata]))
	return x, -abs(x).max(),abs(x).max(),'RdBu_r'
    elif mattype=='dot':
	x = array([[dot(a,b) for a in data] for b in data])
	return x, 0,abs(x).max(),'Reds'




def analyse_top_layer_states(states,weights,output_fn):

    r.r('''source("/Users/khovsep/Dropbox/eblearn/demos/gene_expression/src/plot_heatmap.R")''')
    do_transform = r.r['do_pcatransform']


    toplayerdata = []
    classes = []

    n = len(states[0][0][5])

    outputfile = open('%s_toplayer.tab' % output_fn,'w')

    for ind in xrange(n):
	print >> outputfile, 'inp%d\t' % ind,

    print >> outputfile, 'out'


    for ind in xrange(n):
	print >> outputfile, 'continuous\t',

    print >> outputfile, 'binary'
    for ind in xrange(n):
	print >> outputfile, '\t',

    print >> outputfile, 'class'




    for exampleind,onestate in enumerate(states[0]+states[1]):
	classes.append(onestate[1])

	toplayerdata.append(onestate[5])

	thestate = [onestate[5][x][0] for x in xrange(n)]

	for val in thestate:
	    print >> outputfile, '%g\t' % val,
	print >> outputfile, '%d' % onestate[1]


    outputfile.close()



def computed_combinedfeatures(features,features_used,labels,top0,top1,N,binsize=1,combinetype='dontmiddiff',filterlim=0):

    N0 = len(top0)
    N1 = len(top1)



##  if combinetype == 'ignore_topweights_usemiddiffs':
##	combinedfeatures0 = [[sum(features[top0[0]][0][j][k*binsize:(k+1)*binsize])*1.0/len(features[top0[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[top0[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top0[1:]:
##	    binnedfeatures = [[sum(features[i][0][j][k*binsize:(k+1)*binsize])*1.0/len(features[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures0 = [[combinedfeatures0[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##	combinedfeatures1 = [[sum(features[top1[0]][0][j][k*binsize:(k+1)*binsize])*1.0/len(features[top1[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[top1[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top1[1:]:
##	    binnedfeatures = [[sum(features[i][0][j][k*binsize:(k+1)*binsize])*1.0/len(features[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures1 = [[combinedfeatures1[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##	combinedfeatures0 = [[combinedfeatures0[j][k]*1.0/N0 for k in xrange(len(combinedfeatures0[j]))] for j in xrange(N)]
##	combinedfeatures1 = [[combinedfeatures1[j][k]*1.0/N1 for k in xrange(len(combinedfeatures1[j]))] for j in xrange(N)]
##
##	combinedfeatures0 = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures0]
##	combinedfeatures1 = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures1]
##
##	num_features = size(array(combinedfeatures0),0)
##	size_kernel = size(array(combinedfeatures0),1)
##
##	combinedfeatures_diff = [[combinedfeatures0[j][i]-combinedfeatures1[j][i] for i in xrange(size_kernel)] for j in xrange(num_features)]
##	combinedfeatures_diff_abs = [[abs(x) if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff]
##	combinedfeatures_diff = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff]
##
##
##	# USED PART
##
##	combinedfeatures0_used = [[sum(features_used[top0[0]][0][j][k*binsize:(k+1)*binsize])*1.0/len(features_used[top0[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[top0[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top0[1:]:
##	    binnedfeatures = [[sum(features_used[i][0][j][k*binsize:(k+1)*binsize])*1.0/len(features_used[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures0_used = [[combinedfeatures0_used[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##	combinedfeatures1_used = [[sum(features_used[top1[0]][0][j][k*binsize:(k+1)*binsize])*1.0/len(features_used[top1[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[top1[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top1[1:]:
##	    binnedfeatures = [[sum(features_used[i][0][j][k*binsize:(k+1)*binsize])*1.0/len(features_used[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures1_used = [[combinedfeatures1_used[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##	combinedfeatures0_used = [[combinedfeatures0_used[j][k]*1.0/N0 for k in xrange(len(combinedfeatures0_used[j]))] for j in xrange(N)]
##	combinedfeatures1_used = [[combinedfeatures1_used[j][k]*1.0/N1 for k in xrange(len(combinedfeatures1_used[j]))] for j in xrange(N)]
##
##	combinedfeatures0_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures0_used]
##	combinedfeatures1_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures1_used]
##
##	num_features = size(array(combinedfeatures0_used),0)
##	size_kernel = size(array(combinedfeatures0_used),1)
##
##	combinedfeatures_diff_used = [[combinedfeatures0_used[j][i]-combinedfeatures1_used[j][i] for i in xrange(size_kernel)] for j in xrange(num_features)]
##
##
##	combinedfeatures_diff_abs_used = [[abs(x) if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff_used]
##	combinedfeatures_diff_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff_used]
##
##
##  elif combinetype == 'ignore_topweights_usemiddiffs_plus':
##
##	scaleweight = sign([x[0] for x in labels if x[1] == top0[0]][0])
##
##	combinedfeatures0 = [[sum(features[top0[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[top0[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[top0[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top0[1:]:
##	    scaleweight = sign([x[0] for x in labels if x[1] == i][0])
##
##	    binnedfeatures = [[sum(features[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures0 = [[combinedfeatures0[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##	scaleweight = sign([x[0] for x in labels if x[1] == top1[0]][0])
##
##	combinedfeatures1 = [[sum(features[top1[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[top1[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[top1[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top1[1:]:
##	    scaleweight = sign([x[0] for x in labels if x[1] == i][0])
##
##	    binnedfeatures = [[sum(features[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures1 = [[combinedfeatures1[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##
##	combinedfeatures0 = [[combinedfeatures0[j][k]*1.0/N0 for k in xrange(len(combinedfeatures0[j]))] for j in xrange(N)]
##	combinedfeatures1 = [[combinedfeatures1[j][k]*1.0/N1 for k in xrange(len(combinedfeatures1[j]))] for j in xrange(N)]
##
##	combinedfeatures0 = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures0]
##	combinedfeatures1 = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures1]
##
##
##	num_features = size(array(combinedfeatures0),0)
##	size_kernel = size(array(combinedfeatures0),1)
##
##
##	combinedfeatures_diff = [[combinedfeatures0[j][i]-combinedfeatures1[j][i] for i in xrange(size_kernel)] for j in xrange(num_features)]
##	combinedfeatures_diff_abs = [[abs(x) if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff]
##	combinedfeatures_diff = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff]
##
##
##	# USED PART
##	scaleweight = sign([x[0] for x in labels if x[1] == top0[0]][0])
##
##	combinedfeatures0_used = [[sum(features_used[top0[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[top0[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[top0[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top0[1:]:
##	    scaleweight = sign([x[0] for x in labels if x[1] == i][0])
##
##	    binnedfeatures = [[sum(features_used[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures0_used = [[combinedfeatures0_used[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##	scaleweight = sign([x[0] for x in labels if x[1] == top1[0]][0])
##
##	combinedfeatures1_used = [[sum(features_used[top1[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[top1[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[top1[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top1[1:]:
##	    scaleweight = sign([x[0] for x in labels if x[1] == i][0])
##
##	    binnedfeatures = [[sum(features_used[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures1_used = [[combinedfeatures1_used[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##	combinedfeatures0_used = [[combinedfeatures0_used[j][k]*1.0/N0 for k in xrange(len(combinedfeatures0_used[j]))] for j in xrange(N)]
##	combinedfeatures1_used = [[combinedfeatures1_used[j][k]*1.0/N1 for k in xrange(len(combinedfeatures1_used[j]))] for j in xrange(N)]
##
##	combinedfeatures0_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures0_used]
##	combinedfeatures1_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures1_used]
##
##	num_features = size(array(combinedfeatures0_used),0)
##	size_kernel = size(array(combinedfeatures0_used),1)
##
##	combinedfeatures_diff_used = [[combinedfeatures0_used[j][i]-combinedfeatures1_used[j][i] for i in xrange(size_kernel)] for j in xrange(num_features)]
##
##	combinedfeatures_diff_abs_used = [[abs(x) if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff_used]
##	combinedfeatures_diff_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff_used]
##
##
##  elif combinetype == 'use_topweights_usemiddiffs':
##
##	scaleweight = abs([x[0] for x in labels if x[1] == top0[0]][0])
##
##	combinedfeatures0 = [[sum(features[top0[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[top0[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[top0[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top0[1:]:
##	    scaleweight = abs([x[0] for x in labels if x[1] == i][0])
##
##	    binnedfeatures = [[sum(features[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures0 = [[combinedfeatures0[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##	scaleweight = abs([x[0] for x in labels if x[1] == top1[0]][0])
##
##	combinedfeatures1 = [[sum(features[top1[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[top1[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[top1[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top1[1:]:
##	    scaleweight = abs([x[0] for x in labels if x[1] == i][0])
##
##	    binnedfeatures = [[sum(features[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures1 = [[combinedfeatures1[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##
##	combinedfeatures0 = [[combinedfeatures0[j][k]*1.0/N0 for k in xrange(len(combinedfeatures0[j]))] for j in xrange(N)]
##	combinedfeatures1 = [[combinedfeatures1[j][k]*1.0/N1 for k in xrange(len(combinedfeatures1[j]))] for j in xrange(N)]
##
##	combinedfeatures0 = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures0]
##	combinedfeatures1 = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures1]
##
##
##	num_features = size(array(combinedfeatures0),0)
##	size_kernel = size(array(combinedfeatures0),1)
##
##
##	combinedfeatures_diff = [[combinedfeatures0[j][i]-combinedfeatures1[j][i] for i in xrange(size_kernel)] for j in xrange(num_features)]
##	combinedfeatures_diff_abs = [[abs(x) if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff]
##	combinedfeatures_diff = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff]
##
##
##	# USED PART
##
##
##	scaleweight = abs([x[0] for x in labels if x[1] == top0[0]][0])
##
##	combinedfeatures0_used = [[sum(features_used[top0[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[top0[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[top0[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top0[1:]:
##	    scaleweight = abs([x[0] for x in labels if x[1] == i][0])
##
##	    binnedfeatures = [[sum(features_used[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures0_used = [[combinedfeatures0_used[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##	scaleweight = abs([x[0] for x in labels if x[1] == top1[0]][0])
##
##	combinedfeatures1_used = [[sum(features_used[top1[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[top1[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[top1[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top1[1:]:
##	    scaleweight = abs([x[0] for x in labels if x[1] == i][0])
##
##	    binnedfeatures = [[sum(features_used[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures1_used = [[combinedfeatures1_used[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##	combinedfeatures0_used = [[combinedfeatures0_used[j][k]*1.0/N0 for k in xrange(len(combinedfeatures0_used[j]))] for j in xrange(N)]
##	combinedfeatures1_used = [[combinedfeatures1_used[j][k]*1.0/N1 for k in xrange(len(combinedfeatures1_used[j]))] for j in xrange(N)]
##
##	combinedfeatures0_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures0_used]
##	combinedfeatures1_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures1_used]
##
##	num_features = size(array(combinedfeatures0_used),0)
##	size_kernel = size(array(combinedfeatures0_used),1)
##
##	combinedfeatures_diff_used = [[combinedfeatures0_used[j][i]-combinedfeatures1_used[j][i] for i in xrange(size_kernel)] for j in xrange(num_features)]
##
##	combinedfeatures_diff_abs_used = [[abs(x) if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff_used]
##	combinedfeatures_diff_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff_used]
##
##
##  elif combinetype == 'use_topweights_usemiddiffs_noabs':
##
##	scaleweight = [x[0] for x in labels if x[1] == top0[0]][0]
##
##	combinedfeatures0 = [[sum(features[top0[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[top0[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[top0[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top0[1:]:
##	    scaleweight = [x[0] for x in labels if x[1] == i][0]
##
##	    binnedfeatures = [[sum(features[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures0 = [[combinedfeatures0[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##	scaleweight = [x[0] for x in labels if x[1] == top1[0]][0]
##
##	combinedfeatures1 = [[sum(features[top1[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[top1[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[top1[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top1[1:]:
##	    scaleweight = [x[0] for x in labels if x[1] == i][0]
##
##	    binnedfeatures = [[sum(features[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures1 = [[combinedfeatures1[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##
##	combinedfeatures0 = [[combinedfeatures0[j][k]*1.0/N0 for k in xrange(len(combinedfeatures0[j]))] for j in xrange(N)]
##	combinedfeatures1 = [[combinedfeatures1[j][k]*1.0/N1 for k in xrange(len(combinedfeatures1[j]))] for j in xrange(N)]
##
##	combinedfeatures0 = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures0]
##	combinedfeatures1 = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures1]
##
##
##	num_features = size(array(combinedfeatures0),0)
##	size_kernel = size(array(combinedfeatures0),1)
##
##
##	combinedfeatures_diff = [[combinedfeatures0[j][i]-combinedfeatures1[j][i] for i in xrange(size_kernel)] for j in xrange(num_features)]
##	combinedfeatures_diff_abs = [[abs(x) if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff]
##	combinedfeatures_diff = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff]
##
##
##	# USED PART
##	scaleweight = [x[0] for x in labels if x[1] == top0[0]][0]
##
##	combinedfeatures0_used = [[sum(features_used[top0[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[top0[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[top0[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top0[1:]:
##	    scaleweight = [x[0] for x in labels if x[1] == i][0]
##
##	    binnedfeatures = [[sum(features_used[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures0_used = [[combinedfeatures0_used[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##
##	scaleweight = [x[0] for x in labels if x[1] == top1[0]][0]
##
##	combinedfeatures1_used = [[sum(features_used[top1[0]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[top1[0]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[top1[0]][0][j])*1.0/binsize))] for j in xrange(N)]
##	for i in top1[1:]:
##	    scaleweight = [x[0] for x in labels if x[1] == i][0]
##
##	    binnedfeatures = [[sum(features_used[i][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[i][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[i][0][j])*1.0/binsize))] for j in xrange(N)]
##	    combinedfeatures1_used = [[combinedfeatures1_used[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]
##
##	combinedfeatures0_used = [[combinedfeatures0_used[j][k]*1.0/N0 for k in xrange(len(combinedfeatures0_used[j]))] for j in xrange(N)]
##	combinedfeatures1_used = [[combinedfeatures1_used[j][k]*1.0/N1 for k in xrange(len(combinedfeatures1_used[j]))] for j in xrange(N)]
##
##	combinedfeatures0_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures0_used]
##	combinedfeatures1_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures1_used]
##
##	num_features = size(array(combinedfeatures0_used),0)
##	size_kernel = size(array(combinedfeatures0_used),1)
##
##	combinedfeatures_diff_used = [[combinedfeatures0_used[j][i]-combinedfeatures1_used[j][i] for i in xrange(size_kernel)] for j in xrange(num_features)]
##
##	combinedfeatures_diff_abs_used = [[abs(x) if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff_used]
##	combinedfeatures_diff_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff_used]
##
##
##  else:

    if True:
	print labels
	scaleweight = labels[0][0]
	combinedfeatures_diff = [[sum(features[labels[0][1]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[labels[0][1]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[labels[0][1]][0][j])*1.0/binsize))] for j in xrange(N)]
#	for someind in xrange(10):
#	    scaleweight = labels[8][0]
#	    binnedfeatures = [[sum(features[labels[8][1]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[labels[8][1]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[labels[8][1]][0][j])*1.0/binsize))] for j in xrange(N)]
#	    combinedfeatures_diff = [[combinedfeatures_diff[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]


	for label in labels[1:]:
	    scaleweight = label[0]
	    binnedfeatures = [[sum(features[label[1]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features[label[1]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features[label[1]][0][j])*1.0/binsize))] for j in xrange(N)]
	    combinedfeatures_diff = [[combinedfeatures_diff[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]

	num_features = size(array(combinedfeatures_diff),0)
	size_kernel = size(array(combinedfeatures_diff),1)


	combinedfeatures_diff_abs = [[abs(x) if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff]
	combinedfeatures_diff = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff]



	# USED PART

	scaleweight = labels[0][0]
	combinedfeatures_diff_used = [[sum(features_used[labels[0][1]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[labels[0][1]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[labels[0][1]][0][j])*1.0/binsize))] for j in xrange(N)]
#	for someind in xrange(10):
#	    scaleweight = labels[8][0]
#	    binnedfeatures = [[sum(features_used[labels[8][1]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[labels[8][1]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[labels[8][1]][0][j])*1.0/binsize))] for j in xrange(N)]
#	    combinedfeatures_diff_used = [[combinedfeatures_diff_used[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]


	for label in labels[1:]:
	    scaleweight = label[0]
	    binnedfeatures = [[sum(features_used[label[1]][0][j][k*binsize:(k+1)*binsize])*scaleweight*1.0/len(features_used[label[1]][0][j][k*binsize:(k+1)*binsize]) for k in xrange(int(len(features_used[label[1]][0][j])*1.0/binsize))] for j in xrange(N)]
	    combinedfeatures_diff_used = [[combinedfeatures_diff_used[j][k]+binnedfeatures[j][k] for k in xrange(len(binnedfeatures[j]))] for j in xrange(N)]

	num_features = size(array(combinedfeatures_diff_used),0)
	size_kernel = size(array(combinedfeatures_diff_used),1)


	combinedfeatures_diff_abs_used = [[abs(x) if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff_used]
	combinedfeatures_diff_used = [[x if abs(x) > filterlim else 0 for x in y] for y in combinedfeatures_diff_used]





    return combinedfeatures_diff,combinedfeatures_diff_abs,combinedfeatures_diff_used,combinedfeatures_diff_abs_used






def computed_separate_combinedfeatures(features,labels,N,indstoget,filterlim=0):

    combinedfeatures = []
#   print indstoget
    for labelind in indstoget:

	scaleweight = [x[0] for x in labels if x[1] == labelind][0]

	combinedfeature = features[labelind][0]
#	print len(combinedfeature)
#	print combinedfeature

#	num_features = len(combinedfeature)
#	size_kernel = size(array(combinedfeature),1)

	combinedfeature = [[scaleweight*x if abs(scaleweight*x) > filterlim else 0 for x in combinedfeature[y]] for y in combinedfeature]


	combinedfeatures.append([scaleweight,labelind,combinedfeature])



    return combinedfeatures




#def analyze_top2layers(clustfeats,whichclass,topunitcluster,allfeatures_copy,topxpercent,combinetype,numoffilters,groupsize,outputgroupsize,binsize,outputfn,filterlim,subsample1,pixelsize,plot_top2_analysis=False,plot_all_features=False,input_feature_names=[],states=[]):
#   analyze_top2layers(whichclass,topunitcluster,features,topxpercent_top,    numoffilters,groupsize,outputgroupsize,output_fn_loc,subsample1,pixelsize,plot_top2_analysis,plot_all_features,input_feature_names,newstates)
def analyze_top2layers(whichclass,topunitcluster,allfeatures_copy,topxpercent,numoffilters,groupsize,outputgroupsize,outputfn,	   subsample1,pixelsize,keptpeakmovepercent,keptpeakheightpercent,whichfeaturespeaks,plot_top2_analysis=False,plot_all_features=False,input_feature_names=[],states=[],allmeans=[],allstds=[]):


    r.r('''source("/Users/khovsep/Dropbox/eblearn/demos/gene_expression/src/plot_heatmap.R")''')
    drawheatmaps = r.r['draw_heatmaps']
    qtclust = r.r['do_qtclust']

    allfeatures = allfeatures_copy


    numtracks = len(input_feature_names)

    features = allfeatures[0]

    N = len(features[0][0][0])


    if topunitcluster == 'all':
	topunitcluster = range(N)

    if whichclass == 'both':
	labels = [[features[0][0][0][i]-features[1][0][0][i],i] for i in topunitcluster]
    else:
	labels = [[features[whichclass][0][0][i],i] for i in topunitcluster]


    labels_abs = [[abs(labels[x][0]),labels[x][1]] for x in xrange(len(labels))]


    labels_abs.sort(reverse=True)
    labels.sort()


    ## extract the top up and top down penultimate units

    top0 = [x[1] for x in labels[int(topxpercent*N):] if x[0] > 0]
    top1 = [x[1] for x in labels[:int((1-topxpercent)*N)] if x[0] < 0]

    feature = array([[x[0] if x[1] in top0+top1 else 0 for x in labels]])

    num_features = size(feature,0)
    size_kernel = size(feature,1)

    ## plot top layer as well the states



    N = len(allfeatures[1][0][0])
    N1 = len(allfeatures[1][0][0][0])


    features = dict([(x,allfeatures[1][x]) for x in allfeatures[1] if x in topunitcluster])

    input_features_names_used = [[x+'_' + str(s) for s in range(1,numoffilters+1)] for x in input_feature_names]
    input_features_names_used = [item for sublist in input_features_names_used for item in sublist]


    ## define the various values mapping the second layer from the top to the actual window


    sizeofinputfeatures = len(allfeatures[-1][0][0][0]) #size of filter

    sizeofeachunit = (subsample1-1)*pixelsize+sizeofinputfeatures*pixelsize
    totallengthofwin = N1*subsample1-1 + sizeofinputfeatures
    totallengthofwininbases = (N1-1)*subsample1*pixelsize + sizeofeachunit


    TSSinunits = round((totallengthofwin/2-sizeofeachunit)/(subsample1*pixelsize))
    realTSS = TSSinunits*subsample1*pixelsize+sizeofeachunit





    #define an expanded form of each feature in the second layer from the top, with the expanded view covering base-by-base the entire window

    features_expanded = {}


    for key in features:
	if not key in topunitcluster:
	    continue

	features_expanded[key] = [{},0]
	ind = 0

	for n in xrange(N1):
	    for s in xrange(subsample1):
		for j in xrange(sizeofinputfeatures):
		    for key1 in features[key][0]:
			if features_expanded[key][0].has_key(key1):
			    features_expanded[key][0][key1][ind+j] += features[key][0][key1][n]
			else:
			    features_expanded[key][0][key1] = zeros(totallengthofwin)
			    features_expanded[key][0][key1][ind+j] += features[key][0][key1][n]

		ind += 1



    use_expanded_features = True

    if use_expanded_features:
	features_used = features_expanded
	N1 =	totallengthofwin
    else:
	features_used = features








    xaxis_ticks = list(sorted(-array([x for x in arange(0,totallengthofwin/2+1,1) if x*pixelsize % outputgroupsize == 0])+totallengthofwin/2)) +list(array([x for x in arange(1,totallengthofwin/2+1,1) if x*pixelsize % outputgroupsize == 0])+totallengthofwin/2)
    xaxis_labels = ['%s' % ('Peak Center' if x == 0 else str(x*pixelsize/1000.0) + 'Kb' if abs(x*pixelsize) >= 1000 and abs(x*pixelsize) < 1000000 else (str(x*pixelsize/1000000.0)+'Mb' if abs(x*pixelsize) >= 1000000 else str(x*pixelsize)+'b')) for x in list(array([x for x in arange(0,totallengthofwin/2+1,1) if x*pixelsize % outputgroupsize == 0])-totallengthofwin/2) +list(array([x for x in arange(1,totallengthofwin/2+1,1) if x*pixelsize % outputgroupsize == 0]))]

    xaxis_labels_heatmap = ['%s' % ('' if x*pixelsize % outputgroupsize <> 0 else 'TSS---->' if x == 0 else str(x*pixelsize/1000.0) + 'Kb' if abs(x*pixelsize) >= 1000 and abs(x*pixelsize) < 1000000 else (str(x*pixelsize/1000000.0)+'Mb' if abs(x*pixelsize) >= 1000000 else str(x*pixelsize)+'b')) for x in list(sorted(array([x for x in arange(0,totallengthofwin/2+1,1)])-totallengthofwin/2)) +list(array([x for x in arange(1,totallengthofwin/2-1,1)]))+[totallengthofwin/2]]


    features_combinedtracks_abs = {}
    for key in features_used:
	features_combinedtracks_abs[key] = [{},0]
	ind = 0
	keyind = 0
	for key1 in features_used[key][0]:
	    if ind == 0:
		features_combinedtracks_abs[key][0][keyind] = [abs(x) for x in features_used[key][0][key1]]
	    else:
		features_combinedtracks_abs[key][0][keyind] = [abs(features_used[key][0][key1][x])+features_combinedtracks_abs[key][0][keyind][x] for x in xrange(len(features_used[key][0][key1]))]

	    ind += 1
	    if ind % numoffilters == 0:
		ind = 0
		keyind += 1



    features_maxvals = zeros([N,len(features_used)])
    for ind,key in enumerate(features_used):
	for key1 in features_used[key][0]:
	    features_maxvals[key1,ind] = max([abs(x) for x in features_used[key][0][key1]])

    drawheatmaps(False,True,r.r.matrix(r.r.c(*reshape(features_maxvals,len(features_used)*N,1)), nrow=N),4000,r.r.c(*input_features_names_used),'%s_toplayer_maxabsfeatures_class%s.pdf' % (outputfn,str(whichclass)),r.r.c(*features_used.keys()))




    groupstarts = list(sorted(-array([x for x in arange(0,totallengthofwin/2+1,1) if x * pixelsize % groupsize == 0])+totallengthofwin/2)) +list(array([x for x in arange(1,totallengthofwin/2,1) if x*pixelsize % groupsize == 0])+totallengthofwin/2)

    if groupstarts[0] <> 0:
	groupstarts[0] = 0


    features_maxvals = zeros([len(features_used),N*len(groupstarts)])


    for ind1,key in enumerate(features_used):
	ind = 0
	for groupind,groupstart in enumerate(groupstarts):
	    groupend = groupstarts[groupind+1] if groupind != len(groupstarts)-1 else N1

	    for key1 in features_used[key][0]:
		features_maxvals[ind1,ind] = max([abs(x) for x in features_used[key][0][key1][groupstart:groupend]])
		ind += 1

    drawheatmaps(True,False,r.r.matrix(r.r.c(*reshape(features_maxvals,len(features_used)*N*len(groupstarts),1)), nrow=len(features_used)),4000,r.r.c(*features_used.keys()),'%s_toplayer_maxabsfeatures_bygroups_class%s.pdf' % (outputfn,str(whichclass)))



    combinedfeatures_diff,combinedfeatures_diff_abs,combinedfeatures_diff_used,combinedfeatures_diff_abs_used = computed_combinedfeatures(features,features_used,labels,top0,top1,N)#binsize,combinetype,filterlim,N)
    individual_combinedfeatures = computed_separate_combinedfeatures(features,labels,N,top0+top1)


    alldistinctpeaks = dict([(featureind,[[],[]]) for featureind in whichfeaturespeaks])
    for indiv_featureind,indiv_feature in enumerate(individual_combinedfeatures):
	entirematrix_wholeouterfeature_curves, allinnerfeaturekeys,N4 = create_combined_template(allfeatures[2],indiv_feature[2],subsample1)
	feature = [[sum(entirematrix_wholeouterfeature_curves[key][ind]) for ind in xrange(totallengthofwin)] for key in allinnerfeaturekeys]

	alldistinctpeak = findpeaks_nofilters(feature,outputfn,whichclass,input_feature_names,xaxis_ticks,xaxis_labels,allinnerfeaturekeys,whichfeaturespeaks)
	for featureind in alldistinctpeak:
	    alldistinctpeaks[featureind][0].extend(alldistinctpeak[featureind][0])
	    alldistinctpeaks[featureind][1].extend(alldistinctpeak[featureind][1])




	maxval = abs(array(feature)).max()

	if False:


	    f1 = figure(1)
	    for ind in xrange(len(feature)):


		onefeature = feature[ind]


		ax = f1.add_subplot(len(feature),1,ind+1)
		ylabel(input_feature_names[allinnerfeaturekeys[ind]].replace('_',' '),fontsize=40,rotation='horizontal')

		ax.fill_between(range(len(onefeature)),0,list(onefeature),alpha=0.3,label='template',color='b')
		if alldistinctpeak.has_key(ind):
		    for peak in alldistinctpeak[ind][0]:
			ax.plot([peak[1],peak[1]],[-maxval,maxval],'r')
		    for peak in alldistinctpeak[ind][1]:
			ax.plot([peak[1],peak[1]],[-maxval,maxval],'r')


		ylim(-maxval,maxval)
		ymin, ymax = ylim()
		plot([totallengthofwin/2,totallengthofwin/2],[ymin,ymax],'k',linewidth=2)
		plot([0,totallengthofwin-1],[0,0],'k',linewidth=1)

		xticks(xaxis_ticks,xaxis_labels)

		locs,labels = yticks()
		yticks(locs, map(lambda x: "%g" % x, locs))


		xlim(0,totallengthofwin-1)

		ax.yaxis.grid(True)
		ax.xaxis.grid(True)

		for tick in ax.get_yticklabels():
		    tick.set_color('b')
		    tick.set_fontsize(10)


		for tick in ax.get_xticklabels():
		    tick.set_fontsize(30)



	    f1.set_figwidth(40)
	    f1.set_figheight(60)


	    savefig('%s_curves_entirefeature_overlapping_feature%d_class%s.pdf' % (outputfn,indiv_feature[1],str(whichclass)))
	    close()



    f1 = figure()
    axs = []


    maxy = 0
    for featureind in alldistinctpeaks:
	alldistinctpeaks[featureind][0].sort(reverse=True)
	alldistinctpeaks[featureind][1].sort(reverse=True)

	alldistinctpeaks[featureind] = alldistinctpeaks[featureind][0][:int(len(alldistinctpeaks[featureind][0])*keptpeakheightpercent)] + alldistinctpeaks[featureind][1][:int(len(alldistinctpeaks[featureind][1])*keptpeakheightpercent)]
	ax = f1.add_subplot(len(alldistinctpeaks),1,featureind+1)
	Ys = {}
	for val in alldistinctpeaks[featureind]:
	    if not Ys.has_key(val[1]):
		Ys[val[1]] = 0

	    Ys[val[1]] += val[2]

	X = [x for x in Ys if Ys[x]>0]
	Y = [Ys[x] for x in Ys if Ys[x]>0]
	maxy = max(maxy,max(Y))


#	Y = [x[2] for x in alldistinctpeaks[featureind] if x[2]>0]
	p1 = bar(X,Y,width=1.0,color='r',edgecolor='r',label='ENHANCER features')
#	p1 = bar(X,Y,width=1.0,color='r',edgecolor='r',label='UP features')


	X = [x for x in Ys if Ys[x]<0]
	Y = [-Ys[x] for x in Ys if Ys[x]<0]
	maxy = max(maxy,max(Y))



#	X = [x[1] for x in alldistinctpeaks[featureind] if x[2]<0]
#	Y = [-x[2] for x in alldistinctpeaks[featureind] if x[2]<0]
	p2 = bar(X,Y,width=1.0,color='b',edgecolor='b',label='SILENCER features')
#	p2 = bar(X,Y,width=1.0,color='b',edgecolor='b',label='DOWN features')

#	ylabel('Feature Significance',fontsize=10)
	ylabel(input_feature_names[allinnerfeaturekeys[featureind]].replace('_',' '),fontsize=10,rotation='vertical')

	if featureind == 0:
	    ax.legend(loc=2,markerscale=0.1,frameon=False)
	    leg = gca().get_legend()
	    ltext  = leg.get_texts()  # all the text.Text instance in the legend
	    llines = leg.get_lines()  # all the lines.Line2D instance in the legend
	    frame  = leg.get_frame()  # the patch.Rectangle instance surrounding the legend

	    # see text.Text, lines.Line2D, and patches.Rectangle for more info on
	    # the settable properties of lines, text, and rectangles
    #	    frame.set_facecolor('0.80')	     # set the frame face color to light gray
	    setp(ltext, fontsize='medium')#xx-small')	 # the legend text fontsize
    #	    setp(llines, linewidth=1.5)	     # the legend linewidth


#	ylim(0,maxval)
#	ymin, ymax = ylim()
	xticks(xaxis_ticks,xaxis_labels)



	xlim(0,totallengthofwin-1)
	grid()
#	ax.yaxis.grid(True)
#	ax.xaxis.grid(True)

	for tick in ax.get_yticklabels():
	    tick.set_color('b')
	    tick.set_fontsize(10)


	for tick in ax.get_xticklabels():
	    tick.set_fontsize(9)
	axs.append(ax)


	if True:
	    ax1 = ax.twinx()
    #	    ax1.fill_between(range(len(allmeans[0][featureind])),list(allmeans[0][featureind]-allstds[0][featureind]),list(allmeans[0][featureind]+allstds[0][featureind]),alpha=1,color='r',linewidth=2,label='Median UP')

#	    p3 = ax1.plot(allmeans[0][featureind],'r',label='mean UP')
	    p3 = ax1.plot(allmeans[0][featureind],'r',label='mean ENHANCER')

    #	    xticks(xaxis_ticks,xaxis_labels)
    #	    grid()
    #	    for tick in ax1.get_yticklabels():
    #		tick.set_color('b')
    #		tick.set_fontsize(12)


    #	    for tick in ax1.get_xticklabels():
    #		tick.set_fontsize(15)

    #	    ax1 = ax.twinx()
#	    p4 = ax1.plot(allmeans[1][featureind],'b',label='mean DOWN')
	    p4 = ax1.plot(allmeans[1][featureind],'b',label='mean SILENCER')

    #	    ax1.fill_between(range(len(allmeans[1][featureind])),list(allmeans[1][featureind]-allstds[1][featureind]),list(allmeans[1][featureind]+allstds[1][featureind]),alpha=1,color='b',linewidth=2,label='Median DOWN')

    #	    ax1.plot([0],'r',label=
    #	    ax1.plot(allmeans[1][featureind]-,'b',label='DOWN')
	    xticks(xaxis_ticks,xaxis_labels)
	    grid()
	    for tick in ax1.get_yticklabels():
		tick.set_color('b')
		tick.set_fontsize(10)


	    for tick in ax1.get_xticklabels():
		tick.set_fontsize(10)

	    locs,labels = yticks()
	    yticks([locs[0],locs[-1]], map(lambda x: "%g" % x, [locs[0],locs[-1]]))


	    if featureind == 0:

		ax1.legend(loc=1,markerscale=0.1,frameon=False)
		leg = gca().get_legend()
		ltext  = leg.get_texts()  # all the text.Text instance in the legend
		llines = leg.get_lines()  # all the lines.Line2D instance in the legend
		frame  = leg.get_frame()  # the patch.Rectangle instance surrounding the legend

		# see text.Text, lines.Line2D, and patches.Rectangle for more info on
		# the settable properties of lines, text, and rectangles
	#	frame.set_facecolor('0.80')	 # set the frame face color to light gray
		setp(ltext, fontsize='medium')#xx-small')    # the legend text fontsize
	#	setp(llines, linewidth=1.5)	 # the legend linewidth


	#	legend([p1,p2,p3,p4],['ENHANCING features','SILENCING features','median UP','median DOWN'])




#	title('Significant %s features mined from gene classification model' % input_feature_names[allinnerfeaturekeys[featureind]],fontsize=10)


    for ax in axs:
	axes(ax)
	ylim(0,maxy)
	locs,labels = yticks()
	plot([totallengthofwin/2,totallengthofwin/2],[0,locs[-1]],'k',linewidth=2)
	yticks([locs[0],locs[-1]], map(lambda x: "%g" % x, [locs[0],locs[-1]]))




    f1.set_figwidth(10)
    f1.set_figheight(2*len(alldistinctpeaks))


    savefig('%s_significant_peaks_%g_class%s.pdf' % (outputfn,keptpeakheightpercent,str(whichclass)))
    close()




    store_peaks(alldistinctpeaks,len(feature[0]),pixelsize,outputfn,whichclass,input_feature_names,xaxis_ticks,xaxis_labels,allinnerfeaturekeys)













    if plot_top2_analysis:

	feature = array(combinedfeatures_diff_used)
	num_features = size(feature,0)
	size_kernel = size(feature,1)
	drawheatmaps(False,False,r.r.matrix(r.r.c(*reshape(feature,size_kernel*num_features,1)), nrow=num_features),4000,r.r.c(*input_features_names_used),'%s_2ndlayerfromtop_binned_class%s.pdf' % (outputfn,str(whichclass)),r.r.c(*xaxis_labels_heatmap))


	feature = array(combinedfeatures_diff_abs_used)
	num_features = size(feature,0)
	size_kernel = size(feature,1)
	drawheatmaps(False,False,r.r.matrix(r.r.c(*reshape(feature,size_kernel*num_features,1)), nrow=num_features),4000,r.r.c(*input_features_names_used),'%s_2ndlayerfromtop_binned_class%s_abs.pdf' % (outputfn,str(whichclass)),r.r.c(*xaxis_labels_heatmap))


	features_combinedtracks_diff_abs = []
	ind = 0
	for rowvals in combinedfeatures_diff_abs_used:
	    if ind == 0:
		features_combinedtracks_diff_abs.append(rowvals)
	    else:
		features_combinedtracks_diff_abs[-1] = [rowvals[x] + features_combinedtracks_diff_abs[-1][x] for x in xrange(len(rowvals))]

	    ind += 1
	    if ind % numoffilters == 0:
		ind = 0


	feature = array(features_combinedtracks_diff_abs)
	num_features = size(feature,0)
	size_kernel = size(feature,1)
	drawheatmaps(False,False,r.r.matrix(r.r.c(*reshape(feature,size_kernel*num_features,1)), nrow=num_features),4000,r.r.c(*input_feature_names),'%s_2ndlayerfromtop_binned_combinedtracks_class%s_abs.pdf' % (outputfn,str(whichclass)),r.r.c(*xaxis_labels_heatmap))



	f1 = figure(1)
	ax = f1.add_subplot(111)
	intermat,vmin,vmax,colormap = mycorrcoef(array(features_combinedtracks_diff_abs),'cov')
	ax.imshow(intermat,interpolation='nearest',cmap=get_cmap(colormap),vmin=vmin,vmax=vmax)
	xticks(arange(size(intermat,0)),arange(size(intermat,1)))
	yticks(arange(size(intermat,0)),input_feature_names)
	for tick in ax.get_yticklabels():
	    tick.set_fontsize(28)

	f1.set_figwidth(22)
	f1.set_figheight(20)

	savefig('%s_curves_layer2_covs_corrs_combinedtracks_class%s_abs.pdf' % (outputfn,str(whichclass)))
	close()





	corr_mats_combinedtracks_diff_abs_bygroups = {}
	cov_mats_combinedtracks_diff_abs_bygroups = {}
	dot_mats_combinedtracks_diff_abs_bygroups = {}
	angle_mats_combinedtracks_diff_abs_bygroups = {}


	N1 = len(combinedfeatures_diff_used[0])


	numsubplotcols = 4.0 if len(groupstarts) / 4.0 >= 4 else 3.0 if len(groupstarts) /3.0 >= 3 else 2.0

	pp = PdfPages('%s_curves_layer2_covs_corrs_combinedtracks_class%s_abs_bygroups.pdf' % (outputfn,str(whichclass)))
	f1 = figure(3)
	for groupind,groupstart in enumerate(groupstarts):
	    groupend = groupstarts[groupind+1] if groupind != len(groupstarts)-1 else N1

	    ax = f1.add_subplot(ceil(len(groupstarts)/numsubplotcols),numsubplotcols,groupind+1)
	    cov_mats_combinedtracks_diff_abs_bygroups[groupstart,groupend],vmin,vmax,colormap = mycorrcoef(array([rowvals[groupstart:groupend] for rowvals in features_combinedtracks_diff_abs]),'cov')
	    ax.imshow(cov_mats_combinedtracks_diff_abs_bygroups[groupstart,groupend],interpolation='nearest',cmap=get_cmap(colormap),vmin=vmin,vmax=vmax)
	    xticks(arange(numtracks),arange(numtracks))
	    yticks(arange(numtracks),input_feature_names)
	    for tick in ax.get_yticklabels():
		tick.set_fontsize(25)



	    title('%s -- %s' % (str((groupstart-N1/2)*pixelsize/1000.0) + 'Kb' if abs((groupstart-N1/2)*pixelsize) >= 1000 and abs((groupstart-N1/2)*pixelsize) < 1000000 else (str((groupstart-N1/2)*pixelsize/1000000.0)+'Mb' if abs((groupstart-N1/2)*pixelsize) >= 1000000 else str((groupstart-N1/2)*pixelsize)+'b'),str((groupend-N1/2)*pixelsize/1000.0) + 'Kb' if abs((groupend-N1/2)*pixelsize) >=1000 and abs((groupend-N1/2)*pixelsize) < 1000000 else (str((groupend-N1/2)*pixelsize/1000000.0)+'Mb' if abs((groupend-N1/2)*pixelsize) >= 1000000 else str((groupend-N1/2)*pixelsize)+'b')),fontsize=25)


	if numsubplotcols == 2:
	    f1.set_figwidth(40)
	    f1.set_figheight(30)
	elif numsubplotcols == 3:
	    f1.set_figwidth(32)
	    f1.set_figheight(30)
	elif numsubplotcols == 4:
	    f1.set_figwidth(40)
	    f1.set_figheight(50)


	pp.savefig()
	close()

	pp.close()




	feature_weights = [sum(x) for x in combinedfeatures_diff]


    return combinedfeatures_diff,alldistinctpeaks





def sigmoid(val):
    return 1./(1+exp(-val))





def change_weights(weights,bias,default_input,delta_input):
    newweights = {}

    default_output = bias
    for subfeatureind in weights:
	for valind in xrange(len(weights[subfeatureind])):
	    default_output += weights[subfeatureind][valind]*default_input


    for subfeatureind in weights:
	newweights[subfeatureind] = []
	for valind in xrange(len(weights[subfeatureind])):
	    newweights[subfeatureind].append(default_output + weights[subfeatureind][valind]*delta_input)


    default_output = sigmoid(default_output)

    for subfeatureind in weights:
	for valind in xrange(len(weights[subfeatureind])):
	    newweights[subfeatureind][valind] = log(sigmoid(newweights[subfeatureind][valind])/default_output)


    return newweights





def read_states(states_fn):
    states = {0:[],1:[]}
    layers = ['label','predicted label','input','convolution','subsamples','full','output']

    firstexample = True
    for line in open(states_fn,'r'):
	if line[0] == '=':
	    if not firstexample:
		if numoflines > 2 and correct:
		    states[newstatemap[0]].append([abs(newstatemap[-1][0][0]-newstatemap[-1][1][0])]+newstatemap)
	    firstexample = False
	    newstatemap = []
	    numoflines = 0
	    correct = True
	    continue

	numoflines += 1
	linesparts = line.rstrip().rsplit(': ')
	if linesparts[0] == 'label':
	    newstatemap.append(int(linesparts[1]))
	if linesparts[0] == 'predicted label':
	    #print newstatemap[0],int(linesparts[1])

	    if newstatemap[0] != int(linesparts[1]):
		correct = False
	else:
	    allfeatures = linesparts[1].rsplit(' | ')
	    newstatemap.append([[float(x) for x in onefeaturerow.rsplit()] for onefeaturerow in allfeatures])

    if numoflines > 2 and correct:
	states[newstatemap[0]].append([abs(newstatemap[-1][0][0]-newstatemap[-1][1][0])]+newstatemap)


    states[0].sort(reverse=True)
    states[1].sort(reverse=True)

    return states





def read_features_onefile(plot_all_features,use_new_formula,input_feature_names,weights_fn,outputfn):

    all_conv_features = []
    all_subs_features = []

    r.r('''source("/Users/khovsep/Dropbox/eblearn/demos/gene_expression/src/plot_heatmap.R")''')
    drawheatmaps = r.r['draw_heatmaps']

    layerind = -1
    issubsample = False
    for line in open(weights_fn,'r'):

	if len(line)>6 and line[:6]=='######':

	    if layerind >= 0:
		if use_new_formula:
		    for key in features:
			features[key][0] = change_weights(features[key][0],features[key][1],1,0.00001)

		if plot_all_features and not issubsample:

		    for key in features:

			feature = array([features[key][0][key1] for key1 in features[key][0]])
			num_features = size(feature,0)
			size_kernel = size(feature,1)


			if layerind != 0:
			    plot_bars([features[key][0][key1] for key1 in features[key][0]],features[key][0].keys(),'%s_heatmap_feature%d_layer%d_bars.pdf' % (outputfn,key,layerind),abs(feature).max(),False)
			    drawheatmaps(False,False,r.r.matrix(r.r.c(*reshape(feature,size_kernel*num_features,1)), nrow=num_features),4000,r.r.c(*features[key][0].keys()),'%s_heatmap_feature%d_layer%d.pdf' % (outputfn,key,layerind))
			else:
			    plot_bars([features[key][0][key1] for key1 in features[key][0]],[input_feature_names[i] for i in features[key][0].keys()],'%s_heatmap_feature%d_layer%d_bars.pdf' % (outputfn,key,layerind),abs(feature).max(),False)
			    drawheatmaps(False,False,r.r.matrix(r.r.c(*reshape(feature,size_kernel*num_features,1)), nrow=num_features),4000,r.r.c(*[input_feature_names[i] for i in features[key][0].keys()]),'%s_heatmap_feature%d_layer%d.pdf' % (outputfn,key,layerind))


		if issubsample:
		    all_subs_features.append(features)
		else:
		    all_conv_features.append(features)


	    layerind = line.rstrip('######\n').rsplit()[-1]

#	    print 'Layer %s' % (layerind)
	    if layerind[-1] == 'S':
		issubsample=True
		layerind = int(layerind[:-1])
	    else:
		issubsample = False
		layerind = int(layerind)


	    num_features = 0
	    features = {}
	    firstfeat = True

	    input_features_names_used = []
	    continue

	if issubsample:
	    posvals = line.rstrip().rsplit(': ')

	    features[posvals[0]] = [{0:[float(posvals[1].rsplit()[0])]},float(posvals[1].rsplit()[1])]
	    continue



	if line[0] == '=':
	    lastpart = line.rstrip('=========\n').rsplit()[-1]
	    if lastpart != 'BIASES':
		featureind = int(lastpart)
		biases = False
#					print 'NEW FEAT %d' % (int(weights_line.rstrip('=========\n').rsplit()[-1]))
		if not features.has_key(featureind):
		    features[featureind] = [{},0]
#		 print features.keys()
		num_features += 1
		firstfeat = False

	    else:
		biases = True

	    continue



	if firstfeat:
	    continue


	if biases:
	    posvals = line.rstrip().rsplit(': ')
	    features[int(posvals[0])][1] = float(posvals[1])
	else:
	    posvals = line.rstrip().rsplit(' ')
	    if len(posvals[0])>0:
		size_kernel = len(posvals)-1
		features[featureind][0][int(posvals[0].rstrip(':'))] = [float(x) for x in posvals[1:]]
	    else:
		features[featureind][0][0] = [float(x) for x in posvals[1:]]
		size_kernel = len(posvals)-1


    if use_new_formula:
	for key in features:
	    features[key][0] = change_weights(features[key][0],features[key][1],1,0.00001)

    if plot_all_features and not issubsample:

	for key in features:

	    feature = array([features[key][0][key1] for key1 in features[key][0]])
	    num_features = size(feature,0)
	    size_kernel = size(feature,1)


	    if layerind != 0:
		plot_bars([features[key][0][key1] for key1 in features[key][0]],features[key][0].keys(),'%s_heatmap_feature%d_layer%d_bars.pdf' % (outputfn,key,layerind),abs(feature).max(),False)

		drawheatmaps(False,False,r.r.matrix(r.r.c(*reshape(feature,size_kernel*num_features,1)), nrow=num_features),4000,r.r.c(*features[key][0].keys()),'%s_heatmap_feature%d_layer%d.pdf' % (outputfn,key,layerind))
	    else:
		plot_bars([features[key][0][key1] for key1 in features[key][0]],[input_feature_names[i] for i in features[key][0].keys()],'%s_heatmap_feature%d_layer%d_bars.pdf' % (outputfn,key,layerind),abs(feature).max(),False)

		drawheatmaps(False,False,r.r.matrix(r.r.c(*reshape(feature,size_kernel*num_features,1)), nrow=num_features),4000,r.r.c(*[input_feature_names[i] for i in features[key][0].keys()]),'%s_heatmap_feature%d_layer%d.pdf' % (outputfn,key,layerind))


    if issubsample:
	all_subs_features.append(features)
    else:
	all_conv_features.append(features)






    return all_conv_features,all_subs_features



def create_pallet():
    rgb = [[],[],[]]
    for i in xrange(256*256*256):
	x = '%0.6x' % i
	rgb[0].append(int(x[:2],16)*1.0/255)
	rgb[1].append(int(x[2:4],16)*1.0/255)
	rgb[2].append(int(x[4:],16)*1.0/255)

    return rgb



def create_pallet_strings():
    rgb = []
    for i in xrange(256*256*256):
	rgb.append('#%0.6x' % i)

    return rgb



def convert_rgb_to_palletind(c,m,y):
    r = 1 - c
    g = 1 - m
    b = 1 - y

    val = int('0x' + '%0.2x' % int(ceil(r*255)) + '%0.2x' % int(ceil(g*255)) + '%0.2x' % int(ceil(b*255)),16)
    return val



def calculate_merge_cost(cluster1,cluster2,data1,data2,costtype='negcovar'):
    if costtype == 'negcovar':
	return -abs(corrcoef(array([data1[cluster1[0]:cluster2[1]+1],data2[cluster1[0]:cluster2[1]+1]]))[0,1]/(cluster2[1]-cluster1[0]))





def calculate_similarity(cluster1,cluster2,overlap_prop=True,linkage='single'):

    overlaps = []
    for segment1 in cluster1[0]:
	for segment2 in cluster2[0]:
	    overlaps.append(calculate_overlap(segment1,segment2,overlap_prop))

    if linkage == 'single':
	return min(overlaps)
    elif linkage == 'complete':
	return max(overlaps)
    elif linkage == 'upgma':
	return sum(overlaps)*1.0/(len(cluster1[0])*cluster2[0])




def calculate_overlap(segment1,segment2,overlap_prop=True):


    if segment1[1] <= segment2[1] and segment1[2] >= segment2[2]:
	if overlap_prop:
	    overlap = (segment2[2]-segment2[1]+1)*1.0/(segment1[2]-segment1[1]+1)
	else:
	    overlap = segment2[2]-segment2[1]+1

	return overlap
    elif segment1[1] >= segment2[1] and segment1[2] <= segment2[2]:
	if overlap_prop:
	    overlap = (segment1[2]-segment1[1]+1)*1.0/(segment2[2]-segment2[1]+1)
	else:
	    overlap = segment1[2]-segment1[1]+1

	return overlap
    else:
	gapa = segment1[2] - segment2[1]
	gapb = segment2[2] - segment1[1]

	gap = min(gapa,gapb)
	if gap >= 0:
	    if overlap_prop:
		overlap = (gap+1)*1.0/(max(segment2[2]-segment2[1],segment1[2]-segment1[1])+1)
	    else:
		overlap = gap+1
	else:
	    return 0





def combine_clusters(cluster1,cluster2,mergecost):
    new_cluster = [cluster1[0]+cluster2[0],list(unique(cluster1[1]+cluster2[1])),min(cluster1[2],cluster2[2]),max(cluster1[3],cluster2[3]),mergecost]
    return new_cluster



def separate_overlapping_clusters(clusters):

    separate_clust_groups = [[clusters[0]]]
    curclust = 0

#   trackstokeep = [0,1,3]

    for i in xrange(1,len(clusters)):

#	if len(unique(clusters[i][-1]+trackstokeep))-len(clusters[i][-1]) == len(trackstokeep):
#	    continue

	add_new_group = True

	for clust_group in separate_clust_groups:
	    if clusters[i][0] < clust_group[-1][1]:
		continue
	    else:
		add_new_group = False
		clust_group.append(clusters[i])
		break
	if add_new_group:
	    separate_clust_groups.append([clusters[i]])


    return separate_clust_groups






def find_overlapping_clusters(segments,numoftracks,interactiongroups,min_overlap=0,maxnumclusters=12):

    clusters = {}

    for i in xrange(len(segments)):
	clusters[i] = [[segments[i]],segments[i][0],segments[i][1],segments[i][2],1.0]


    merge_cost = {}
    for i in clusters:
	for j in clusters:
	    if len(unique(clusters[i][1] + clusters[j][1])) == 3 and not merge_cost.has_key((i,j)) and not merge_cost.has_key((j,i)):
		merge_cost[i,j] = calculate_similarity(clusters[i],clusters[j])


    newclustind = len(clusters)
#   print clusters.values()

    while len(merge_cost) > 0 and max(merge_cost.values()) > min_overlap and len(clusters)>maxnumclusters:
	clust1,clust2 = max(merge_cost.iterkeys(),key=lambda k:merge_cost[k])

	clusters[newclustind] = combine_clusters(clusters[clust1],clusters[clust2],merge_cost[clust1,clust2])
	#print clusters[clust1],clusters[clust2],clusters[newclustind],max(merge_cost.values())

	del(clusters[clust1])
	del(clusters[clust2])

	del(merge_cost[clust1,clust2])

	for clustkey in merge_cost.keys():
	    if clustkey[0] == clust1 or clustkey[0] == clust2:
		del(merge_cost[clustkey])
		merge_cost[clustkey[1],newclustind] = calculate_similarity(clusters[clustkey[1]],clusters[newclustind])
	    if clustkey[1] == clust1 or clustkey[1] == clust2:
		del(merge_cost[clustkey])
		merge_cost[clustkey[0],newclustind] = calculate_similarity(clusters[clustkey[0]],clusters[newclustind])

	newclustind += 1


    clusterswithtrackcombos = []
    for cluster in clusters.values():
	clusterswithtrackcombos.append([cluster[2],cluster[3],cluster[4],cluster[1]])



    clusterswithtrackcombos.sort()


#   separate_inds_for = [[0,1],[0,3],[9,10],[0,5],[0,7]]
    unoverlapping_clusters = []
    for ind in interactiongroups:
	newclustercombos = [cluster for cluster in clusterswithtrackcombos if len(unique(ind + cluster[-1]))-len(cluster[-1])==0]#len(ind)]
	print [cluster[-1] for cluster in newclustercombos]
	unoverlapping_clusters.append(separate_overlapping_clusters(newclustercombos))



#   unoverlapping_clusters = separate_overlapping_clusters(clusterswithtrackcombos)

#   print unoverlapping_clusters

    clusters_by_trackpairs = {}
    for cluster in clusterswithtrackcombos:
	tracks = sorted(cluster[-1])
	for i in xrange(len(tracks)-1):
	    for j in xrange(1,len(tracks)):
		if not clusters_by_trackpairs.has_key((tracks[i],tracks[j])):
		    clusters_by_trackpairs[tracks[i],tracks[j]] = [[cluster[0],cluster[1],cluster[2]]]
		else:
		    clusters_by_trackpairs[tracks[i],tracks[j]].append([cluster[0],cluster[1],cluster[2]])



    return unoverlapping_clusters,clusters_by_trackpairs





def bottomup(data1,data2,max_error=0,maxnumclusters=12):
    clusters = []

    for i in xrange(len(data1)):
	clusters.append([i,i])#,0])

    merge_cost = []
    for i in xrange(len(clusters)-1):
	merge_cost.append(calculate_merge_cost(clusters[i],clusters[i+1],data1,data2))

    while len(merge_cost) > 0 and min(merge_cost) < max_error or len(clusters)>maxnumclusters:

	i = argmin(merge_cost)
	clusters[i] = [clusters[i][0],clusters[i+1][1]]#,-calculate_merge_cost(clusters[i],clusters[i+1],data1,data2)]
	clusters.pop(i+1)
	merge_cost.pop(i)
	if i < len(merge_cost):
	    merge_cost[i] = calculate_merge_cost(clusters[i],clusters[i+1],data1,data2)
	if i > 0:
	    merge_cost[i-1] = calculate_merge_cost(clusters[i-1],clusters[i],data1,data2)


#   print [x[1]-x[0] for x in clusters]
    return clusters





def call_ggobi(data,rownames,groupstarts):
    r.r('''source("/Users/khovsep/Dropbox/eblearn/demos/gene_expression/src/plot_heatmap.R")''')

    plot_with_ggobi = r.r['plot_with_ggobi']

    N4 = len(data[0])
    colnames = []
    classes = []
    for i in xrange(len(data)):
	newdata = []
	if i == 0:
	    colnames = ['%d' % (posmarker-N4/2) for posmarker in xrange(len(data[0]))]
	    classind = 0
	    for groupind,groupstart in enumerate(groupstarts):
		groupend = groupstarts[groupind+1] if groupind != len(groupstarts)-1 else N4

		classes.extend([groupind+1 for j in xrange(groupstart,groupend)])




    plot_with_ggobi(r.r.matrix(r.r.c(*reshape(array(data),len(data)*len(colnames),1)), nrow=len(data)),r.r.c(*colnames),r.r.c(*rownames),r.r.c(*classes))
    raw_input("Press Enter to Stop ggobi...")




#def compute_profile_proximity_plots(alldistinctpeaks,states,N,pixelsize,outputfn,whichclass,input_feature_names,xaxis_ticks,xaxis_labels,allinnerfeaturekeys):
#   TFpeaks = alldistinctpeaks[0]
#   window = 100
#   for peak in TFpeaks:
#	peak_region = [peak[1]-window/(pixelsize*1.0),peak[1]+window/(pixelsize*1.0)]


def store_peaks(alldistinctpeaks,N,pixelsize,outputfn,whichclass,input_feature_names,xaxis_ticks,xaxis_labels,allinnerfeaturekeys):

    outputfile = open('%s_peaks_class%s.dat' % (outputfn,whichclass),'w')
    for ind in alldistinctpeaks:
	print >> outputfile, '=== %s Peaks ===' % input_feature_names[allinnerfeaturekeys[ind]]
	for peak in alldistinctpeaks[ind]:
	    print >> outputfile, '%d %f' % ((peak[1]-N/2)*pixelsize,peak[2])





def findpeaks(feature,keptpeakmovepercent,keptpeakheightpercent,outputfn,whichclass,factornames,xaxis_ticks,xaxis_labels,allinnerfeaturekeys,whichfeaturespeaks):
    n = len(feature)

    alldistinctpeaks = {}
#   pp = PdfPages('%s_histograms_of_peakmovements_class%s.pdf' % (outputfn,str(whichclass)))

#   f1 = figure(1)
    for featureind,onefeature in enumerate(feature):
	if featureind not in whichfeaturespeaks:
	    continue

	N = len(onefeature)

	absfeature = abs(array(onefeature))
	diffs = absfeature[1:]-absfeature[:-1]

	peaklocs = [i+1 for i in xrange(N-1) if i < N-2 and diffs[i]>=0 and diffs[i+1] < 0]
#	valleylocs = [0]+[i+1 for i in xrange(N-1) if i < N-2 and diffs[i]<=0 and diffs[i+1] > 0]+[N-1]




	peakvals = [[absfeature[i],i,onefeature[i]] for i in peaklocs]

	allmoves = []
	currentpeak = peakvals[0]
	for peakind,peak in enumerate(peakvals[1:]):
	    minval = min(absfeature[currentpeak[1]:peak[1]+1])
	    moveleft = abs(currentpeak[0]-minval)
	    moveright = abs(peak[0]-minval)
	    allmoves.extend([moveleft,moveright])
	    currentpeak = peak



	allmoves.sort()
	movelimit = allmoves[int(len(allmoves)*keptpeakmovepercent)]

	peakvals = sorted(peakvals,reverse=True)

	chosenpeaks = peakvals[:int(len(peakvals)*keptpeakheightpercent)]
	chosenpeaks.sort(key=lambda chosenpeaks: chosenpeaks[1])


	groupsofpeaks = [[chosenpeaks[0]]]

	for peakind,peak in enumerate(chosenpeaks[1:]):
	    minval = min(absfeature[groupsofpeaks[-1][0][1]:peak[1]+1])
	    moveleft = abs(groupsofpeaks[-1][0][0]-minval)
	    moveright = abs(peak[0]-minval)
	    if moveleft < movelimit or moveright < movelimit:
		groupsofpeaks[-1].append(peak)
	    else:
		groupsofpeaks.append([peak])


	distinct_peaks = [peakgroup[i] for peakgroup in groupsofpeaks for i in xrange(len(peakgroup)) if peakgroup[i][0]==max([peak[0] for peak in peakgroup])]
	alldistinctpeaks[featureind] = distinct_peaks

#	f2 = figure(2)
#	hist(allmoves,bins=10,normed=True)
#	title('Sizes of peaks for %s' % (factornames[allinnerfeaturekeys[featureind]]))
#	pp.savefig(f2)
#	close()


#	ax = f1.add_subplot(n,1,featureind+1)
#	ax.plot(range(len(onefeature)),[0]+list(diffs),label='diffs',color='b',linewidth=3)

#	ymin,ymax = ylim()

#	for peak in distinct_peaks:
#	    ax.plot([peak[1],peak[1]],[ymin,ymax],'r')

#	ylabel(factornames[allinnerfeaturekeys[featureind]].replace('_',' '),fontsize=40,rotation='horizontal')


#	ymin, ymax = ylim()
#	plot([N/2,N/2],[ymin,ymax],'k',linewidth=2)
#	plot([0,N-1],[0,0],'k',linewidth=1)

#	xticks(xaxis_ticks,xaxis_labels)
#
#	locs,labels = yticks()
#	yticks(locs, map(lambda x: "%g" % x, locs))


#	xlim(0,N-1)

#	ax.yaxis.grid(True)
#	ax.xaxis.grid(True)

##	for tick in ax.get_yticklabels():
##	    tick.set_color('b')
##	    tick.set_fontsize(10)
##
##
##	for tick in ax.get_xticklabels():
##	    tick.set_fontsize(30)
##
##
##  f1.set_figwidth(40)
##  f1.set_figheight(60)
##
##  pp.savefig(f1)
##  close()
##  pp.close()

    return alldistinctpeaks




def findpeaks_nofilters(feature,outputfn,whichclass,factornames,xaxis_ticks,xaxis_labels,allinnerfeaturekeys,whichfeaturespeaks):
    n = len(feature)

    alldistinctpeaks = {}
#   pp = PdfPages('%s_histograms_of_peakmovements_class%s.pdf' % (outputfn,str(whichclass)))

#   f1 = figure(1)
    for featureind,onefeature in enumerate(feature):
	if featureind not in whichfeaturespeaks:
	    continue

	N = len(onefeature)

	absfeature = abs(array(onefeature))
	diffs = absfeature[1:]-absfeature[:-1]

	peaklocs = [i+1 for i in xrange(N-1) if i < N-2 and diffs[i]>=0 and diffs[i+1] < 0]
#		valleylocs = [0]+[i+1 for i in xrange(N-1) if i < N-2 and diffs[i]<=0 and diffs[i+1] > 0]+[N-1]




	peakvals = [[[absfeature[i],i,onefeature[i]] for i in peaklocs if onefeature[i]>0],[[absfeature[i],i,onefeature[i]] for i in peaklocs if onefeature[i]<0]]

#	allmoves = []
#	currentpeak = peakvals[0]
#	for peakind,peak in enumerate(peakvals[1:]):
#	    minval = min(absfeature[currentpeak[1]:peak[1]+1])
#	    moveleft = abs(currentpeak[0]-minval)
#	    moveright = abs(peak[0]-minval)
#	    allmoves.extend([moveleft,moveright])
#	    currentpeak = peak



#	allmoves.sort()
#	movelimit = allmoves[int(len(allmoves)*keptpeakmovepercent)]

#	peakvals = sorted(peakvals,reverse=True)

#	chosenpeaks = peakvals#peakvals[:int(len(peakvals)*keptpeakheightpercent)]
#	chosenpeaks.sort(key=lambda chosenpeaks: chosenpeaks[1])

#	groupsofpeaks = [[chosenpeaks[0]]]

#	for peakind,peak in enumerate(chosenpeaks[1:]):
#	    minval = min(absfeature[groupsofpeaks[-1][0][1]:peak[1]+1])
#	    moveleft = abs(groupsofpeaks[-1][0][0]-minval)
#	    moveright = abs(peak[0]-minval)
#	    if moveleft < movelimit or moveright < movelimit:
#		groupsofpeaks[-1].append(peak)
#	    else:
#		groupsofpeaks.append([peak])


#	distinct_peaks = chosenpeaks#[peakgroup[i] for peakgroup in groupsofpeaks for i in xrange(len(peakgroup)) if peakgroup[i][0]==max([peak[0] for peak in peakgroup])]
	alldistinctpeaks[featureind] = peakvals#chosenpeaks#distinct_peaks

	if False:
	    f2 = figure(2)
	    hist(allmoves,bins=10,normed=True)
	    title('Sizes of peaks for %s' % (factornames[allinnerfeaturekeys[featureind]]))
	    pp.savefig(f2)
	    close()


	    ax = f1.add_subplot(n,1,featureind+1)
	    ax.plot(range(len(onefeature)),[0]+list(diffs),label='diffs',color='b',linewidth=3)

	    ymin,ymax = ylim()

	    for peak in distinct_peaks:
		ax.plot([peak[1],peak[1]],[ymin,ymax],'r')

	    ylabel(factornames[allinnerfeaturekeys[featureind]].replace('_',' '),fontsize=40,rotation='horizontal')


	    ymin, ymax = ylim()
	    plot([N/2,N/2],[ymin,ymax],'k',linewidth=2)
	    plot([0,N-1],[0,0],'k',linewidth=1)

	    xticks(xaxis_ticks,xaxis_labels)

	    locs,labels = yticks()
	    yticks(locs, map(lambda x: "%g" % x, locs))


	    xlim(0,N-1)

	    ax.yaxis.grid(True)
	    ax.xaxis.grid(True)

	    for tick in ax.get_yticklabels():
		tick.set_color('b')
		tick.set_fontsize(10)


	    for tick in ax.get_xticklabels():
		tick.set_fontsize(30)


#   f1.set_figwidth(40)
#   f1.set_figheight(60)

#   pp.savefig(f1)
#   close()
#   pp.close()

    return alldistinctpeaks








##def compute_profile_relationship_plots(feature,input):
##  Hists = [1,3,5]
##
##  TFprofiles = sorted(feature[0])
##  N = len(TFprofiles)
##
##
##  diffs = TFprofiles[1:]-TFprofiles[:-1]
##  bands
##
##  nonzeroleftbounds = [i for i in xrange(N) if TFprofiles[i] != 0 and (TFprofiles[i-1] == 0 or i == 0)]
##  nonzerorightbounds = [i for i in xrange(N) if TFprofiles[i] != 0 and (TFprofiles[i+1] == 0 or i == N-1)]
##
##  nonzeroperiods = [[nonzeroleftbounds[i],nonzerorightbounds[i]] for i in xrange(len(nonzeroleftbounds))]
##
##  nonzeropeaks = [max(abs(TFprofiles[nonzeroperiods[i][0]:nonzeroperiods[i][1]+1])) for i in xrange(len(nonzeroperiods))]
##
##
##
##
##
##
##
##  TFprofiles
##
##
##
##  numbinsTFprofs = 10.0
##  TFprofbins = arange(TFprofiles[0],TFprofiles[-1],(TFprofiles[-1]-TFprofiles[0])/numbinsTFprofs)
##  numsubplotcols = 4.0 if len(TFprofbins) / 4.0 >= 4 else 3.0 if len(TFprofbins) /3.0 >= 3 else 2.0
##  numsubplotrows = ceil(len(TFprofbins)/numsubplotcols)
##
##  f1 = figure(1)
##  TFprofiles = feature[0]
##
##
##  for TFindbin,TFprofileleft in enumerate(TFprofbins):
##	TFprofileright = TFprofbins[TFindbin+1] if TFindbin < len(TFprofbins)-1 else sorted(TFprofiles)[-1]
##	print TFprofileleft,TFprofileright
##	relevant_prof_inds = [i for i in xrange(N4) if TFprofiles[i]>=TFprofileleft and TFprofiles[i]<TFprofileright]
##	print [(x-N4/2)*pixelsize for x in relevant_prof_inds]
##	ax = f1.add_subplot(numsubplotrows,numsubplotcols,TFindbin+1)
##
##	for Histind,Hist in enumerate(Hists):
##	    Histprofs = feature[Hist]
##	    TFpeaktoHistpeaks = {}
##	    for TFind,TFpos in enumerate(relevant_prof_inds):
##		for Histpos in xrange(N4):
##		    if not TFpeaktoHistpeaks.has_key(((Histpos-TFpos)*pixelsize)):
##			TFpeaktoHistpeaks[((Histpos-TFpos)*pixelsize)] = []
##		    TFpeaktoHistpeaks[((Histpos-TFpos)*pixelsize)].append(Histprofs[Histpos])
##
##	    x = sorted(TFpeaktoHistpeaks.keys())
##	    y = [mean(TFpeaktoHistpeaks[key]) for key in x]
##	    ax.plot(x,y,lines[Hist],label=input_feature_names[Hist])
##
##	legend()
##	ax.set_xlabel('Relative distance from %s profile' % input_feature_names[0])
##	ax.set_ylabel('Average Profile')
##	title('Average Histone profiles relative to %s profiles in [%g,%g] range' % (input_feature_names[0],TFprofileleft,TFprofileright))
##
##
##  f1.set_figwidth(40)
##  f1.set_figheight(60)
##
##
##  savefig('%s_TFpeakstoHistpeaks_relationship_class%s.pdf' % (outputfn,str(whichclass)))
##  close()
##
##
##
##



def create_combined_template(lowerfeatures,upperfeatures,subsample):

    N2 = len(upperfeatures[0])	# number of units in 3rd layer from top (number of pixels in the heatmaps)
    top = range(len(upperfeatures))

    entirematrix_wholeouterfeature_curves = {}
    firstfeature_coarse = {}


    for i in top:
	outerfeature = lowerfeatures[i][0]

	numouterfeatures = len(outerfeature)
	for key in outerfeature:
	    N1 = len(outerfeature[key])	  # size of the filter

	N4 = (N2*subsample-1)+N1

	ind = 0

	multweights = upperfeatures[i]

	for n in xrange(N2):
	    themultweight = multweights[n]
	    for s in xrange(subsample):
		for j in xrange(N1):
		    for key1 in outerfeature:
			if entirematrix_wholeouterfeature_curves.has_key(key1):
			    entirematrix_wholeouterfeature_curves[key1][ind+j].append(outerfeature[key1][j]*themultweight)
			else:
			    entirematrix_wholeouterfeature_curves[key1] = [[] for k in xrange(N4)]
			    entirematrix_wholeouterfeature_curves[key1][ind+j].append(outerfeature[key1][j]*themultweight)
		ind += 1



    allinnerfeaturekeys = []
    for key in entirematrix_wholeouterfeature_curves:
	if key not in allinnerfeaturekeys:
	    allinnerfeaturekeys.append(key)

    allinnerfeaturekeys.sort()

    return entirematrix_wholeouterfeature_curves,allinnerfeaturekeys,N4





#def create_feature_curves_overlapping_withtoplayer_usediffs(topfeatures_diff,whichclass,allfeatures_copy,topxpercent,groupsize,outputgroupsize,interactiongroups,binsize,subsample1,pixelsize,outputfn,input_feature_names,sortingmetric,do_average_over_subfeatures,states=[]):
def create_feature_curves_overlapping_withtoplayer_usediffs(topfeatures_diff,whichclass,allfeatures_copy,groupsize,outputgroupsize,interactiongroups,subsample1,pixelsize,keptpeakmovepercent,keptpeakheightpercent,whichfeaturespeaks,outputfn,input_feature_names,alldistinctpeaks,states=[]):


    r.r('''source("/Users/khovsep/Dropbox/eblearn/demos/gene_expression/src/plot_heatmap.R")''')
    drawheatmaps = r.r['draw_heatmaps']




    cols = ['b','g','r','k','y','c','m']
    marks = ['o','s','v','d','*','+','h']


#   lines = ['-','--',':']
    marks = [[mark,col] for col in cols for mark in marks]
#   lines = [col+line for line in lines for col in cols]

    marksnew = ['o','s','v','d','*','+','h']


    x = range(20)
    shuffle(x)
    cols = cm.Paired(x)
    marks = marks + [[mark,col] for col in cols for mark in marksnew]

#   print marks
    lines = ['k-','b-','b--','r-','r--','g-','g--','m-','m--','c-','c--']


    numtracks = len(input_feature_names)

    allfeatures = allfeatures_copy

    entirematrix_wholeouterfeature_curves, allinnerfeaturekeys,N4 = create_combined_template(allfeatures[2],topfeatures_diff,subsample1)

##  N2 = len(topfeatures_diff[0])  # number of units in 3rd layer from top (number of pixels in the heatmaps)
##
##  top = range(len(topfeatures_diff))
##
##
##  features_coarse = allfeatures[2]
##
##  entirematrix_wholeouterfeature_curves = {}
##
##  firstfeature_coarse = {}
##
##  for i in top:
##
###	    print '=================', i, '===================='
##	outerfeature = features_coarse[i][0]
##
##	numouterfeatures = len(outerfeature)
##	for key in outerfeature:
##	    N1 = len(outerfeature[key])	  # size of the filter
##
##	N4 = (N2*subsample1-1)+N1
##
##	ind = 0
##
##	multweights = topfeatures_diff[i]
##
##	for n in xrange(N2):
##	    themultweight = multweights[n]
##	    for s in xrange(subsample1):
##		for j in xrange(N1):
##		    for key1 in outerfeature:
##			if entirematrix_wholeouterfeature_curves.has_key(key1):
##			    entirematrix_wholeouterfeature_curves[key1][ind+j].append(outerfeature[key1][j]*themultweight)
##			else:
##			    entirematrix_wholeouterfeature_curves[key1] = [[] for k in xrange(N4)]
##			    entirematrix_wholeouterfeature_curves[key1][ind+j].append(outerfeature[key1][j]*themultweight)
##		ind += 1
##
##
##
##  allinnerfeaturekeys = []
##  for key in entirematrix_wholeouterfeature_curves:
##	if key not in allinnerfeaturekeys:
##	    allinnerfeaturekeys.append(key)
##
##  allinnerfeaturekeys.sort()

#   N4 = (N2*subsample1-1)+N1


    feature = [[sum(entirematrix_wholeouterfeature_curves[key][ind]) for ind in xrange(N4)] for key in allinnerfeaturekeys]
    maxval = abs(array(feature)).max()

    groupstarts = list(sorted(-array([x for x in arange(0,N4/2+1,1) if x * pixelsize % groupsize == 0])+N4/2)) +list(array([x for x in arange(1,N4/2,1) if x*pixelsize % groupsize == 0])+N4/2)

    if groupstarts[0] <> 0:
	groupstarts[0] = 0


    xaxis_labels_heatmap = ['%s' % (str(x*pixelsize/1000.0) + 'Kb' if abs(x*pixelsize) >= 1000 and abs(x*pixelsize) < 1000000 else (str(x*pixelsize/1000000.0)+'Mb' if abs(x*pixelsize) >= 1000000 else str(x*pixelsize)+'b')) for x in [(groupstart-N4/2) for groupstart in groupstarts]]



    features_sumvals = zeros([len(feature),len(groupstarts)])

    for groupind,groupstart in enumerate(groupstarts):
	groupend = groupstarts[groupind+1] if groupind != len(groupstarts)-1 else N4
	for ind,rowvals in enumerate(feature):
	    features_sumvals[ind,groupind] = sum(rowvals[groupstart:groupend])

    drawheatmaps(False,False,r.r.matrix(r.r.c(*reshape(features_sumvals,len(feature)*len(groupstarts),1)), nrow=len(feature)),4000,r.r.c(*input_feature_names),'%s_alllayers_sumabsfeatures_bygroups_class%s.pdf' % (outputfn,str(whichclass)),r.r.c(*xaxis_labels_heatmap))



    features_maxvals = zeros([len(feature),len(groupstarts)])

    for groupind,groupstart in enumerate(groupstarts):
	groupend = groupstarts[groupind+1] if groupind != len(groupstarts)-1 else N4
	for ind,rowvals in enumerate(feature):
	    t = max(rowvals[groupstart:groupend])
	    t1 = min(rowvals[groupstart:groupend])
	    if abs(t) > abs(t1):
		features_maxvals[ind,groupind] = t
	    else:
		features_maxvals[ind,groupind] = t1

    drawheatmaps(False,False,r.r.matrix(r.r.c(*reshape(features_maxvals,len(feature)*len(groupstarts),1)), nrow=len(feature)),4000,r.r.c(*input_feature_names),'%s_alllayers_maxabsfeatures_bygroups_class%s.pdf' % (outputfn,str(whichclass)),r.r.c(*xaxis_labels_heatmap))




    f1 = figure(1)

    ax = f1.add_subplot(111)
    intermat,vmin,vmax,colormap = mycorrcoef(array(feature),'cov')
    ax.imshow(intermat,interpolation='nearest',cmap=get_cmap(colormap),vmin=vmin,vmax=vmax)
    xticks(arange(numtracks),arange(numtracks))
    yticks(arange(numtracks),input_feature_names)
    for tick in ax.get_yticklabels():
	tick.set_fontsize(28)

    f1.set_figwidth(22)
    f1.set_figheight(20)

    savefig('%s_curves_alllayers_covs_corrs_combinedtracks_abs_class%s.pdf' % (outputfn,str(whichclass)))
    close()



    corr_mats_combinedtracks_diff_abs_bygroups = {}
    cov_mats_combinedtracks_diff_abs_bygroups = {}
    dot_mats_combinedtracks_diff_abs_bygroups = {}
    angle_mats_combinedtracks_diff_abs_bygroups = {}


    numsubplotcols = 4.0 if len(groupstarts) / 4.0 >= 4 else 3.0 if len(groupstarts) /3.0 >= 3 else 2.0

    pp = PdfPages('%s_curves_alllayers_covs_corrs_combinedtracks_class%s_abs_bygroups.pdf' % (outputfn,str(whichclass)))

    f1 = figure(3)
    for groupind,groupstart in enumerate(groupstarts):
	groupend = groupstarts[groupind+1] if groupind != len(groupstarts)-1 else N4

	ax = f1.add_subplot(ceil(len(groupstarts)/numsubplotcols),numsubplotcols,groupind+1)
	cov_mats_combinedtracks_diff_abs_bygroups[groupstart,groupend],vmin,vmax,colormap = mycorrcoef(array([rowvals[groupstart:groupend] for rowvals in feature]),'cov')
	ax.imshow(cov_mats_combinedtracks_diff_abs_bygroups[groupstart,groupend],interpolation='nearest',cmap=get_cmap(colormap),vmin=vmin,vmax=vmax)
	xticks(arange(numtracks),arange(numtracks))
	yticks(arange(numtracks),input_feature_names)
	title('%s -- %s' % (str((groupstart-N4/2)*pixelsize/1000.0) + 'Kb' if abs((groupstart-N4/2)*pixelsize) >= 1000 and abs((groupstart-N4/2)*pixelsize) < 1000000 else (str((groupstart-N4/2)*pixelsize/1000000.0)+'Mb' if abs((groupstart-N4/2)*pixelsize) >= 1000000 else str((groupstart-N4/2)*pixelsize)+'b'),str((groupend-N4/2)*pixelsize/1000.0) + 'Kb' if abs((groupend-N4/2)*pixelsize) >=1000 and abs((groupend-N4/2)*pixelsize) < 1000000 else (str((groupend-N4/2)*pixelsize/1000000.0)+'Mb' if abs((groupend-N4/2)*pixelsize) >= 1000000 else str((groupend-N4/2)*pixelsize)+'b')),fontsize=25)
	for tick in ax.get_yticklabels():
	    tick.set_fontsize(25)



    if numsubplotcols == 2:
	f1.set_figwidth(40)
	f1.set_figheight(30)
    elif numsubplotcols == 3:
	f1.set_figwidth(32)
	f1.set_figheight(30)
    elif numsubplotcols == 4:
	f1.set_figwidth(40)
	f1.set_figheight(50)


    pp.savefig()
    close()



    pp.close()


    if True:

	allclusters = []
	if False:
	    f1 = figure(1)
	    f2 = figure(2)
	plotind = 0

	for i in xrange(numtracks-1):
	    plotind += i
	    for j in xrange(i+1,numtracks):
		plotind += 1

		if i==j:
		    if i == 0:
			ax = f1.add_subplot(numtracks-1,numtracks-1,plotind)
			title(input_feature_names[j],fontsize=20)
			ylabel(input_feature_names[i],fontsize=20)


		    continue


		groups = bottomup(array(feature[i]),array(feature[j]),max_error=-0.7,maxnumclusters=30)
		allclusters.extend([[[i,j],x[0],x[1]] for x in groups])


		if False:
		    ax = f1.add_subplot(numtracks-1,numtracks-1,plotind)
		    ax2 = f2.add_subplot(numtracks-1,numtracks-1,plotind)
		    newdata = []


		    for groupind,group in enumerate(groups):
			groupstart = group[0]
			groupend = group[1]

			newfeature = [array(feature[j])[groupstart:groupend+1],array(feature[i])[groupstart:groupend+1]]
			ax.plot(newfeature[0],newfeature[1],marker=marks[groupind][0],c=marks[groupind][1],label='%s--%s' % (str((groupstart-N4/2)*pixelsize/1000.0) + 'Kb' if abs((groupstart-N4/2)*pixelsize) >= 1000 and abs((groupstart-N4/2)*pixelsize) < 1000000 else (str((groupstart-N4/2)*pixelsize/1000000.0)+'Mb' if abs((groupstart-N4/2)*pixelsize) >= 1000000 else str((groupstart-N4/2)*pixelsize)+'b'),str((groupend-N4/2)*pixelsize/1000.0) + 'Kb' if abs((groupend-N4/2)*pixelsize) >=1000 and abs((groupend-N4/2)*pixelsize) < 1000000 else (str((groupend-N4/2)*pixelsize/1000000.0)+'Mb' if abs((groupend-N4/2)*pixelsize) >= 1000000 else str((groupend-N4/2)*pixelsize)+'b')),markersize=3)

			ax2.plot([],[],marker=marks[groupind][0],c=marks[groupind][1],label='%s--%s' % (str((groupstart-N4/2)*pixelsize/1000.0) + 'Kb' if abs((groupstart-N4/2)*pixelsize) >= 1000 and abs((groupstart-N4/2)*pixelsize) < 1000000 else (str((groupstart-N4/2)*pixelsize/1000000.0)+'Mb' if abs((groupstart-N4/2)*pixelsize) >= 1000000 else str((groupstart-N4/2)*pixelsize)+'b'),str((groupend-N4/2)*pixelsize/1000.0) + 'Kb' if abs((groupend-N4/2)*pixelsize) >=1000 and abs((groupend-N4/2)*pixelsize) < 1000000 else (str((groupend-N4/2)*pixelsize/1000000.0)+'Mb' if abs((groupend-N4/2)*pixelsize) >= 1000000 else str((groupend-N4/2)*pixelsize)+'b')))
			legend(markerscale=0.5,numpoints=1,ncol=2)
			leg = gca().get_legend()
			ltext  = leg.get_texts()
			setp(ltext, fontsize='xx-small')

		    if i == 0:
			figure(1)
			title(input_feature_names[j],fontsize=20)
			figure(2)
			title(input_feature_names[j],fontsize=20)


		    if j == i+1:
			figure(1)
			ylabel(input_feature_names[i],fontsize=20)
			figure(2)
			ylabel(input_feature_names[i],fontsize=20)



	if False:
	    f1.set_figwidth(50)
	    f1.set_figheight(50)
	    f1.savefig('%s_scattermatrix_tracks_class%s.pdf' % (outputfn,str(whichclass)))
	    close()

	    f2.set_figwidth(50)
	    f2.set_figheight(50)
	    f2.savefig('%s_scattermatrix_tracks_clusters_class%s.pdf' % (outputfn,str(whichclass)))
	    close()


	combined_clusters_sep,combined_clusters = find_overlapping_clusters(allclusters,numtracks,interactiongroups,min_overlap=0.8)

	if False:
	    f1 = figure(1)
	    f2 = figure(2)
	    plotind = 0

	    for i in xrange(numtracks-1):
		plotind += i
		for j in xrange(i+1,numtracks):
		    plotind += 1

		    ax = f1.add_subplot(numtracks-1,numtracks-1,plotind)
		    ax2 = f2.add_subplot(numtracks-1,numtracks-1,plotind)
		    newdata = []


		    groupind = -1
		    for group in combined_clusters[i,j]:

			groupind += 1

			groupstart = group[0]
			groupend = group[1]

			newfeature = [array(feature[j])[groupstart:groupend+1],array(feature[i])[groupstart:groupend+1]]
			ax.plot(newfeature[0],newfeature[1],marker=marks[groupind][0],c=marks[groupind][1],label='%s--%s' % (str((groupstart-N4/2)*pixelsize/1000.0) + 'Kb' if abs((groupstart-N4/2)*pixelsize) >= 1000 and abs((groupstart-N4/2)*pixelsize) < 1000000 else (str((groupstart-N4/2)*pixelsize/1000000.0)+'Mb' if abs((groupstart-N4/2)*pixelsize) >= 1000000 else str((groupstart-N4/2)*pixelsize)+'b'),str((groupend-N4/2)*pixelsize/1000.0) + 'Kb' if abs((groupend-N4/2)*pixelsize) >=1000 and abs((groupend-N4/2)*pixelsize) < 1000000 else (str((groupend-N4/2)*pixelsize/1000000.0)+'Mb' if abs((groupend-N4/2)*pixelsize) >= 1000000 else str((groupend-N4/2)*pixelsize)+'b')),markersize=3)

			ax2.plot([],[],marker=marks[groupind][0],c=marks[groupind][1],label='%s--%s' % (str((groupstart-N4/2)*pixelsize/1000.0) + 'Kb' if abs((groupstart-N4/2)*pixelsize) >= 1000 and abs((groupstart-N4/2)*pixelsize) < 1000000 else (str((groupstart-N4/2)*pixelsize/1000000.0)+'Mb' if abs((groupstart-N4/2)*pixelsize) >= 1000000 else str((groupstart-N4/2)*pixelsize)+'b'),str((groupend-N4/2)*pixelsize/1000.0) + 'Kb' if abs((groupend-N4/2)*pixelsize) >=1000 and abs((groupend-N4/2)*pixelsize) < 1000000 else (str((groupend-N4/2)*pixelsize/1000000.0)+'Mb' if abs((groupend-N4/2)*pixelsize) >= 1000000 else str((groupend-N4/2)*pixelsize)+'b')))
			legend(markerscale=0.5,numpoints=1,ncol=2)
			leg = gca().get_legend()
			ltext  = leg.get_texts()
			setp(ltext, fontsize='xx-small')

		    if i == 0:
			figure(1)
			title(input_feature_names[j],fontsize=20)
			figure(2)
			title(input_feature_names[j],fontsize=20)


		    if j == i+1:
			figure(1)
			ylabel(input_feature_names[i],fontsize=20)
			figure(2)
			ylabel(input_feature_names[i],fontsize=20)




	    f1.set_figwidth(50)
	    f1.set_figheight(50)
	    f1.savefig('%s_scattermatrix_combined_tracks_class%s.pdf' % (outputfn,str(whichclass)))
	    close()

	    f2.set_figwidth(50)
	    f2.set_figheight(50)
	    f2.savefig('%s_scattermatrix_tracks_combined_clusters_class%s.pdf' % (outputfn,str(whichclass)))
	    close()





    xaxis_ticks = list(sorted(-array([x for x in arange(0,N4/2+1,1) if x*pixelsize % outputgroupsize == 0])+N4/2)) +list(array([x for x in arange(1,N4/2+1,1) if x*pixelsize % outputgroupsize == 0])+N4/2)
    xaxis_labels = ['%s' % ('' if x == 0 else str(x*pixelsize/1000.0) + 'Kb' if abs(x*pixelsize) >= 1000 and abs(x*pixelsize) < 1000000 else (str(x*pixelsize/1000000.0)+'Mb' if abs(x*pixelsize) >= 1000000 else str(x*pixelsize)+'b')) for x in list(array([x for x in arange(0,N4/2+1,1) if x*pixelsize % outputgroupsize == 0])-N4/2) +list(array([x for x in arange(1,N4/2+1,1) if x*pixelsize % outputgroupsize == 0]))]


    print xaxis_ticks
    print xaxis_labels

    if True:

	#sep_inds = [0,1,3]
#	interactiongroups
#	sep_inds = [[0,1],[0,3],[9,10],[0,5],[0,7]]


	for sepind,combined_clust_per_ind in enumerate(combined_clusters_sep):
	    numofrows = len(combined_clust_per_ind)


	    f1 = figure(1)
	    for groupind,sep_group in enumerate(combined_clust_per_ind):

		maxval = -inf
		minval = inf
		ax = f1.add_subplot(numofrows+1,1,groupind+1)
		for cluster in sep_group:
		    for trackind,track in enumerate(cluster[-1]):
			onefeature = feature[track][cluster[0]:cluster[1]+1]
			maxval = max(maxval,max(onefeature))
			minval = min(minval,min(onefeature))

		for cluster in sep_group:
		    for trackind,track in enumerate(cluster[-1]):
			onefeature = array(feature[track][cluster[0]:cluster[1]+1])

			maxlocal = max(onefeature)
			minlocal = min(onefeature)

			onefeature = (onefeature-minlocal)/(maxlocal-minlocal)*(maxval-minval)+minval

			ax.plot(range(cluster[0],cluster[1]+1),list(onefeature),lines[track])

		ymin, ymax = ylim()
		plot([N4/2,N4/2],[ymin,ymax],'k',linewidth=1)
		xticks(xaxis_ticks,xaxis_labels)
		xlim(0,N4-1)

		ax.yaxis.grid(True)
		ax.xaxis.grid(True)

		for tick in ax.get_xticklabels():
		    tick.set_fontsize(15)


	    ax = f1.add_subplot(len(combined_clust_per_ind)+1,1,len(combined_clust_per_ind)+1)
	    for trackind,tracklabel in enumerate(input_feature_names):
		ax.plot([],[],lines[trackind],label=tracklabel.replace('_',' '),linewidth=1)

	    legend()

	    f1.set_figwidth(80)
	    f1.set_figheight((numofrows+1)*4)


	    savefig('%s_curves_with_interactions_for%s_class%s.pdf' % (outputfn,'+'.join([input_feature_names[x] for x in interactiongroups[sepind]]),str(whichclass)))
	    close()





#   alldistinctpeaks1 = findpeaks(feature,keptpeakmovepercent,keptpeakheightpercent,outputfn,whichclass,input_feature_names,xaxis_ticks,xaxis_labels,allinnerfeaturekeys, whichfeaturespeaks)

    #compute_profile_proximity_plots(alldistinctpeaks,states,outputfn,whichclass,input_feature_names,xaxis_ticks,xaxis_labels,allinnerfeaturekeys)
#   store_peaks(alldistinctpeaks,N4,pixelsize,outputfn,whichclass,input_feature_names,xaxis_ticks,xaxis_labels,allinnerfeaturekeys)

    maxval = abs(array(feature)).max()
#   if states != []:
#	for exampleind,onestate in enumerate(states[0]+states[1]):
#	    maxvalstates = abs(array(onestate[2])).max()

    f1 = figure(1)
    for ind in xrange(len(feature)):


	onefeature = feature[ind]


	ax = f1.add_subplot(len(feature),1,ind+1)
#		ax.fill_between(range(len(onefeature)),0,onestate[2][ind],alpha=1,label='data',color='y')
#		ylim(-maxvalstates,maxvalstates)
	ylabel(input_feature_names[allinnerfeaturekeys[ind]].replace('_',' '),fontsize=40,rotation='horizontal')

#		ax.yaxis.grid(True)
#		ax.xaxis.grid(True)

#		for tick in ax.get_yticklabels():
#		    tick.set_color('y')
#		    tick.set_fontsize(10)
#		for tick in ax.get_xticklabels():
#		    tick.set_fontsize(30)
#
#	ax2 = ax.twinx()
	ax.fill_between(range(len(onefeature)),0,list(onefeature),alpha=0.3,label='template',color='b')
#	if alldistinctpeaks1.has_key(ind):
#	    for peak in alldistinctpeaks1[ind]:
#		ax.plot([peak[1],peak[1]],[0,maxval],'r')

#	if alldistinctpeaks.has_key(ind):
#	    for peak in alldistinctpeaks[ind]:
#		ax.plot([peak[1],peak[1]],[-maxval,0],'r')


#	ax.plot(range(len(onefeature)),onefeature_smooth,label='template',color='b',linewidth=2)

	ylim(-maxval,maxval)
	ymin, ymax = ylim()
	plot([N4/2,N4/2],[ymin,ymax],'k',linewidth=2)
	plot([0,N4-1],[0,0],'k',linewidth=1)

	xticks(xaxis_ticks,xaxis_labels)

	locs,labels = yticks()
	yticks(locs, map(lambda x: "%g" % x, locs))


	xlim(0,N4-1)

	ax.yaxis.grid(True)
	ax.xaxis.grid(True)

	for tick in ax.get_yticklabels():
	    tick.set_color('b')
	    tick.set_fontsize(10)


	for tick in ax.get_xticklabels():
	    tick.set_fontsize(30)



    f1.set_figwidth(40)
    f1.set_figheight(60)


    savefig('%s_curves_entirefeature_overlapping_class%s.pdf' % (outputfn,str(whichclass)))
    close()





    if False:
    #	TFpeaktoHistpeaks = [[] for i in range(3)]
	Hists = [1,3,5]
	pp = PdfPages('%s_TFpeakstoHistpeaks_relationship_class%s.pdf' % (outputfn,str(whichclass)))

	for Histind,Hist in enumerate(Hists):
	    Histprofs = feature[Histind]
	    TFpeaktoHistpeaks = {}
	    for TFind,TFprofile in enumerate(feature[0]):
		#TFpeaktoHistpeaks[Histind] = array([[TFprofile,(Histpos-TFind)*pixelsize,Histprofs[Histpos]] for Histpos in xrange(N4)])
		for Histpos in xrange(N4):
		    if not TFpeaktoHistpeaks.has_key((TFprofile,(Histpos-TFind)*pixelsize)):
			TFpeaktoHistpeaks[TFprofile,(Histpos-TFind)*pixelsize] = []
		    TFpeaktoHistpeaks[TFprofile,(Histpos-TFind)*pixelsize].append(Histprofs[Histpos])


	    keys = TFpeaktoHistpeaks.keys()
	    numbinsx = 10
	    numbinsy = 50
	    x = sorted(unique([key[0] for key in keys]))
	    y = sorted(unique([key[1] for key in keys]))

	    xbins = arange(x[0],x[-1],(x[-1]-x[0])/numbinsx)
	    ybins = arange(y[0],y[-1],(y[-1]-y[0])/numbinsy)

	    X,Y = meshgrid(xbins,ybins)
	    print X
	    print Y
	    Z = zeros((len(ybins),len(xbins)))
	    for i,xi in enumerate(x):
		for j,yj in enumerate(y):
		    if TFpeaktoHistpeaks.has_key((xi,yj)):
			Z[j,i] = mean(TFpeaktoHistpeaks[xi,yj])


	    for i,xi in enumerate(x):
		for j,yj in enumerate(y):
		    if TFpeaktoHistpeaks.has_key((xi,yj)):
			Z[j,i] = mean(TFpeaktoHistpeaks[xi,yj])



	    print Z
	    f1 = figure(1)
	    ax = p3.Axes3D(f1)
	    ax.plot_surface(X,Y,Z)
	    ax.set_xlabel('height of %s profile' % input_feature_names[0])
	    ax.set_ylabel('distance from %s profile' % input_feature_names[0])
	    ax.set_zlabel('height of %s profile' % input_feature_names[Hist])
	    f1.add_axes(ax)
	    pp.savefig()
	    close()



	pp.close()



    Hists = [1,3,5]
#   pp = PdfPages('%s_TFpeakstoHistpeaks_relationship_class%s.pdf' % (outputfn,str(whichclass)))



    TFprofiles = sorted(feature[0])
    numbinsTFprofs = 10.0
    TFprofbins = arange(TFprofiles[0],TFprofiles[-1],(TFprofiles[-1]-TFprofiles[0])/numbinsTFprofs)
    numsubplotcols = 4.0 if len(TFprofbins) / 4.0 >= 4 else 3.0 if len(TFprofbins) /3.0 >= 3 else 2.0
    numsubplotrows = ceil(len(TFprofbins)/numsubplotcols)

    f1 = figure(1)
    TFprofiles = feature[0]


    for TFindbin,TFprofileleft in enumerate(TFprofbins):
	TFprofileright = TFprofbins[TFindbin+1] if TFindbin < len(TFprofbins)-1 else sorted(TFprofiles)[-1]
	print TFprofileleft,TFprofileright
	relevant_prof_inds = [i for i in xrange(N4) if TFprofiles[i]>=TFprofileleft and TFprofiles[i]<TFprofileright]
	print [(x-N4/2)*pixelsize for x in relevant_prof_inds]
	ax = f1.add_subplot(numsubplotrows,numsubplotcols,TFindbin+1)

	for Histind,Hist in enumerate(Hists):
	    Histprofs = feature[Hist]
	    TFpeaktoHistpeaks = {}
	    for TFind,TFpos in enumerate(relevant_prof_inds):
		for Histpos in xrange(N4):
		    if not TFpeaktoHistpeaks.has_key(((Histpos-TFpos)*pixelsize)):
			TFpeaktoHistpeaks[((Histpos-TFpos)*pixelsize)] = []
		    TFpeaktoHistpeaks[((Histpos-TFpos)*pixelsize)].append(Histprofs[Histpos])

	    x = sorted(TFpeaktoHistpeaks.keys())
	    y = [mean(TFpeaktoHistpeaks[key]) for key in x]
	    ax.plot(x,y,lines[Hist],label=input_feature_names[Hist])

	legend()
	ax.set_xlabel('Relative distance from %s profile' % input_feature_names[0])
	ax.set_ylabel('Average Profile')
	title('Average Histone profiles relative to %s profiles in [%g,%g] range' % (input_feature_names[0],TFprofileleft,TFprofileright))


    f1.set_figwidth(40)
    f1.set_figheight(60)


    savefig('%s_TFpeakstoHistpeaks_relationship_class%s.pdf' % (outputfn,str(whichclass)))
    close()


#   compute_profile_relationship_plots(feature,input)




    if False:
	feature = [[mean(entirematrix_wholeouterfeature_curves[key][ind]) for ind in xrange(N4)] for key in allinnerfeaturekeys]
	feature_above = [[mean(entirematrix_wholeouterfeature_curves[key][ind])+std(entirematrix_wholeouterfeature_curves[key][ind]) for ind in xrange(N4)] for key in allinnerfeaturekeys]
	feature_below = [[mean(entirematrix_wholeouterfeature_curves[key][ind])-std(entirematrix_wholeouterfeature_curves[key][ind]) for ind in xrange(N4)] for key in allinnerfeaturekeys]



	maxval = array(feature_above).max()
	minval = array(feature_below).min()




	f1 = figure(1)
	for ind in xrange(len(feature)):


	    onefeature = feature[ind]
	    onefeature_below = feature_below[ind]
	    onefeature_above = feature_above[ind]

	    ax = f1.add_subplot(len(feature),1,ind+1)
    #				    ax.fill_between(range(len(onefeature)),0,onestate[2][ind],alpha=1,label='data',color='y')
    #				    ylim(-maxvalstates,maxvalstates)
	    ylabel(input_feature_names[allinnerfeaturekeys[ind]].replace('_',' '),fontsize=40,rotation='horizontal')

    #				    ax.yaxis.grid(True)
    #				    ax.xaxis.grid(True)

    #				    for tick in ax.get_yticklabels():
    #					    tick.set_color('y')
    #					    tick.set_fontsize(10)
    #				    for tick in ax.get_xticklabels():
    #					    tick.set_fontsize(30)
    #
#	    ax2 = ax.twinx()
	    ax.plot(range(len(onefeature)),list(onefeature),label='template',color='b',linewidth=3)
	    ax.fill_between(range(len(onefeature)),list(onefeature_below),list(onefeature_above),alpha=0.3,label='template',color='b')


	    ylim(minval,maxval)
	    ymin, ymax = ylim()
	    plot([N4/2,N4/2],[ymin,ymax],'k',linewidth=2)
	    plot([0,N4-1],[0,0],'k',linewidth=1)

	    xticks(xaxis_ticks,xaxis_labels)

	    locs,labels = yticks()
	    yticks(locs, map(lambda x: "%g" % x, locs))


	    xlim(0,N4-1)

	    ax.yaxis.grid(True)
	    ax.xaxis.grid(True)

	    for tick in ax.get_yticklabels():
		tick.set_color('b')
		tick.set_fontsize(10)


	    for tick in ax.get_xticklabels():
		tick.set_fontsize(30)



	f1.set_figwidth(40)
	f1.set_figheight(60)


	savefig('%s_curves_entirefeature_meanstdbars_overlapping_class%s.pdf' % (outputfn,str(whichclass)))
	close()






    feature = [[median(entirematrix_wholeouterfeature_curves[key][ind]) for ind in xrange(N4)] for key in allinnerfeaturekeys]
    feature_above = [[sorted(entirematrix_wholeouterfeature_curves[key][ind])[int(len(entirematrix_wholeouterfeature_curves[key][ind])*0.75)] for ind in xrange(N4)] for key in allinnerfeaturekeys]
    feature_below = [[sorted(entirematrix_wholeouterfeature_curves[key][ind])[int(len(entirematrix_wholeouterfeature_curves[key][ind])*0.25)] for ind in xrange(N4)] for key in allinnerfeaturekeys]



    maxval = array(feature_above).max()
    minval = array(feature_below).min()




    f1 = figure(1)
    for ind in xrange(len(feature)):


	onefeature = feature[ind]
	onefeature_below = feature_below[ind]
	onefeature_above = feature_above[ind]

	ax = f1.add_subplot(len(feature),1,ind+1)
#								ax.fill_between(range(len(onefeature)),0,onestate[2][ind],alpha=1,label='data',color='y')
#								ylim(-maxvalstates,maxvalstates)
	ylabel(input_feature_names[allinnerfeaturekeys[ind]].replace('_',' '),fontsize=40,rotation='horizontal')

#								ax.yaxis.grid(True)
#								ax.xaxis.grid(True)

#								for tick in ax.get_yticklabels():
#										tick.set_color('y')
#										tick.set_fontsize(10)
#								for tick in ax.get_xticklabels():
#										tick.set_fontsize(30)
#
#	ax2 = ax.twinx()
	ax.plot(range(len(onefeature)),list(onefeature),label='template',color='b',linewidth=3)
	ax.fill_between(range(len(onefeature)),list(onefeature_below),list(onefeature_above),alpha=0.3,label='template',color='b')


	ylim(minval,maxval)
	ymin, ymax = ylim()
	plot([N4/2,N4/2],[ymin,ymax],'k',linewidth=2)
	plot([0,N4-1],[0,0],'k',linewidth=1)

	xticks(xaxis_ticks,xaxis_labels)

	locs,labels = yticks()
	yticks(locs, map(lambda x: "%g" % x, locs))


	xlim(0,N4-1)

	ax.yaxis.grid(True)
	ax.xaxis.grid(True)

	for tick in ax.get_yticklabels():
	    tick.set_color('b')
	    tick.set_fontsize(10)


	for tick in ax.get_xticklabels():
	    tick.set_fontsize(30)



    f1.set_figwidth(40)
    f1.set_figheight(60)


    savefig('%s_curves_entirefeature_medianiqr_overlapping_class%s.pdf' % (outputfn,str(whichclass)))
    close()




def get_summaries(alldata):
    numfeatures = len(alldata[0])
    n1 = size(alldata[0][0],0)
    n2 = size(alldata[1][0],0)
    N = size(alldata[0][0],1)
    means = [[],[]]
#   means[0] = [[median(alldata[0][i][:,j]) for j in xrange(N)] for i in xrange(numfeatures)]
#   means[1] = [[median(alldata[1][i][:,j]) for j in xrange(N)] for i in xrange(numfeatures)]

    stds = [[],[]]
    stds[0] = [zeros(N) for i in xrange(numfeatures)]
    stds[1] = [zeros(N) for i in xrange(numfeatures)]

    means[0] = [alldata[0][i].mean(0) for i in xrange(numfeatures)]
    means[1] = [alldata[1][i].mean(0) for i in xrange(numfeatures)]

#   stds = [[],[]]
#   stds[0] = [alldata[0][i].std(0) for i in xrange(numfeatures)]
#   stds[1] = [alldata[1][i].std(0) for i in xrange(numfeatures)]



#   medians[1] = [alldata[1][i].median(0) for i in xrange(numfeatures)]


#   print len(summaries[0]),size(summaries[0][0])
#   print len(summaries[1]),size(summaries[1][0])

    return means,stds


def read_fulldata(fulldatafn):

    fulldatafile = open(fulldatafn,'r')
    data = [[],[]]
    for line in fulldatafile:
	if line[0] == '#':
	    parts = line.rstrip().rsplit(' ')
	    if parts[1] == 'numtracks':
		numfeatures = int(parts[3])
		data[0] = [[] for i in xrange(numfeatures)]
		data[1] = [[] for i in xrange(numfeatures)]


	    if parts[1] == 'featuresize':
		featuresize = int(parts[3])



	    continue

	parts = line.rstrip().rsplit(' ')
	features = [[float(x) for x in parts[1+i*featuresize:(i+1)*featuresize+1]] for i in xrange(numfeatures)]
	if parts[0] == '1':
	    for i in xrange(numfeatures):
		data[0][i].append(features[i])
	else:
	    for i in xrange(numfeatures):
		data[1][i].append(features[i])

    for i in xrange(numfeatures):
	data[0][i] = array(data[0][i])
    for i in xrange(numfeatures):
	data[1][i] = array(data[1][i])


    return data


def main(args):

    if len(args) == 0:
	print '=====INPUTS=====\ntopxpercent_top = float(args[0])\ntopxpercent_middle = float(args[1])\nbinsize = int(args[2])\nclustfeats = args[3] == "1"\nplot_all_features = args[4] == "1"\nplot_top2_analysis = args[5] == "1"\nuse_new_formula = args[6] == "1"\ndo_average_over_subfeatures = args[7] == "1"\ntypeofoutput = args[8]\ncombinetype = args[9]\ngroupsize = int(args[10])\noutputgroupsize = int(args[11])\npixelsize = int(args[12])\nsubsample1 = int(args[13])\nsortingmetric = args[14]\nfilterlim = float(args[15])\nstates_output_groupsize = int(args[16])\ninput_feature_names = args[17].rsplit(",")\noutput_fn = args[18]\nweights_fns = args[19]\nstates_fn = args[20]\ntopunitclusters = [[int(x) for x in y.rsplit(",")] if y != "all" else "all" for y in args[21].rsplit(":")]\n================='
	exit()


##
##  topxpercent_top = float(args[0])
##  topxpercent_middle = float(args[1])
##  binsize = int(args[2])
##  clustfeats = args[3] == '1'
##  plot_all_features = args[4] == '1'
##  plot_top2_analysis = args[5] == '1'
##  use_new_formula = args[6] == '1'
##  do_average_over_subfeatures = args[7] == '1'
##  typeofoutput = args[8]
##
##  whichclass = args[9] if args[9] == 'both' else int(args[9])
##  combinetype = args[10]
##  groupsize = int(args[11])
##  outputgroupsize = int(args[12])
##  pixelsize = int(args[13])
##  subsample1 = int(args[14])
##  numoffilters = int(args[15])
##  sortingmetric = args[16]
##  filterlim = float(args[17])
##  states_output_groupsize = int(args[18])
##  input_feature_names = args[19].rsplit(',')
##  output_fn = args[20]
##  weights_fns = args[21]
##  states_fn = args[22]
##  topunitclusters = [[int(x) for x in y.rsplit(',')] if y != 'all' else 'all' for y in args[23].rsplit(':')]
##  interactiongroups = [[int(x) for x in y.rsplit(',')] if y != 'all' else 'all' for y in args[24].rsplit(':')]



    topxpercent_top = float(args[0])
    plot_all_features = args[1] == '1'
    plot_top2_analysis = args[2] == '1'
    use_new_formula = args[3] == '1'

    whichclass = args[4] if args[4] == 'both' else int(args[4])
    groupsize = int(args[5])
    outputgroupsize = int(args[6])
    pixelsize = int(args[7])
    subsample1 = int(args[8])
    numoffilters = int(args[9])
    keptpeakmovepercent = float(args[10])
    keptpeakheightpercent = float(args[11])
    whichfeaturespeaks = [int(x) for x in args[12].rsplit(',')]
    input_feature_names = args[13].rsplit(',')
    output_fn = args[14]
    weights_fns = args[15]
    states_fn = args[16]
    fulldatafn = args[17]
    topunitclusters = [[int(x) for x in y.rsplit(',')] if y != 'all' else 'all' for y in args[18].rsplit(':')]
    interactiongroups = [[int(x) for x in y.rsplit(',')] if y != 'all' else 'all' for y in args[19].rsplit(':')]




#   print topunitclusters

    folders = output_fn.rsplit('/')[:-1]


    for j in xrange(1,len(folders)+1):
	folder = '/'.join(folders[:j])
	if not (folder == '' or os.path.exists(folder)):
	    os.mkdir(folder)


    topunitclusters = [sorted(x) if x != 'all' else 'all' for x in topunitclusters]

    print 'Reading features'
    features,features_subs = read_features_onefile(plot_all_features,use_new_formula,input_feature_names,weights_fns,output_fn)
    if states_fn != '0':
	print 'Reading states'
	states = read_states(states_fn)
	print 'Read %d correctly classified UP and %d DOWN genes' % (len(states[0]),len(states[1]))
	newstates = states


#	newstates = [[],[]]
#	newstates[0] = states[0][:states_output_groupsize] + states[0][(len(states[0])-states_output_groupsize)/2:(len(states[0])+states_output_groupsize)/2] + states[0][-states_output_groupsize:]
#	newstates[1] = states[1][:states_output_groupsize] + states[1][(len(states[1])-states_output_groupsize)/2:(len(states[1])+states_output_groupsize)/2] + states[1][-states_output_groupsize:]

    else:
	newstates = []

    alldata = read_fulldata(fulldatafn)
    means,stds = get_summaries(alldata)


#   print states
#   print newstates

    features.reverse()
#   output_fn += str(binsize) + '_' + str(topxpercent_top) + '_' + combinetype
    output_fn += str(topxpercent_top)



    for topunitclusterind,topunitcluster in enumerate(topunitclusters):
	print 'Cluster #', topunitclusterind , ': ', topunitcluster
	output_fn_loc = output_fn + 'geneclustind_%d' % topunitclusterind

	combinedfeatures_diff,alldistinctpeaks = analyze_top2layers(whichclass,topunitcluster,features,topxpercent_top,numoffilters,groupsize,outputgroupsize,output_fn_loc,subsample1,pixelsize,keptpeakmovepercent,keptpeakheightpercent,whichfeaturespeaks,plot_top2_analysis,plot_all_features,input_feature_names,newstates,means,stds)#,states)
	#combinedfeatures_diff = analyze_top2layers(clustfeats,whichclass,topunitcluster,features,topxpercent_top,combinetype,numoffilters,groupsize,outputgroupsize,binsize,output_fn_loc,filterlim,subsample1,pixelsize,plot_top2_analysis,plot_all_features,input_feature_names,newstates)#,states)



	#output_fn_loc += '_' + sortingmetric + str(topxpercent_middle) + ('_newform' if use_new_formula else '')
	output_fn_loc += '_' + ('_newform' if use_new_formula else '')

	#create_feature_curves_overlapping_withtoplayer_usediffs(combinedfeatures_diff,whichclass,features,topxpercent_middle,groupsize,outputgroupsize,interactiongroups,binsize,subsample1,pixelsize,output_fn_loc,input_feature_names,sortingmetric,do_average_over_subfeatures,newstates)


#	create_feature_curves_overlapping_withtoplayer_usediffs(combinedfeatures_diff,whichclass,features,groupsize,outputgroupsize,interactiongroups,subsample1,pixelsize,keptpeakmovepercent,keptpeakheightpercent,whichfeaturespeaks,output_fn_loc,input_feature_names,alldistinctpeaks,newstates)




if __name__ == "__main__": main(sys.argv[1:])








































from re import *
import sys
from operator import itemgetter
import os
import matplotlib
matplotlib.use("Agg")
from pylab import *
from numpy import *

from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D




def read_peaks(peaksfn):
    peaks = {}
    ind = -1
    for line in open(peaksfn,'r'):
	if line[0] == '=':
	    factorname = line.rstrip().rsplit(' ')[1]
	    ind += 1
	    peaks[factorname,ind] = []
	    continue

	parts = line.rstrip().rsplit(' ')
	peaks[factorname,ind].append([int(parts[0]),float(parts[1])])

    return peaks




def read_genes(genefn):
    genes = []
    for line in open(genefn,'r'):
	parts = line.rstrip().rsplit()
	genes.append([parts[1],int(parts[2]),parts[3],parts[4]])
    return genes


def read_tracks(trackfns):
    tracks = []

    for trackfn in trackfns:
	print 'Reading track file: %s' % trackfn
	if trackfn == '0':
	    tracks.append([])
	else:
	    onetrack = {}
	    starts = {}
	    spans = {}

	    trackfile = open(trackfn,'r')
	    trackfile.next()
	    trackfile.next()
	    for line in trackfile:

		if line[0:3] == 'chr':
		    if False:
			if line[0:4] not in ['chr1','chr2']:
			    break
		    parts = line.rstrip().replace('=',' ').rsplit(' ')
		    chrom = parts[0]
		    print 'Reading chromosom %s' % chrom

		    onetrack[chrom] = zeros(int(parts[6]))
		    spans[chrom] = int(parts[2])
		    print 'Span = %d' % spans[chrom]
		    starts[chrom] = int(parts[4])
		    ind = 0
		    continue

		onetrack[chrom][ind] = int(line.rstrip())
		ind += 1

	    maxval = -inf
	    for chrom in onetrack.keys():
		maxval = max(max(onetrack[chrom]),maxval)

	    for chrom in onetrack.keys():
		onetrack[chrom] = onetrack[chrom]*1.0/maxval

	    tracks.append([onetrack,spans,starts])


    return tracks


def read_trackpeaks(trackpeakfns):
    tracks = []

    for trackfn in trackpeakfns:
	print 'Reading peaks from %s' % trackfn
	if trackfn == '0':
	    tracks.append([])
	else:
	    onetrack = {}

	    for line in open(trackfn,'r').readlines()[1:]:
		parts = line.rstrip().rsplit('\t')
		if False:
		    if parts[0] not in ['chr1','chr2']:
			break
		if onetrack.has_key(parts[0]):
		    onetrack[parts[0]].append([(int(parts[1])+int(parts[2]))/2,int(parts[1]),int(parts[2]),float(parts[4])])
		else:
		    onetrack[parts[0]] = [[(int(parts[1])+int(parts[2]))/2,int(parts[1]),int(parts[2]),float(parts[4])]]

	    onetrack_rightsort = {}
	    for key in onetrack.keys():
		onetrack[key].sort()
		onetrack_rightsort[key] = sorted(onetrack[key],key=itemgetter(1),reverse=True)


	    tracks.append([onetrack,onetrack_rightsort])

    return tracks




def return_count_agg(loca,locb,counts,dataaggregation):
    vals = counts[min(max(loca,0),len(counts)):min(max(locb,0),len(counts))]
    if len(vals) == 0:
	return 0

    else:
	if dataaggregation == 'average':
	    return vals.mean()
	elif dataaggregation == 'max':
	    return vals.max()
	elif dataaggregation == 'sum':
	    return vals.sum()



def fast_forward_peaks(loc,startind,peaks):

#   print startind,len(peaks)
    peakcenter = (peaks[startind][1]+peaks[startind][0])/2
    smallestdistancetopeak = abs(peakcenter - loc)
    smallestdistancepeakind = startind
    for ind in range(startind+1,len(peaks)):
	peakcenter = (peaks[ind][1]+peaks[ind][0])/2
	distancetopeak = abs(peakcenter-loc)
	if distancetopeak < smallestdistancetopeak:
	    smallestdistancetopeak = distancetopeak
	    smallestdistancepeakind = ind
	else:
	    return [smallestdistancetopeak,smallestdistancepeakind,peakcenter]

    return [smallestdistancetopeak,smallestdistancepeakind,peakcenter]




def fast_forward_featurepeaks(loc,startind,featurepeaks):

    smallestdistancetopeak = abs(featurepeaks[startind][0] - loc)
    smallestdistancepeakind = startind
    for ind in range(startind+1,len(featurepeaks)):
#	peakcenter = (peaks[ind][1]+peaks[ind][0])/2
	distancetopeak = abs(featurepeaks[ind][0]-loc)
	if distancetopeak <= smallestdistancetopeak:
	    smallestdistancetopeak = distancetopeak
	    smallestdistancepeakind = ind
	else:
	    return [smallestdistancetopeak,smallestdistancepeakind]

    return [smallestdistancetopeak,smallestdistancepeakind]






def fast_forward_genes(loc,startind,genes):

    smallestdistancetopeak = abs(genes[startind][0] - loc)
    smallestdistancepeakind = startind
    for ind in range(startind+1,len(genes)):
	distancetopeak = abs(genes[ind][0]-loc)
	if distancetopeak <= smallestdistancetopeak:
	    smallestdistancetopeak = distancetopeak
	    smallestdistancepeakind = ind
	else:
	    return [smallestdistancetopeak,smallestdistancepeakind]

    return [smallestdistancetopeak,smallestdistancepeakind]






def store_feature_matrixes(pixel_matrixes,skipped_pixels,genes,feature_fn, pixel_size, window_size, genefn, tracksfn):

    featurefile = open(feature_fn,'w')
    print >> featurefile, '# windowsize = %d\n# pixelsize = %d\n# numtracks = %d\n# featuresize = %d' % (window_size,pixel_size,len(pixel_matrixes),size(pixel_matrixes[0],1))
    print >> featurefile, '# genes = %s\n# tracks = %s' % (genefn, '; '.join(tracksfn))
    print >> featurefile, '# skipped_pixels = %s' % (' | '.join([','.join([str(x) for x in skipped_pixels_pertrack]) for skipped_pixels_pertrack in skipped_pixels]))
    for geneind,gene in enumerate(genes):
	print >> featurefile, gene[3],
	for trackind,track_matrix in enumerate(pixel_matrixes):
	    for ind,feature in enumerate(track_matrix[geneind,:]):
		if feature == int(feature):
		    print >> featurefile, int(feature),

		else:
		    print >> featurefile, feature,

	print >> featurefile

    featurefile.close()




def create_agregate_profiles_old(genes,tracks,trackpeaks,peaks,peakwidth,window_size,refFactorID):

    refPeaks = peaks[refFactorID]

    allpeakprofiles = []
    for peak in refPeaks:

	offset = peak[0]
	locations = []
	for gene in genes:
	    tss = gene[1]
	    chrom = gene[0]
	    orient = gene[2]

	    if orient == '+':
		locations.append(['+',chrom,tss + offset])
	    else:
		locations.append(['-',chrom,tss - offset])

	profiles = [[] for i in range(numtracks)]
	for trackind in range(numtracks):
	    if tracks[trackind] == []:
    #		print 'adding blanks in track %d' % trackind
		profiles[trackind] = []
		continue

	    counts = tracks[trackind][0][chrom]
	    span = tracks[trackind][1][chrom]
	    start = tracks[trackind][2][chrom]



	    profiles[trackind] = []#zeros([len(locations),(window_size*2)/span])

	    trackpeakind_pos = 0
	    trackpeakind_neg = 0

#	    curlocind = 0
	    deltas = []
	    for locind,location in enumerate(locations):

		chrom = location[1]
		orient = location[0]
		basepos = location[2]
		if not tracks[trackind][0].has_key(chrom) or not trackpeaks[refFactorID[1]][0].has_key(chrom):
#		    print 'no chrom %s in track %d or in track %d' % (chrom,trackind,refFactorID[1])
#		    pixel_matrixes[trackind][locind,:] = -1
		    continue


		peakarea = [basepos - peakwidth,basepos + peakwidth]

		if orient == '+':

		    winstart = basepos - window_size
		    winend = basepos + window_size

		    [trackpeakind_pos,pixval] = fast_forward_peaks_pos(peakarea[0],peakarea[1],trackpeakind_pos,trackpeaks[refFactorID[1]][0][chrom])
#		    pixval = 1
		    if pixval != 0:
			oneprofile = zeros((window_size*2)/span)
			for pixind,[loca,locb] in enumerate([[a,a+span] for a in range(winstart,winend,span)]):
			    pixval = return_count_agg((loca-start)/span,(locb-start)/span,counts,'sum')
			    oneprofile[pixind] = pixval


			profiles[trackind].append(oneprofile)
		else:

		    winstart = basepos + window_size
		    winend = basepos - window_size

		    [trackpeakind_neg,pixval] = fast_forward_peaks_pos(peakarea[0],peakarea[1],trackpeakind_neg,trackpeaks[refFactorID[1]][1][chrom])
#		    pixval = 1
		    if pixval != 0:
			oneprofile = zeros((window_size*2)/span)
			for pixind,[loca,locb] in enumerate([[a-span,a] for a in range(winstart,winend,-span)]):
			    pixval = return_count_agg((loca-start)/span,(locb-start)/span,counts,'sum')
			    oneprofile[pixind] = pixval

			profiles[trackind].append(oneprofile)

	    profiles[trackind] = array(profiles[trackind])


	allpeakprofiles.append(profiles)

    return allpeakprofiles




def create_agregate_profiles_old1(genes,tracks,trackpeaks,peaks,peakwidth,window_size,refFactorID,use_orientation):

    refPeaks = peaks[refFactorID]

    allpeakprofiles = []
    smallestdistancestopeak = [[] for i in xrange(len(refPeaks))]

    allchroms = list(unique([x[0] for x in genes]))

    for peakind,peak in enumerate(refPeaks):

	offset = peak[0]
	locations = []
	profiles = [[] for i in range(numtracks)]

	trackpeakind = dict([(chrom,0) for chrom in allchroms])

	for gene in genes:
	    tss = gene[1]
	    chrom = gene[0]
	    orient = gene[2]


	    basepos = tss + offset if orient == '+' else tss - offset

	    if not trackpeaks[refFactorID[1]][0].has_key(chrom):
#		print 'no chrom %s in reference track %s' % (chrom,refFactorID[0])
		continue


#	    print chrom,tss,orient,


	    #peakarea = [basepos - peakwidth,basepos + peakwidth]


	    [smallestdistancetopeak,trackpeakind[chrom],peakcenter] = fast_forward_peaks(basepos,trackpeakind[chrom],trackpeaks[refFactorID[1]][0][chrom])
#	    print basepos,chrom,peak[0],peak[1],smallestdistancetopeak
	    smallestdistancestopeak[peakind].append(smallestdistancetopeak)

	    if smallestdistancetopeak > peakwidth:
#		print 'no %s peaks in the current peak area' % refFactorID[0]
		continue


	    winstart_untransformed = peakcenter - window_size
	    winend_untransformed = peakcenter + window_size


	    for trackind in range(numtracks):
		if tracks[trackind] == [] or not tracks[trackind][0].has_key(chrom):
#		    print 'adding blanks in track %d' % trackind
		    continue


		counts = tracks[trackind][0][chrom]
		span = tracks[trackind][1][chrom]
		start = tracks[trackind][2][chrom]

		oneprofile = zeros(window_size/span*2)
		winstart = (winstart_untransformed-start)/span
		winend = (winend_untransformed-start)/span

		if orient == '+' or not use_orientation:
		    vals = counts[min(max(winstart,0),len(counts)):min(max(winend,0),len(counts))]
		    if len(vals) != 0:
			profilestartind = max(0,-winstart)
			profileendind = len(oneprofile)-max(0,winend-(len(counts)))-1
#			print len(vals), profilestartind,profileendind, winend,winstart,len(counts),len(oneprofile),type(oneprofile),type(vals)


			oneprofile[profilestartind:profileendind+1] = vals

			profiles[trackind].append(oneprofile)
		else:
		    vals = array(list(reversed(counts[min(max(winstart,0),len(counts)):min(max(winend,0),len(counts))])))
		    if len(vals) != 0:

			profilestartind = max(0,winend-(len(counts)))
			profileendind = len(oneprofile)-max(0,-winstart)-1

#			print len(vals), profilestartind,profileendind, winend,winstart,len(counts),len(oneprofile),type(oneprofile),type(vals)
			oneprofile[profilestartind:profileendind+1] = vals

			profiles[trackind].append(oneprofile)

#		print 'len of profiles[%d]' % (trackind),len(profiles[trackind])


	profiles = [array(profiles[trackind]) for trackind in range(numtracks)]
#	print 'len of profiles', len(profiles)
	allpeakprofiles.append(profiles)

    return allpeakprofiles,smallestdistancestopeak





def create_agregate_profiles(genes,tracks,trackpeaks,peaks,peakwidth,window_size,refFactorID,use_orientation,geneclass_widths):

    refPeaks = peaks[refFactorID]

    allpeakprofiles = []
    smallestdistancestopeak = [[] for i in xrange(len(refPeaks))]

    allchroms = list(unique([x[0] for x in genes]))




    all_ref_peaks = {}
    for peakind,peak in enumerate(refPeaks):

	offset = peak[0]
	for gene in genes:
	    tss = gene[1]
	    chrom = gene[0]
	    orient = gene[2]
	    if not all_ref_peaks.has_key(chrom):
		all_ref_peaks[chrom] = []


	    basepos = tss + offset if orient == '+' else tss - offset

	    all_ref_peaks[chrom].append([basepos,int(sign(int(gene[3]))),int(sign(peak[1])),orient])

    profiles_geneclasspeakclass = [[[] for i in range(numtracks)],[[] for i in range(numtracks)]]
    profiles_peakclass = [[[] for i in range(numtracks)],[[] for i in range(numtracks)]]
    profiles_geneclass = [[[[] for i in range(numtracks)],[[] for i in range(numtracks)]] for i in geneclass_widths]

#   profiles_geneclasspeakclass = [[],[]]
#   profiles_peakclass = [[],[]]
#   profiles_geneclass = [[[],[]] for i in geneclass_widths]



    for chrom in trackpeaks[refFactorID[1]][0]:
	trackpeakind = 0
	all_ref_peaks[chrom].sort()
#	print all_ref_peaks[chrom]

	for peakind,realpeak in enumerate(trackpeaks[refFactorID[1]][0][chrom]):

	    if not all_ref_peaks.has_key(chrom):
		continue

	    peakcenter = realpeak[0]#(realpeak[1]+realpeak[1])/2
	    [smallestdistancetopeak,trackpeakind] = fast_forward_featurepeaks(peakcenter,trackpeakind,all_ref_peaks[chrom])

#	    print peakcenter,smallestdistancetopeak,trackpeakind,all_ref_peaks[chrom][trackpeakind][0],all_ref_peaks[chrom][trackpeakind+1][0]
	    orient = all_ref_peaks[chrom][trackpeakind][-1]

#	    smallestdistancestopeak[peakind].append(smallestdistancetopeak)

	    if smallestdistancetopeak > peakwidth:
		continue


	    if all_ref_peaks[chrom][trackpeakind][1] == all_ref_peaks[chrom][trackpeakind][2]:
		profiles = profiles_geneclasspeakclass[0] if all_ref_peaks[chrom][trackpeakind][1]>0 and all_ref_peaks[chrom][trackpeakind][2]>0 else profiles_geneclasspeakclass[1]
		profiles = compute_profiles_sub(profiles,peakcenter,window_size,tracks,chrom,orient,use_orientation)


	    profiles = profiles_peakclass[0] if all_ref_peaks[chrom][trackpeakind][2]>0 else profiles_peakclass[1]
	    profiles = compute_profiles_sub(profiles,peakcenter,window_size,tracks,chrom,orient,use_orientation)


    for trackind in range(numtracks):
	profiles_peakclass[0][trackind] = array(profiles_peakclass[0][trackind])
	profiles_peakclass[1][trackind] = array(profiles_peakclass[1][trackind])

	profiles_geneclasspeakclass[0][trackind] = array(profiles_geneclasspeakclass[0][trackind])
	profiles_geneclasspeakclass[1][trackind] = array(profiles_geneclasspeakclass[1][trackind])


    print 'Gene & Feature types:'
    print '  ENHANCER: %d profiles' % size(profiles_geneclasspeakclass[0][0],0)
    print '  INHIBITOR: %d profiles' % size(profiles_geneclasspeakclass[1][0],0)

    print 'Feature types:'
    print '  ENHANCER: %d profiles' % size(profiles_peakclass[0][0],0)
    print '  INHIBITOR: %d profiles' % size(profiles_peakclass[1][0],0)





    genes_by_chrom = {}
    for gene in genes:
	if not genes_by_chrom.has_key(gene[0]):
	    genes_by_chrom[gene[0]] = []
	genes_by_chrom[gene[0]].append([gene[1],gene[2],gene[3]])


    for chrom in trackpeaks[refFactorID[1]][0]:
	genes_by_chrom[chrom].sort()
	trackpeakind = 0
	for peakind,realpeak in enumerate(trackpeaks[refFactorID[1]][0][chrom]):
	    peakcenter = realpeak[0]#(realpeak[0]+realpeak[1])/2
	    [smallestdistancetopeak,trackpeakind] = fast_forward_genes(peakcenter,trackpeakind,genes_by_chrom[chrom])

	    for geneclass_width_ind,geneclass_width in enumerate(geneclass_widths):
		if smallestdistancetopeak > geneclass_width:
		    continue

		orient = genes_by_chrom[chrom][trackpeakind][1]
		profiles = profiles_geneclass[geneclass_width_ind][0] if int(genes_by_chrom[chrom][trackpeakind][2])>0 else profiles_geneclass[geneclass_width_ind][1]
		profiles = compute_profiles_sub(profiles,peakcenter,window_size,tracks,chrom,orient,use_orientation)


    for geneclass_width_ind,geneclass_width in enumerate(geneclass_widths):
	for trackind in range(numtracks):
	    profiles_geneclass[geneclass_width_ind][0][trackind] = array(profiles_geneclass[geneclass_width_ind][0][trackind])
	    profiles_geneclass[geneclass_width_ind][1][trackind] = array(profiles_geneclass[geneclass_width_ind][1][trackind])


	print 'Gene types (%d window around TSS):' % geneclass_width
	print '	 ENHANCER: %d profiles' % size(profiles_geneclass[geneclass_width_ind][0][0],0)
	print '	 INHIBITOR: %d profiles' % size(profiles_geneclass[geneclass_width_ind][1][0],0)





    return profiles_geneclasspeakclass,profiles_peakclass,profiles_geneclass,smallestdistancestopeak




def compute_profiles_sub(profiles,peakcenter,windows_size,tracks,chrom,orient,use_orientation):
    winstart_untransformed = peakcenter - window_size
    winend_untransformed = peakcenter + window_size

    numtracks = len(tracks)
    span = 10
    for trackind in range(numtracks):
	oneprofile = zeros(window_size/span*2)


	if tracks[trackind] == [] or not tracks[trackind][0].has_key(chrom):
	    profiles[trackind].append(oneprofile)
	    continue

	counts = tracks[trackind][0][chrom]
	span = tracks[trackind][1][chrom]
	start = tracks[trackind][2][chrom]

	winstart = (winstart_untransformed-start)/span
	winend = (winend_untransformed-start)/span

	if orient == '+' or not use_orientation:
	    vals = counts[min(max(winstart,0),len(counts)):min(max(winend,0),len(counts))]
	    if len(vals) != 0:
		profilestartind = max(0,-winstart)
		profileendind = len(oneprofile)-max(0,winend-(len(counts)))-1

		oneprofile[profilestartind:profileendind+1] = vals

	    profiles[trackind].append(oneprofile)

	else:
	    vals = array(list(reversed(counts[min(max(winstart,0),len(counts)):min(max(winend,0),len(counts))])))
	    if len(vals) != 0:

		profilestartind = max(0,winend-(len(counts)))
		profileendind = len(oneprofile)-max(0,-winstart)-1

		oneprofile[profilestartind:profileendind+1] = vals

	    profiles[trackind].append(oneprofile)



    return profiles



#def store_peakprofiles(outputfn,allpeakprofiles_geneclasspeakclass,allpeakprofiles_peakclass,allpeakprofiles_geneclass):



def plot_all_peakprofiles(peaks,allpeakprofiles,smallestdistancestopeak,profile_summaries,refFactorID,tracknames,outputfn,window_size,span=10,outputgroupsize = 2500,bin=4000):

    lines = ['k-','b-','b--','r-','r--','g-','g--','m-','m--','c-','c--']


    pp = PdfPages('%s_ref%s.pdf' % (outputfn,refFactorID[0]))






#   print 'len of allpeakprofiles', len(allpeakprofiles)
    for profileind,peakprofile in enumerate(allpeakprofiles):
	f1 = figure(1)
	numtracks = len(peakprofile)

#	print 'len of peakprofile' , len(peakprofile)
#	print peakprofile[0],type(peakprofile[0])

	for trackind in xrange(numtracks):

	    if len(peakprofile[trackind]) == 0:
		continue


	    N4 = size(peakprofile[trackind],1)
    #	    print N4
	    xaxis_ticks = list(sorted(-array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)
	    xaxis_labels = ['%s' % ('' if x == 0 else str(x*span/1000.0) + 'Kb' if abs(x*span) >= 1000 and abs(x*span) < 1000000 else (str(x*span/1000000.0)+'Mb' if abs(x*span) >= 1000000 else str(x*span)+'b')) for x in list(array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])-N4/2) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0]))]

    #	    print xaxis_ticks
    #	    print xaxis_labels

	    averageprofile = peakprofile[trackind].mean(0)
#	    stdevprofile = peakprofile[trackind].std(0)

	    #medianprofile = peakprofile[trackind].median(0)
	    #profile_above = [sorted(peakprofile[trackind][:,i])[int(size(peakprofile[trackind],0))*0.75] for i in xrange(size(peakprofile[trackind],1))]
	    #profile_below = [sorted(peakprofile[trackind][:,i])[int(size(peakprofile[trackind],0))*0.25] for i in xrange(size(peakprofile[trackind],1))]



	    ax = f1.add_subplot(111)
	    ax.plot(averageprofile,lines[trackind],label=tracknames[trackind])



	    #ax = f1.add_subplot(212)
	    #ax.plot(averageprofile,lines[trackind],label=tracknames[trackind])
	    xticks(xaxis_ticks,xaxis_labels)




	xlabel('Offset from %s peaks (bases)' % refFactorID[0])
	ylabel('Average profiles')
	legend()


	title('TF & epigenetic enrichment profiles relative to %s peak [%d,%g]' % (refFactorID[0],peaks[refFactorID][profileind][0],peaks[refFactorID][profileind][1]))
	pp.savefig()
	close()

#	f1 = figure(1)
#	hist(log(array([x for x in smallestdistancestopeak[profileind] if x != 0])),bins=20)
#	title('Histogram of distances of %s peaks to %s features for peak [%d,%g]' % (refFactorID[0],refFactorID[0],peaks[refFactorID][profileind][0],peaks[refFactorID][profileind][1]))
#	f1.set_figwidth(40)
#	f1.set_figheight(50)
#	pp.savefig()
#	close()





    f1 = figure()
    ax = f1.add_subplot(111, projection='3d')
    n = 100


    max_peakheights = max([abs(peaks[refFactorID][profileind][1]) for profileind in xrange(len(allpeakprofiles))])
    min_peakdiststotss = min([peaks[refFactorID][profileind][0] for profileind in xrange(len(allpeakprofiles))])
    min_peaktolocdists = min([mean(smallestdistancestopeak[profileind]) for profileind in xrange(len(allpeakprofiles))])

    peakdiststotss = [peaks[refFactorID][profileind][0] for profileind in xrange(len(allpeakprofiles)) if peaks[refFactorID][profileind][1]>0]
    peakheights = [abs(peaks[refFactorID][profileind][1]) for profileind in xrange(len(allpeakprofiles)) if peaks[refFactorID][profileind][1]>0]
    peaktolocdists = [mean(smallestdistancestopeak[profileind]) for profileind in xrange(len(allpeakprofiles)) if peaks[refFactorID][profileind][1]>0]

    n = len(peakdiststotss)


    ax.scatter(peakdiststotss, peakheights, peaktolocdists,c='r')
    ax.plot(peakdiststotss,ones(n)*max_peakheights, peaktolocdists,'k.')#,markersize=0.6)
    ax.plot(ones(n)*min_peakdiststotss,peakheights, peaktolocdists,'k.')#,markersize=0.6)
    ax.plot(peakdiststotss,peakheights, ones(n)*min_peaktolocdists,'k.')#,markersize=0.6)



    peakdiststotss = [peaks[refFactorID][profileind][0] for profileind in xrange(len(allpeakprofiles)) if peaks[refFactorID][profileind][1]<0]
    peakheights = [abs(peaks[refFactorID][profileind][1]) for profileind in xrange(len(allpeakprofiles)) if peaks[refFactorID][profileind][1]<0]
    peaktolocdists = [mean(smallestdistancestopeak[profileind]) for profileind in xrange(len(allpeakprofiles)) if peaks[refFactorID][profileind][1]<0]


    n = len(peakdiststotss)


    ax.scatter(peakdiststotss, peakheights, peaktolocdists,c='b')
    ax.plot(peakdiststotss,ones(n)*max_peakheights, peaktolocdists,'k.')#,markersize=0.6)
    ax.plot(ones(n)*min_peakdiststotss,peakheights, peaktolocdists,'k.')#,markersize=0.6)
    ax.plot(peakdiststotss,peakheights, ones(n)*min_peaktolocdists,'k.')#,markersize=0.6)


    ax.set_xlabel('distance to tss')
    ax.set_ylabel('relevance')
    ax.set_zlabel('proximity of peak to feature')


    f1.set_figwidth(15)
    f1.set_figheight(15)


    pp.savefig()
    close


    f1 = figure()

    numtracks = len(peakprofile)

    if numtracks == 1:
	numsubplotcols = 1
	numsubplotrows = 1
    else:
	numsubplotcols = 4.0 if numtracks / 4.0 >= 4 else 3.0 if numtracks/3.0 >= 3 else 2.0
	numsubplotrows = ceil(numtracks/numsubplotcols)


    for trackind in xrange(numtracks):
	ax = f1.add_subplot(numsubplotrows,numsubplotcols,trackind+1,projection='3d')
	peakdiststotss = [profile_summary[0] for profile_summary in profile_summaries if profile_summary[2][trackind] != -1]
	peakheights = [profile_summary[1] for profile_summary in profile_summaries if profile_summary[2][trackind] != -1]
	binsummary = [profile_summary[2][trackind] for profile_summary in profile_summaries if profile_summary[2][trackind] != -1]

	n = len(peakdiststotss)

#	ax.scatter(peakdiststotss, peakheights, binsummary,'r')
	ax.plot(peakdiststotss,ones(n)*max(peakheights), binsummary,'r.')#,markersize=0.6)
	params=polyfit(peakdiststotss,binsummary,2)
	binsummary_hat=polyval(params,peakdiststotss)
	ax.plot(peakdiststotss,ones(n)*max(peakheights), binsummary_hat,'r--',linewidth=0.6)#,markersize=0.6)
	params=polyfit(peakdiststotss,binsummary,1)
	binsummary_hat=polyval(params,peakdiststotss)
	ax.plot(peakdiststotss,ones(n)*max(peakheights), binsummary_hat,'g--',linewidth=0.6)#,markersize=0.6)


	ax.plot(ones(n)*min(peakdiststotss),peakheights, binsummary,'g.')#,markersize=0.6)
	params=polyfit(peakheights,binsummary,2)
	binsummary_hat=polyval(params,sorted(peakheights))
	ax.plot(ones(n)*min(peakdiststotss),sorted(peakheights), binsummary_hat,'r--',linewidth=0.6)#,markersize=0.6)
	params=polyfit(peakheights,binsummary,1)
	binsummary_hat=polyval(params,sorted(peakheights))
	ax.plot(ones(n)*min(peakdiststotss),sorted(peakheights), binsummary_hat,'g--',linewidth=0.6)#,markersize=0.6)


	ax.plot(peakdiststotss,peakheights, ones(n)*min(binsummary),'b.')#,markersize=0.6)



	ax.set_xlabel('distance to tss')
	ax.set_ylabel('relevance')
	ax.set_zlabel('%s height [within %d of peak]' % (tracknames[trackind],bin))


    f1.set_figwidth(20)
    f1.set_figheight(20)


    pp.savefig()
    close




    pp.close()





def plot_all_peakprofiles_twoclasses(allpeakprofiles_geneclasspeakclass,allpeakprofiles_peakclass,allpeakprofiles_geneclass,smallestdistancestopeak,refFactorID,tracknames,outputfn,window_size,geneclass_widths,span=10,outputgroupsize = 2500):

#   lines = ['k-','k--','b-','b--','r-','r--','g-','g--','m-','m--','c-','c--','y-','y--']

#   cols = ['k','b','r','g','y','c','m']
#   lines = ['-','--']
#   lines = ['k-','b-','b--','r-','r--','g-','g--','m-','m--','c-','c--']
    lines = ['k-','b-','g-','r-','m-','c-','y-','k--','r--','b--','g--']


#   styles =
    x = range(cm.Paired.N)
    shuffle(x)


    pp = PdfPages('%s_ref%s.pdf' % (outputfn,refFactorID[0]))

    f1 = figure(1)
    ax = f1.add_subplot(111)
    numtracks = len(allpeakprofiles_geneclasspeakclass[0])
    ind = 0

    #cols.extend(cm.Paired(x[:numtracks-len(cols)]))
    #print cols
#   print allpeakprofiles_geneclasspeakclass

    for trackind in xrange(numtracks):
	if len(allpeakprofiles_geneclasspeakclass[0][trackind]) > 0:
	    ax.plot(allpeakprofiles_geneclasspeakclass[0][trackind]-allpeakprofiles_geneclasspeakclass[1][trackind],lines[ind],label='ENHANCER - IHIBITOR %s' % tracknames[trackind])
#	    ax.plot(allpeakprofiles_geneclasspeakclass[0][trackind]-allpeakprofiles_geneclasspeakclass[1][trackind],c=cols[ind],linestyle='-',label='ENHANCER - IHIBITOR`%s' % tracknames[trackind])

	ind += 1

#   ind = 0
#   numtracks = len(allpeakprofiles_geneclasspeakclass[1])
#   for trackind in xrange(numtracks):
#	if len(allpeakprofiles_geneclasspeakclass[1][trackind]) > 0:
#	    ax.plot(allpeakprofiles_geneclasspeakclass[1][trackind],c=cols[ind],linestyle='--',label='INHIBITOR %s' % tracknames[trackind])
#	ind += 1


    N4 = len(allpeakprofiles_geneclasspeakclass[1][0])
    xaxis_ticks = list(sorted(-array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)
    xaxis_labels = ['%s' % ('' if x == 0 else str(x*span/1000.0) + 'Kb' if abs(x*span) >= 1000 and abs(x*span) < 1000000 else (str(x*span/1000000.0)+'Mb' if abs(x*span) >= 1000000 else str(x*span)+'b')) for x in list(array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])-N4/2) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0]))]

    ymin, ymax = ylim()
    plot([N4/2,N4/2],[ymin,ymax],'k',linewidth=1)
    xticks(xaxis_ticks,xaxis_labels)
    xlim(0,N4-1)

    ax.yaxis.grid(True)
    ax.xaxis.grid(True)

    for tick in ax.get_xticklabels():
	tick.set_fontsize(15)



    xlabel('Offset from %s peaks (bases)' % refFactorID[0])
    ylabel('Average profiles')
    legend(ncol=2)



    title('TF & epigenetic enrichment profiles relative to %s peaks for classes of peaks based on Gene & Feature types' % (refFactorID[0]))
    f1.set_figwidth(25)
    f1.set_figheight(25)

    pp.savefig()
    close()


    f1 = figure(1)
    ax = f1.add_subplot(111)
    numtracks = len(allpeakprofiles_peakclass[0])
    ind = 0
    for trackind in xrange(numtracks):
	if len(allpeakprofiles_peakclass[0][trackind]) > 0:
	    ax.plot(allpeakprofiles_peakclass[0][trackind]-allpeakprofiles_peakclass[1][trackind],lines[ind],label='ENHANCER - INHIBITOR %s' % tracknames[trackind])
#	    ax.plot(allpeakprofiles_peakclass[0][trackind]-allpeakprofiles_peakclass[1][trackind],c=cols[ind],linestyle='-',label='ENHANCER - INHIBITOR %s' % tracknames[trackind])

	ind += 1

#   ind = 0
#   numtracks = len(allpeakprofiles_peakclass[1])
#   for trackind in xrange(numtracks):
#	if len(allpeakprofiles_peakclass[1][trackind]) > 0:
#	    ax.plot(allpeakprofiles_peakclass[1][trackind],c=cols[ind],linestyle='--',label='INHIBITOR %s' % tracknames[trackind])
#	ind += 1


    N4 = len(allpeakprofiles_peakclass[1][0])
    xaxis_ticks = list(sorted(-array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)
    xaxis_labels = ['%s' % ('' if x == 0 else str(x*span/1000.0) + 'Kb' if abs(x*span) >= 1000 and abs(x*span) < 1000000 else (str(x*span/1000000.0)+'Mb' if abs(x*span) >= 1000000 else str(x*span)+'b')) for x in list(array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])-N4/2) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0]))]

    ymin, ymax = ylim()
    plot([N4/2,N4/2],[ymin,ymax],'k',linewidth=1)
    xticks(xaxis_ticks,xaxis_labels)
    xlim(0,N4-1)

    ax.yaxis.grid(True)
    ax.xaxis.grid(True)

    for tick in ax.get_xticklabels():
	tick.set_fontsize(15)


    xlabel('Offset from %s peaks (bases)' % refFactorID[0])
    ylabel('Average profiles')
    legend(ncol=2)



    title('TF & epigenetic enrichment profiles relative to %s peaks for classes of peaks based on Feature types' % (refFactorID[0]))
    f1.set_figwidth(25)
    f1.set_figheight(25)

    pp.savefig()
    close()



    for geneclass_width_ind in xrange(len(geneclass_widths)):
#	print allpeakprofiles_geneclass[geneclass_width_ind]
	f1 = figure(1)
	ax = f1.add_subplot(111)
	numtracks = len(allpeakprofiles_geneclass[geneclass_width_ind][0])

	#print 'numtracks', numtracks
	ind = 0
	for trackind in xrange(numtracks):

	    if len(allpeakprofiles_geneclass[geneclass_width_ind][0][trackind]) > 0:
		ax.plot(allpeakprofiles_geneclass[geneclass_width_ind][0][trackind]-allpeakprofiles_geneclass[geneclass_width_ind][1][trackind],lines[ind],label='ENHANCER - INHIBITOR %s' % tracknames[trackind])
#		ax.plot(allpeakprofiles_geneclass[geneclass_width_ind][0][trackind],c=cols[ind],linestyle='-',label='ENHANCER - INHIBITOR %s' % tracknames[trackind])

	    ind += 1

#	numtracks = len(allpeakprofiles_geneclass[geneclass_width_ind][1])
#	ind = 0

#	for trackind in xrange(numtracks):
#	    if len(allpeakprofiles_geneclass[geneclass_width_ind][1][trackind]) > 0:
#		ax.plot(allpeakprofiles_geneclass[geneclass_width_ind][1][trackind],c=cols[ind],linestyle='--',label='INHIBITOR %s' % tracknames[trackind])
#	    ind += 1


	N4 = len(allpeakprofiles_geneclass[geneclass_width_ind][1][0])
	xaxis_ticks = list(sorted(-array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)
	xaxis_labels = ['%s' % ('' if x == 0 else str(x*span/1000.0) + 'Kb' if abs(x*span) >= 1000 and abs(x*span) < 1000000 else (str(x*span/1000000.0)+'Mb' if abs(x*span) >= 1000000 else str(x*span)+'b')) for x in list(array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])-N4/2) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0]))]

	ymin, ymax = ylim()
	plot([N4/2,N4/2],[ymin,ymax],'k',linewidth=1)
	xticks(xaxis_ticks,xaxis_labels)
	xlim(0,N4-1)

	ax.yaxis.grid(True)
	ax.xaxis.grid(True)

	for tick in ax.get_xticklabels():
	    tick.set_fontsize(15)

	xlabel('Offset from %s peaks (bases)' % refFactorID[0])
	ylabel('Average profiles')
	legend(ncol=2)


	title('TF & epigenetic enrichment profiles relative to %s peaks for classes of peaks based on Gene class in region of %d around TSS' % (refFactorID[0],geneclass_widths[geneclass_width_ind]))
	f1.set_figwidth(25)
	f1.set_figheight(25)

	pp.savefig()
	close()





    f1 = figure(1)
    ax = f1.add_subplot(111)
    numtracks = len(allpeakprofiles_geneclasspeakclass[0])
    ind = 0

    #cols.extend(cm.Paired(x[:numtracks-len(cols)]))
    #print cols
#	print allpeakprofiles_geneclasspeakclass

    for trackind in xrange(numtracks):
	if len(allpeakprofiles_geneclasspeakclass[0][trackind]) > 0:
	    ax.plot(log(allpeakprofiles_geneclasspeakclass[0][trackind]/allpeakprofiles_geneclasspeakclass[1][trackind]),lines[ind],label='log(ENHANCER/IHIBITOR) %s' % tracknames[trackind])
#			ax.plot(allpeakprofiles_geneclasspeakclass[0][trackind]-allpeakprofiles_geneclasspeakclass[1][trackind],c=cols[ind],linestyle='-',label='ENHANCER - IHIBITOR`%s' % tracknames[trackind])

	ind += 1

#	ind = 0
#	numtracks = len(allpeakprofiles_geneclasspeakclass[1])
#	for trackind in xrange(numtracks):
#		if len(allpeakprofiles_geneclasspeakclass[1][trackind]) > 0:
#			ax.plot(allpeakprofiles_geneclasspeakclass[1][trackind],c=cols[ind],linestyle='--',label='INHIBITOR %s' % tracknames[trackind])
#		ind += 1


    N4 = len(allpeakprofiles_geneclasspeakclass[1][0])
    xaxis_ticks = list(sorted(-array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)
    xaxis_labels = ['%s' % ('' if x == 0 else str(x*span/1000.0) + 'Kb' if abs(x*span) >= 1000 and abs(x*span) < 1000000 else (str(x*span/1000000.0)+'Mb' if abs(x*span) >= 1000000 else str(x*span)+'b')) for x in list(array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])-N4/2) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0]))]

    ymin, ymax = ylim()
    plot([N4/2,N4/2],[ymin,ymax],'k',linewidth=1)
    xticks(xaxis_ticks,xaxis_labels)
    xlim(0,N4-1)

    ax.yaxis.grid(True)
    ax.xaxis.grid(True)

    for tick in ax.get_xticklabels():
	tick.set_fontsize(15)



    xlabel('Offset from %s peaks (bases)' % refFactorID[0])
    ylabel('Average profiles')
    legend(ncol=2)



    title('TF & epigenetic enrichment profiles relative to %s peaks for classes of peaks based on Gene & Feature types' % (refFactorID[0]))
    f1.set_figwidth(25)
    f1.set_figheight(25)

    pp.savefig()
    close()


    f1 = figure(1)
    ax = f1.add_subplot(111)
    numtracks = len(allpeakprofiles_peakclass[0])
    ind = 0
    for trackind in xrange(numtracks):
	if len(allpeakprofiles_peakclass[0][trackind]) > 0:
	    ax.plot(log(allpeakprofiles_peakclass[0][trackind]/allpeakprofiles_peakclass[1][trackind]),lines[ind],label='log(ENHANCER/IHIBITOR) %s' % tracknames[trackind])
#			ax.plot(allpeakprofiles_peakclass[0][trackind]-allpeakprofiles_peakclass[1][trackind],c=cols[ind],linestyle='-',label='ENHANCER - INHIBITOR %s' % tracknames[trackind])

	ind += 1

#	ind = 0
#	numtracks = len(allpeakprofiles_peakclass[1])
#	for trackind in xrange(numtracks):
#		if len(allpeakprofiles_peakclass[1][trackind]) > 0:
#			ax.plot(allpeakprofiles_peakclass[1][trackind],c=cols[ind],linestyle='--',label='INHIBITOR %s' % tracknames[trackind])
#		ind += 1


    N4 = len(allpeakprofiles_peakclass[1][0])
    xaxis_ticks = list(sorted(-array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)
    xaxis_labels = ['%s' % ('' if x == 0 else str(x*span/1000.0) + 'Kb' if abs(x*span) >= 1000 and abs(x*span) < 1000000 else (str(x*span/1000000.0)+'Mb' if abs(x*span) >= 1000000 else str(x*span)+'b')) for x in list(array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])-N4/2) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0]))]

    ymin, ymax = ylim()
    plot([N4/2,N4/2],[ymin,ymax],'k',linewidth=1)
    xticks(xaxis_ticks,xaxis_labels)
    xlim(0,N4-1)

    ax.yaxis.grid(True)
    ax.xaxis.grid(True)

    for tick in ax.get_xticklabels():
	tick.set_fontsize(15)


    xlabel('Offset from %s peaks (bases)' % refFactorID[0])
    ylabel('Average profiles')
    legend(ncol=2)



    title('TF & epigenetic enrichment profiles relative to %s peaks for classes of peaks based on Feature types' % (refFactorID[0]))
    f1.set_figwidth(25)
    f1.set_figheight(25)

    pp.savefig()
    close()



    for geneclass_width_ind in xrange(len(geneclass_widths)):
#		print allpeakprofiles_geneclass[geneclass_width_ind]
	f1 = figure(1)
	ax = f1.add_subplot(111)
	numtracks = len(allpeakprofiles_geneclass[geneclass_width_ind][0])

	#print 'numtracks', numtracks
	ind = 0
	for trackind in xrange(numtracks):

	    if len(allpeakprofiles_geneclass[geneclass_width_ind][0][trackind]) > 0:
		ax.plot(log(allpeakprofiles_geneclass[geneclass_width_ind][0][trackind]/allpeakprofiles_geneclass[geneclass_width_ind][1][trackind]),lines[ind],label='log(ENHANCER/IHIBITOR) %s' % tracknames[trackind])
#				ax.plot(allpeakprofiles_geneclass[geneclass_width_ind][0][trackind]-allpeakprofiles_geneclass[geneclass_width_ind][1][trackind],lines[ind],label='ENHANCER - INHIBITOR %s' % tracknames[trackind])

#				ax.plot(allpeakprofiles_geneclass[geneclass_width_ind][0][trackind],c=cols[ind],linestyle='-',label='ENHANCER - INHIBITOR %s' % tracknames[trackind])

	    ind += 1

#		numtracks = len(allpeakprofiles_geneclass[geneclass_width_ind][1])
#		ind = 0

#		for trackind in xrange(numtracks):
#			if len(allpeakprofiles_geneclass[geneclass_width_ind][1][trackind]) > 0:
#				ax.plot(allpeakprofiles_geneclass[geneclass_width_ind][1][trackind],c=cols[ind],linestyle='--',label='INHIBITOR %s' % tracknames[trackind])
#			ind += 1


	N4 = len(allpeakprofiles_geneclass[geneclass_width_ind][1][0])
	xaxis_ticks = list(sorted(-array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0])+N4/2)
	xaxis_labels = ['%s' % ('' if x == 0 else str(x*span/1000.0) + 'Kb' if abs(x*span) >= 1000 and abs(x*span) < 1000000 else (str(x*span/1000000.0)+'Mb' if abs(x*span) >= 1000000 else str(x*span)+'b')) for x in list(array([x for x in arange(0,N4/2+1,1) if x*span % outputgroupsize == 0])-N4/2) +list(array([x for x in arange(1,N4/2+1,1) if x*span % outputgroupsize == 0]))]

	ymin, ymax = ylim()
	plot([N4/2,N4/2],[ymin,ymax],'k',linewidth=1)
	xticks(xaxis_ticks,xaxis_labels)
	xlim(0,N4-1)

	ax.yaxis.grid(True)
	ax.xaxis.grid(True)

	for tick in ax.get_xticklabels():
	    tick.set_fontsize(15)

	xlabel('Offset from %s peaks (bases)' % refFactorID[0])
	ylabel('Average profiles')
	legend(ncol=2)


	title('TF & epigenetic enrichment profiles relative to %s peaks for classes of peaks based on Gene class in region of %d around TSS' % (refFactorID[0],geneclass_widths[geneclass_width_ind]))
	f1.set_figwidth(25)
	f1.set_figheight(25)

	pp.savefig()
	close()







    pp.close()




def compute_binned_profile_summaries(allpeakprofiles,peaks,refFactorID,bin = 4000,span=10):
    profile_summaries = []


    for profileind,peakprofile in enumerate(allpeakprofiles):
	numtracks = len(peakprofile)
	profilesummary = -ones(numtracks)
	for trackind in xrange(numtracks):

	    if len(peakprofile[trackind]) == 0:
		continue

	    N4 = size(peakprofile[trackind],1)

	    averageprofile = peakprofile[trackind].mean(0)
	    average_profile_bin = sum(averageprofile[N4/2-bin/span:N4/2+bin/span])
	    profilesummary[trackind] = average_profile_bin


	profile_summaries.append([peaks[refFactorID][profileind][0],peaks[refFactorID][profileind][1],profilesummary])

    return profile_summaries






def compute_binned_profile_summaries_twoclasses(allpeakprofiles):
    profile_summaries = [[],[]]

    numtracks = len(allpeakprofiles[0])
    profile_summaries[0] = [[] for i in xrange(numtracks)]
    for trackind in xrange(numtracks):

	if len(allpeakprofiles[0][trackind]) == 0:
	    continue

	N4 = size(allpeakprofiles[0][trackind],1)



	averageprofile = array([x for x in allpeakprofiles[0][trackind] if sum(x) != 0]).mean(0)
	profile_summaries[0][trackind] = averageprofile


    numtracks = len(allpeakprofiles[1])
    profile_summaries[1] = [[] for i in xrange(numtracks)]
    for trackind in xrange(numtracks):

	if len(allpeakprofiles[1][trackind]) == 0:
	    continue

	N4 = size(allpeakprofiles[1][trackind],1)

	averageprofile = array([x for x in allpeakprofiles[1][trackind] if sum(x) != 0]).mean(0)

	profile_summaries[1][trackind] = averageprofile

    return profile_summaries


def create_classification_datasets(allpeakprofiles_geneclasspeakclass,allpeakprofiles_peakclass,allpeakprofiles_geneclass,refFactorID,peaksfn,genefn,trackfns,trackpeakfns,tracknames,peakwidth,window_size,geneclass_widths,use_orientation,pixel_size,outputfn,span=10):




    profiles = allpeakprofiles_geneclasspeakclass
    outputfile = open('%s_%s_geneclasspeakclass_%d_%d_%s.dat' % (outputfn,refFactorID[0],peakwidth,window_size,use_orientation),'w')

    numtracks = len(profiles[0])

    N4 = size(profiles[0][0],1)
    binsize = pixel_size/span

    featuresize = N4/(binsize)


    print >> outputfile, '# windowsize = %d\n# pixelsize = %d\n# numtracks = %d\n# featuresize = %d' % (window_size*2,pixel_size,numtracks,featuresize)
    print >> outputfile, '# genes = %s\n# featurepeaks = %s\n# peaks = %s\n# tracks = %s' % (genefn, peaksfn, trackpeakfns[0],'; '.join(trackfns))


    N = size(profiles[0][0],0)
    N4 = size(profiles[0][0],1)
    for profind in xrange(N):
	print >> outputfile, 1,
	for trackind in xrange(numtracks):
	    pixeldata = [sum(profiles[0][trackind][profind,i:i+binsize]) for i in xrange(0,N4-binsize+1,binsize)]
	    for x in pixeldata:
		print >> outputfile, x,
	print >> outputfile

    N = size(profiles[1][0],0)
    N4 = size(profiles[1][0],1)
    for profind in xrange(N):
	print >> outputfile, -1,
	for trackind in xrange(numtracks):
	    pixeldata = [sum(profiles[1][trackind][profind,i:i+binsize]) for i in xrange(0,N4-binsize+1,binsize)]
	    for x in pixeldata:
		print >> outputfile, x,
	print >> outputfile
    outputfile.close()


    profiles = allpeakprofiles_peakclass
    outputfile = open('%s_%s_peakclass_%d_%d_%s.dat' % (outputfn,refFactorID[0],peakwidth,window_size,use_orientation),'w')

    numtracks = len(profiles[0])

    N4 = size(profiles[0][0],1)
    featuresize = N4/(binsize)


    print >> outputfile, '# windowsize = %d\n# pixelsize = %d\n# numtracks = %d\n# featuresize = %d' % (window_size*2,pixel_size,numtracks,featuresize)
    print >> outputfile, '# genes = %s\n# featurepeaks = %s\n# peaks = %s\n# tracks = %s' % (genefn, peaksfn, trackpeakfns[0],'; '.join(trackfns))

    binsize = pixel_size/span

    N = size(profiles[0][0],0)
    N4 = size(profiles[0][0],1)
    for profind in xrange(N):
	print >> outputfile, 1,
	for trackind in xrange(numtracks):
	    pixeldata = [sum(profiles[0][trackind][profind,i:i+binsize]) for i in xrange(0,N4-binsize+1,binsize)]
	    for x in pixeldata:
		print >> outputfile, x,
	print >> outputfile

    N = size(profiles[1][0],0)
    N4 = size(profiles[1][0],1)
    for profind in xrange(N):
	print >> outputfile, -1,
	for trackind in xrange(numtracks):
	    pixeldata = [sum(profiles[1][trackind][profind,i:i+binsize]) for i in xrange(0,N4-binsize+1,binsize)]
	    for x in pixeldata:
		print >> outputfile, x,
	print >> outputfile
    outputfile.close()


    for geneclasswidth_ind in xrange(len(geneclass_widths)):

	profiles = allpeakprofiles_geneclass[geneclasswidth_ind]
	outputfile = open('%s_%s_geneclass_%d_%d_%s.dat' % (outputfn,refFactorID[0],window_size,geneclass_widths[geneclasswidth_ind],use_orientation),'w')

	numtracks = len(profiles[0])

	N4 = size(profiles[0][0],1)
	featuresize = N4/(binsize)

	print >> outputfile, '# windowsize = %d\n# pixelsize = %d\n# numtracks = %d\n# featuresize = %d\n# geneclass_width = %d' % (window_size*2,pixel_size,numtracks,featuresize,geneclass_widths[geneclasswidth_ind])
	print >> outputfile, '# genes = %s\n# featurepeaks = %s\n# peaks = %s\n# tracks = %s' % (genefn, peaksfn, trackpeakfns[0],'; '.join(trackfns))

	binsize = pixel_size/span

	N = size(profiles[0][0],0)
	N4 = size(profiles[0][0],1)
	for profind in xrange(N):
	    print >> outputfile, 1,
	    for trackind in xrange(numtracks):
		pixeldata = [sum(profiles[0][trackind][profind,i:i+binsize]) for i in xrange(0,N4-binsize+1,binsize)]
		for x in pixeldata:
		    print >> outputfile, x,
	    print >> outputfile

	N = size(profiles[1][0],0)
	N4 = size(profiles[1][0],1)
	for profind in xrange(N):
	    print >> outputfile, -1,
	    for trackind in xrange(numtracks):
		pixeldata = [sum(profiles[1][trackind][profind,i:i+binsize]) for i in xrange(0,N4-binsize+1,binsize)]
		for x in pixeldata:
		    print >> outputfile, x,
	    print >> outputfile
	outputfile.close()

if __name__=='__main__':

    args = sys.argv[1:]

    if len(args) == 0:
	print '=====INPUTS=====\nwindow_size = int(args[0])\npixel_size = args[1]\ndataaggregation = args[2]\nnormalize_data = args[3] == "1" or args[3] == "yes" or args[3] == "on"\noutputfn = args[4]\ngenesfn = args[5]\ntrackfns = args[6:]\n================='
	exit()



    peakwidth = int(args[0])
    window_size = int(args[1])
    geneclass_widths = [int(x) for x in args[2].rsplit(',')]
    use_orientation = args[3] == '1'
    pixel_size = int(args[4])
    refFactorIDs = [(x.rsplit(',')[0],int(x.rsplit(',')[1])) for x in args[5].rsplit(':')]
    outputfn = args[6]
    peaksfn = args[7]
    genesfn = args[8]
    tracknames = args[9].rsplit(',')
    numtracks = int(args[10])
    trackfns = args[11:11+numtracks]
    print trackfns

    trackpeakfns = args[11+numtracks:]
    print trackpeakfns


    print 'Reading Peaks'
    peaks = read_peaks(peaksfn)
    print 'Reading Genes'
    genes = read_genes(genesfn)
    print 'Reading Tracks'
    tracks = read_tracks(trackfns)
    print 'Reading Track Peaks'
    trackpeaks = read_trackpeaks(trackpeakfns)


    for refFactorID in refFactorIDs:
	print refFactorID
	print 'Building profiles'
	allpeakprofiles_geneclasspeakclass,allpeakprofiles_peakclass,allpeakprofiles_geneclass,smallestdistancestopeak = create_agregate_profiles(genes,tracks,trackpeaks,peaks,peakwidth,window_size,refFactorID,use_orientation,geneclass_widths)
	create_classification_datasets(allpeakprofiles_geneclasspeakclass,allpeakprofiles_peakclass,allpeakprofiles_geneclass,refFactorID,peaksfn,genesfn,trackfns,trackpeakfns,tracknames,peakwidth,window_size,geneclass_widths,use_orientation,pixel_size,outputfn)


	allpeakprofiles_geneclasspeakclass = compute_binned_profile_summaries_twoclasses(allpeakprofiles_geneclasspeakclass)
	allpeakprofiles_peakclass = compute_binned_profile_summaries_twoclasses(allpeakprofiles_peakclass)
	for geneclass_width_ind in xrange(len(geneclass_widths)):
	    allpeakprofiles_geneclass[geneclass_width_ind] = compute_binned_profile_summaries_twoclasses(allpeakprofiles_geneclass[geneclass_width_ind])



#	store_peakprofiles(outputfn,allpeakprofiles_geneclasspeakclass,allpeakprofiles_peakclass,allpeakprofiles_geneclass)

#	profile_summaries = compute_binned_profile_summaries(allpeakprofiles,peaks,refFactorID)

#	print 'Plotting profiles'
#	plot_all_peakprofiles(peaks,allpeakprofiles,smallestdistancestopeak,profile_summaries,refFactorID,tracknames,outputfn,window_size)


	plot_all_peakprofiles_twoclasses(allpeakprofiles_geneclasspeakclass,allpeakprofiles_peakclass,allpeakprofiles_geneclass,smallestdistancestopeak,refFactorID,tracknames,outputfn,window_size,geneclass_widths)
#	plot_all_peakprofiles_twoclasses(allpeakprofiles_geneclasspeakclass,smallestdistancestopeak,refFactorID,tracknames,outputfn,window_size)




















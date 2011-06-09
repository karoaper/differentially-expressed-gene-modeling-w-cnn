from numpy import *
from re import *
import sys
from operator import itemgetter
import os

def read_genes(genefn):
    genes = []
    for line in open(genefn,'r'):
	parts = line.rstrip().rsplit()
	genes.append([parts[1],int(parts[2]),parts[3],parts[4]])
    return genes


def read_tracks(trackfns):
    tracks = []

    for trackfn in trackfns:
	if trackfn == '0':
	    tracks.append([])
	else:
	    onetrack = {}

	    for line in open(trackfn,'r').readlines()[1:]:
		parts = line.rstrip().rsplit('\t')
		if onetrack.has_key(parts[0]):
		    onetrack[parts[0]].append([int(parts[1]),int(parts[2]),float(parts[4])])
		else:
		    onetrack[parts[0]] = [[int(parts[1]),int(parts[2]),float(parts[4])]]

	    onetrack_rightsort = {}
	    for key in onetrack.keys():
		onetrack[key].sort()
		onetrack_rightsort[key] = sorted(onetrack[key],key=itemgetter(1),reverse=True)

	    tracks.append([onetrack,onetrack_rightsort])

    return tracks


def fast_forward_peaks_pos(loca,locb,startind,peaks,dataagggregation,overlap_prop):
    loclen = locb-loca

    peakoverlap = 0
    retind = startind
    juststarted = True
    for ind in range(startind,len(peaks)):
	gapc = locb - peaks[ind][1]
	if gapc < 0 and juststarted:
	    retind = ind
	    juststarted = False

	if dataaggregation == 'binary':

	    if peaks[ind][0] <= loca and peaks[ind][1] >= locb:
		return [retind,1]
	    elif peaks[ind][0] >= loca and peaks[ind][1] <= locb:
		return [retind,1]
	    else:
		gapa = locb - peaks[ind][0]
		gapb = peaks[ind][1] - loca

		if gapa < 0:
		    return [retind,0]
		else:
		    gap = min(gapa,gapb)
		    if gap >= 0:
			return [retind,1]

	elif dataaggregation == 'binary_overhalf':

	    if peaks[ind][0] <= loca and peaks[ind][1] >= locb:
		return [retind,1]
	    elif peaks[ind][0] >= loca and peaks[ind][1] <= locb:
		peakoverlap += peaks[ind][1]-peaks[ind][0]
	    else:
		gapa = locb - peaks[ind][0]
		gapb = peaks[ind][1] - loca

		if gapa < 0:
		    return [retind,int(peakoverlap*1.0/loclen >= 0.5)]
		else:
		    gap = min(gapa,gapb)
		    if gap >= 0:
			peakoverlap += gap

	elif dataaggregation == 'average':
	    if peaks[ind][0] <= loca and peaks[ind][1] >= locb:
		if overlap_prop:
		    peakoverlap += loclen*1.0/(peaks[ind][1]-peaks[ind][0])*peaks[ind][2]*loclen
		else:
		    peakoverlap += loclen*peaks[ind][2]
		return [retind,peakoverlap*1.0/loclen]
	    elif peaks[ind][0] >= loca and peaks[ind][1] <= locb:
		peakoverlap += (peaks[ind][1]-peaks[ind][0])*peaks[ind][2]
	    else:
		gapa = locb - peaks[ind][0]
		gapb = peaks[ind][1] - loca

		if gapa < 0:
		    return [retind,peakoverlap*1.0/loclen]
		else:
		    gap = min(gapa,gapb)
		    if gap >= 0:
			if overlap_prop:
			    peakoverlap += gap*1.0/(peaks[ind][1]-peaks[ind][0])*peaks[ind][2]*gap
			else:
			    peakoverlap += peaks[ind][2]*gap

	elif dataaggregation == 'sum':
	    if peaks[ind][0] <= loca and peaks[ind][1] >= locb:
		if overlap_prop:
		    peakoverlap += loclen*1.0/(peaks[ind][1]-peaks[ind][0])*peaks[ind][2]
		else:
		    peakoverlap += peaks[ind][2]
		return [retind,peakoverlap]
	    elif peaks[ind][0] >= loca and peaks[ind][1] <= locb:
		peakoverlap += peaks[ind][2]
	    else:
		gapa = locb - peaks[ind][0]
		gapb = peaks[ind][1] - loca

		if gapa < 0:
		    return [retind,peakoverlap]
		else:
		    gap = min(gapa,gapb)
		    if gap >= 0:
			if overlap_prop:
			    peakoverlap += gap*1.0/(peaks[ind][1]-peaks[ind][0])*peaks[ind][2]
			else:
			    peakoverlap += peaks[ind][2]

    return [len(peaks),0]


def fast_forward_peaks_neg(loca,locb,startind,peaks,dataagggregation,overlap_prop):
    loclen = locb-loca

    peakoverlap = 0
    retind = startind
    juststarted = True

    for ind in range(startind,len(peaks)):
	gapc = loca - peaks[ind][0]
	if gapc > 0 and juststarted:
	    retind = ind
	    juststarted = False

	if dataaggregation == 'binary':

	    if peaks[ind][0] <= loca and peaks[ind][1] >= locb:
		return [retind,1]
	    elif peaks[ind][0] >= loca and peaks[ind][1] <= locb:
		return [retind,1]
	    else:
		gapa = locb - peaks[ind][0]
		gapb = peaks[ind][1] - loca

		if gapb < 0:
		    return [retind,0]
		else:
		    gap = min(gapa,gapb)
		    if gap >= 0:
			return [retind,1]

	elif dataaggregation == 'binary_overhalf':

	    if peaks[ind][0] <= loca and peaks[ind][1] >= locb:
		return [retind,1]
	    elif peaks[ind][0] >= loca and peaks[ind][1] <= locb:
		peakoverlap += peaks[ind][1]-peaks[ind][0]
	    else:
		gapa = locb - peaks[ind][0]
		gapb = peaks[ind][1] - loca

		if gapb < 0:
		    return [retind,int(peakoverlap*1.0/loclen >= 0.5)]
		else:
		    gap = min(gapa,gapb)
		    if gap >= 0:
			peakoverlap += gap

	elif dataaggregation == 'average':
	    if peaks[ind][0] <= loca and peaks[ind][1] >= locb:
		if overlap_prop:
		    peakoverlap += loclen*1.0/(peaks[ind][1]-peaks[ind][0])*peaks[ind][2]*loclen
		else:
		    peakoverlap += loclen*peaks[ind][2]
		return [retind,peakoverlap*1.0/loclen]
	    elif peaks[ind][0] >= loca and peaks[ind][1] <= locb:
		peakoverlap += (peaks[ind][1]-peaks[ind][0])*peaks[ind][2]
	    else:
		gapa = locb - peaks[ind][0]
		gapb = peaks[ind][1] - loca

		if gapb < 0:
		    return [retind,peakoverlap*1.0/loclen]
		else:
		    gap = min(gapa,gapb)
		    if gap >= 0:
			if overlap_prop:
			    peakoverlap += gap*1.0/(peaks[ind][1]-peaks[ind][0])*peaks[ind][2]*gap
			else:
			    peakoverlap += peaks[ind][2]*gap

	elif dataaggregation == 'sum':
	    if peaks[ind][0] <= loca and peaks[ind][1] >= locb:
		if overlap_prop:
		    peakoverlap += loclen*1.0/(peaks[ind][1]-peaks[ind][0])*peaks[ind][2]
		else:
		    peakoverlap += peaks[ind][2]
		return [retind,peakoverlap]
	    elif peaks[ind][0] >= loca and peaks[ind][1] <= locb:
		peakoverlap += peaks[ind][2]
	    else:
		gapa = locb - peaks[ind][0]
		gapb = peaks[ind][1] - loca

		if gapb < 0:
		    return [retind,peakoverlap]
		else:
		    gap = min(gapa,gapb)
		    if gap >= 0:
			if overlap_prop:
			    peakoverlap += gap*1.0/(peaks[ind][1]-peaks[ind][0])*peaks[ind][2]
			else:
			    peakoverlap += peaks[ind][2]

    return [len(peaks),0]


def create_track_pixel_matrixes(genes,tracks,winsize,pixelsize,dataagggregation,overlap_prop,normalize_data):

    numgenes = len(genes)
    numtracks = len(tracks)

    feature_length = len(range(-winsize/2,winsize/2,pixelsize))

    print 'Using pixel_size = %d and winsize = %d, resulting in feature_size = %d (numtracks = %d)' % (pixelsize,winsize,feature_length,numtracks)
    pixel_matrixes = [zeros([numgenes,feature_length]) for i in range(numtracks)]

    for geneind,gene in enumerate(genes):
	if geneind*1.0/len(genes) * 100 % 25 < 1 and (geneind+1)*1.0/len(genes) * 100 % 25 > 1:
	    print 'Completed %f%% out of total %d genes' % (floor(geneind*1.0/len(genes) * 100),len(genes))

	tss = gene[1]
	chrom = gene[0]
	orient = gene[2]
	if orient == '+':
	    winstart = tss - winsize/2
	    winend = tss + winsize/2
	else:
	    winstart = tss + winsize/2
	    winend = tss - winsize/2

	for trackind in range(numtracks):
	    if tracks[trackind] == []:
		print 'adding blanks in track %d' % trackind
		pixel_matrixes[trackind][geneind,:] = 0

	    else:

		if not tracks[trackind][0].has_key(chrom):
		    print 'no chrom %s in track %d' % (chrom,trackind)
		    pixel_matrixes[trackind][geneind,:] = 0
		    continue

		newind = 0
		if gene[2] == '+':
		    peaks = tracks[trackind][0][chrom]
		    for pixind,[loca,locb] in enumerate([[a,a+pixelsize] for a in range(winstart,winend,pixelsize)]):
			[newind,pixval] = fast_forward_peaks_pos(loca,locb,newind,peaks,dataagggregation,overlap_prop)
			pixel_matrixes[trackind][geneind,pixind] = pixval
		else:
		    peaks = tracks[trackind][1][chrom]
		    for pixind,[loca,locb] in enumerate([[a-pixelsize,a] for a in range(winstart,winend,-pixelsize)]):
			[newind,pixval] = fast_forward_peaks_neg(loca,locb,newind,peaks,dataagggregation,overlap_prop)
			pixel_matrixes[trackind][geneind,pixind] = pixval


    if not dataaggregation == 'binary' and not dataaggregation == 'binary_overhalf' and normalize_data:
	for trackind in range(numtracks):
	    maxpeak = pixel_matrixes[trackind].max()
	    pixel_matrixes[trackind] = pixel_matrixes[trackind]*1.0/maxpeak if maxpeak != 0 else 0


    return pixel_matrixes



def get_pixel_size(genes,tracks):


    min_peak_len = infty
    min_gap_len = infty

    allpeaks = []

    print 'Min/Average peak and gap lengths per track are:'
    for trackind,track in enumerate(tracks):
	gap_lens_track = []
	peak_lens_track = []
	for chrom in track[0].keys():
	    chromdata = track[0][chrom]
	    gap_lens_track.extend([chromdata[i+1][0]-chromdata[i][1] for i in range(len(chromdata)-1) if chromdata[i+1][0]-chromdata[i][1] > 0])
	    peak_lens_track.extend([chromdata[i][1]-chromdata[i][0] for i in range(len(chromdata))])

	min_gap_len_track = min(gap_lens_track)
	min_peak_len_track = min(peak_lens_track)
	aver_gap_len_track = mean(gap_lens_track)
	aver_peak_len_track = mean(peak_lens_track)
	print 'track %d: %d/%f, %d/%f' % (trackind,min_peak_len_track,aver_peak_len_track,min_gap_len_track,aver_gap_len_track)
	min_peak_len = min(min_peak_len,min_peak_len_track)
	min_gap_len = min(min_gap_len,min_gap_len_track)


    return [min_gap_len,min_peak_len]




def store_feature_matrixes(pixel_matrixes,genes,feature_fn, pixel_size, window_size, genefn, tracksfn,trackstobuild):

    featurefile = open(feature_fn,'w')
    print >> featurefile, '# windowsize = %d\n# pixelsize = %d\n# numtracks = %d\n# featuresize = %d' % (window_size,pixel_size,len(pixel_matrixes),size(pixel_matrixes[0],1))
    print >> featurefile, '# genes = %s\n# tracks = %s' % (genefn, '; '.join([tracksfn[trackcombo[0]] + '-' + tracksfn[trackcombo[1]] if len(trackcombo)>1 else tracksfn[trackcombo[0]] for trackcombo in trackstobuild]))
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



#def compute_track_diff(track1,track2):
#   numgenes = size(track1,0)
#   feature_length = size(track1,1)
#
#   diff_track = zeros([numgenes,feature_length])
#
#   for geneind in xrage(numgenes):
#   diff_track = track1 - track2


def build_combined_matrixes(pixel_matrixes,trackstobuild):

    numgenes = size(pixel_matrixes[0],0)
    numtracks = len(trackstobuild)

    feature_length = size(pixel_matrixes[0],1)

    new_pixel_matrixes = []#zeros([numgenes,feature_length]) for i in range(numtracks)]

    for ind,trackcombo in enumerate(trackstobuild):
	if len(trackcombo) == 1:
	    new_pixel_matrixes.append(pixel_matrixes[trackcombo[0]])
	else:
#	    new_pixel_matrixes.append(compute_track_diff(pixel_matrixes[trackcombo[0]],pixel_matrixes[trackcombo[1]]))
	    new_pixel_matrixes.append(pixel_matrixes[trackcombo[0]]- pixel_matrixes[trackcombo[1]])


    return new_pixel_matrixes


if __name__=='__main__':

    args = sys.argv[1:]

    if len(args) == 0:
	print '=====INPUTS=====\nwindow_size = int(args[0])\npixel_size = args[1]\ndataaggregation = args[2]\noverlap_prop = args[3] == "1" or args[3] == "yes" or args[3] == "on"\nnormalize_data = args[4] == "1" or args[4] == "yes" or args[4] == "on"\noutputfn = args[5]\ngenesfn = args[6]\ntrackstobuild = [[int(x) for x in y.rsplit(",")] for y in args[7].rsplit(":")]\ntrackfns = args[8:]\n================='
	exit()



    window_size = int(args[0])
    pixel_size = args[1]
    dataaggregation = args[2]
    overlap_prop = args[3] == '1' or args[3] == 'yes' or args[3] == 'on'
    normalize_data = args[4] == '1' or args[4] == 'yes' or args[4] == 'on'
    outputfn = args[5]
    genesfn = args[6]
    trackstobuild = [[int(x) for x in y.rsplit(',')] for y in args[7].rsplit(':')]
    trackfns = args[8:]


    print 'Reading Genes'
    genes = read_genes(genesfn)
    print 'Reading Tracks'
    tracks = read_tracks(trackfns)
    print 'Building pixel matrixes'

    if pixel_size == 'mingap':
	pixel_size = get_pixel_size(genes,tracks)[0]
    elif pixel_size == 'minpeak':
	pixel_size = get_pixel_size(genes,tracks)[1]
    elif pixel_size == 'minpeakgap':
	pixel_size = min(get_pixel_size(genes,tracks))
    else:
	pixel_size = int(pixel_size)

    pixel_matrixes = create_track_pixel_matrixes(genes,tracks,window_size,pixel_size,dataaggregation,overlap_prop,normalize_data)
    pixel_matrixes = build_combined_matrixes(pixel_matrixes,trackstobuild)

    print 'Storing pixel matrixes'

    folder = '/'.join(outputfn.rsplit('/')[:-1])

    if not (folder == '' or os.path.exists(folder)):
	os.mkdir(folder)

    store_feature_matrixes(pixel_matrixes,genes,outputfn,pixel_size,window_size,genesfn,trackfns,trackstobuild)






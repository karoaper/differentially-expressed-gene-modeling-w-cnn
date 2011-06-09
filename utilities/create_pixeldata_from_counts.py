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
		    parts = line.rstrip().replace('=',' ').rsplit(' ')
		    chrom = parts[0]
		    print 'Reading chromosom %s' % chrom

		    onetrack[chrom] = zeros(int(parts[6]),dtype='int16')
		    spans[chrom] = int(parts[2])
		    starts[chrom] = int(parts[4])
		    ind = 0
		    continue

		onetrack[chrom][ind] = int(line.rstrip())
		ind += 1

	    tracks.append([onetrack,spans,starts])

    return tracks





def return_count_agg(loca,locb,counts,dataagggregation):
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



def create_track_pixel_matrixes(genes,tracks,winsize,pixelsize,dataagggregation,normalize_data):

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


		counts = tracks[trackind][0][chrom]
		span = tracks[trackind][1][chrom]
		start = tracks[trackind][2][chrom]

		newind = 0
		if gene[2] == '+':
		    for pixind,[loca,locb] in enumerate([[a,a+pixelsize] for a in range(winstart,winend,pixelsize)]):
			pixval = return_count_agg((loca-start)/span,(locb-start)/span,counts,dataagggregation)
			pixel_matrixes[trackind][geneind,pixind] = pixval
		else:
		    peaks = tracks[trackind][1][chrom]
		    for pixind,[loca,locb] in enumerate([[a-pixelsize,a] for a in range(winstart,winend,-pixelsize)]):
			pixval = return_count_agg((loca-start)/span,(locb-start)/span,counts,dataagggregation)
			pixel_matrixes[trackind][geneind,pixind] = pixval



    pixels_to_skip = []
    if normalize_data:
	for trackind in range(numtracks):
	    maxval = pixel_matrixes[trackind].max()
	    maxvals = pixel_matrixes[trackind].max(0)
	    pixels_to_skip.append([i for i in xrange(size(maxvals)) if maxvals[i] == 0])
	    pixels_to_keep = [i for i in xrange(size(maxvals)) if maxvals[i] != 0]

	    pixel_matrixes[trackind] = pixel_matrixes[trackind]*1.0/maxval if maxval != 0 else 0
	    pixel_matrixes[trackind] = pixel_matrixes[trackind][:,pixels_to_keep]



    for trackind in range(numtracks):
	print 'Number of void pixels in track %d = %d' % (trackind,len(pixels_to_skip[trackind]))

    return pixel_matrixes,pixels_to_skip



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




if __name__=='__main__':




    args = sys.argv[1:]

    if len(args) == 0:
	print '=====INPUTS=====\nwindow_size = int(args[0])\npixel_size = args[1]\ndataaggregation = args[2]\nnormalize_data = args[3] == "1" or args[3] == "yes" or args[3] == "on"\noutputfn = args[4]\ngenesfn = args[5]\ntrackfns = args[6:]\n================='
	exit()



    window_size = int(args[0])
    pixel_size = args[1]
    dataaggregation = args[2]
    normalize_data = args[3] == '1' or args[3] == 'yes' or args[3] == 'on'
    outputfn = args[4]
    genesfn = args[5]
    trackfns = args[6:]

    print 'Reading Genes'
    genes = read_genes(genesfn)
    print 'Reading Tracks'
    tracks = read_tracks(trackfns)
    print 'Building pixel matrixes'

    pixel_size = int(pixel_size)

    pixel_matrixes,skipped_pixels = create_track_pixel_matrixes(genes,tracks,window_size,pixel_size,dataaggregation,normalize_data)
    print 'Storing pixel matrixes'

    folder = '/'.join(outputfn.rsplit('/')[:-1])

    if not (folder == '' or os.path.exists(folder)):
	os.mkdir(folder)

    store_feature_matrixes(pixel_matrixes,skipped_pixels,genes,outputfn,pixel_size,window_size,genesfn,trackfns)










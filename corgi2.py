#!/usr/bin/env python
import sys
import re
import os
import json
import argparse
import numpy as np

# absolute path to this script
SIM_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
sys.path.append(SIM_PATH+'/py/')

from mappability_corgi import MappabilityTrack
from misc_corgi import *

# parse input args
parser = argparse.ArgumentParser(description='CORGi v2.0a')
parser.add_argument('-i',  type=str, required=False, metavar='<str>',                     help="input_readname_dict.json")
parser.add_argument('-of', type=str, required=True,  metavar='<str>',                     help="* output_clipped_sequences.fa")
parser.add_argument('-ov', type=str, required=True,  metavar='<str>',                     help="* output_svs.dat")
#
parser.add_argument('--min-mapq', type=int, required=False,  metavar='<int>', default=3,  help="minimum MAPQ")
parser.add_argument('--min-clip', type=int, required=False,  metavar='<int>', default=50, help="minimum softclip size (S)")
parser.add_argument('--min-ins',  type=int, required=False,  metavar='<int>', default=50, help="minimum ins size (I)")
parser.add_argument('--min-del',  type=int, required=False,  metavar='<int>', default=50, help="minimum del size (D)")
#
args = parser.parse_args()

INPUT_READ_DICT = args.i
OUTPUT_READS    = args.of		# interesting_read_segments.fa
OUTPUT_SV_DATA  = args.ov		# sv_json.dat

MIN_MAPQ   = args.min_mapq
MIN_SIZE_S = args.min_clip
MIN_SIZE_I = args.min_ins
MIN_SIZE_D = args.min_del

# bed dir
EXCL_BED = [MappabilityTrack(SIM_PATH + 'bed/hg38_centromere.bed.gz',      bed_buffer=10000),
            MappabilityTrack(SIM_PATH + 'bed/hg38_gap.bed.gz',             bed_buffer=1000),
            MappabilityTrack(SIM_PATH + 'bed/hg38_simpleRepeats.bed.gz',   bed_buffer=50),
            MappabilityTrack(SIM_PATH + 'bed/hg38_microsatellites.bed.gz', bed_buffer=10)]

# read in input readname:data dict
READNAME_DICT = {}
if INPUT_READ_DICT != None:
	f = open(INPUT_READ_DICT,'r')
	READNAME_DICT = json.load(f)
	f.close()
# here's the dict for next round:
NEXT_READNAME_DICT = {}

# various params
allow_hard = True
minXmatch  = 50

CLIP_FULLALN_IDENTITY_THRESH = 0.90

if not sys.stdin.isatty():
	input_stream = sys.stdin
else:
	print 'No input.'
	exit(1)

# cigar ops of interest
CLIP_CHAR  = 'S'
REF_CHAR   = 'MX=D'
READ_CHAR  = 'MX=I'
MATCH_CHAR = 'M='
XMTCH_CHAR = 'XID'
if allow_hard:
	CLIP_CHAR += 'H'

last_printed  = (None,0)
print_every   = 1000000
reads_written = 0

nSkip_done   = 0
nSkip_boring = 0
nSkip_unmap  = 0
nSkip_suppl  = 0
nSkip_bed    = 0
nLowQual     = 0

leaf_svs = {}

OUTPUT_READS_FILE_HANDLE = open(OUTPUT_READS,'w')

for line in input_stream:
	if len(line) and line[0] != '#':
		splt  = line.strip().split('\t')
		cigar = splt[5]
		flag  = int(splt[1])
		ref   = splt[2]
		pos   = int(splt[3])
		rdat  = splt[9]

		rnm   = splt[0]
		if rnm in READNAME_DICT:
			rnm = READNAME_DICT[rnm]

		rname_splt = [n.split('_')[-5:] for n in rnm.split('&')[1:]]
		rname_orig = rnm.split('&')[0]
		depth      = rnm.count('&')
		#dchar      = '_'*(depth == 0) + '&'*(depth >= 1)
		dchar      = '&'
		if depth == 0:
			rnm = rnm.replace('/','-')		# downstream tools don't like backslashes in read names, lets replace with '-'
			rnm = rnm.replace('_','-')		# reserve underscores for our tool-specific delimiters
			rnm = rnm.replace('&','-')		# reserve amprisands too

		#
		myAnchors = [n[0] for n in rname_splt]
		myType    = [n[1] for n in rname_splt]

		orientation = 'FWD'
		if flag&16:
			orientation = 'REV'

		#
		#
		#
		if flag&2048:
			nSkip_suppl += 1
			continue	# skip supplementary alignments (hmm?)

		if 'X' in myType:
			continue	# skip Xs for now. they're so rare and I'm not sure it's worth including them

		if ref in CHR_BLACKLIST:
			continue
		#
		#
		#

		if cigar == '*':
			nSkip_unmap += 1
			if depth == 0:
				continue	# skip unmapped reads at depth 0
			else:
				# output svs that are anchored but connect to unknown sequence
				#print splt
				if '*' not in leaf_svs:
					leaf_svs['*'] = []
				leaf_svs['*'].append(get_sv_traceback(myAnchors, myType, UNMAPPED_ANCHOR, rname_orig))
				continue

		mapq = int(splt[4])
		if mapq < MIN_MAPQ:
			# treat low mapq at depth > 0 as if they were unmapped
			if depth > 0:
				if '*' not in leaf_svs:
					leaf_svs['*'] = []
				my_unmap_anchor = UNMAPPED_ANCHOR
				if myType[-1] == 'I':
					my_unmap_anchor = '*:1-'+str(len(rdat))+':FWD'
				leaf_svs['*'].append(get_sv_traceback(myAnchors, myType, my_unmap_anchor, rname_orig))
			nLowQual += 1
			continue

		####if not(ref in CHR_WHITELIST):
		####	continue
		####if ref in BREAK_CHR:
		####	break

		letters = re.split(r"\d+",cigar)[1:]
		numbers = [int(n) for n in re.findall(r"\d+",cigar)]
		if len(letters) == 0 or len(numbers) == 0:
			print 'invalid line?'
			print splt
			continue

		events = []
		adj    = 0
		radj   = 0
		for i in xrange(len(letters)):
			if letters[i] == 'D' and numbers[i] >= MIN_SIZE_D:
				pass
			elif letters[i] == 'I' and numbers[i] >= MIN_SIZE_I:
				events.append((adj, 'I', radj, numbers[i]))
			elif letters[i] == 'X' and numbers[i] >= minXmatch:
				events.append((adj, 'X', radj, numbers[i]))

			if letters[i] in REF_CHAR:
				adj += numbers[i]
			if letters[i] in READ_CHAR:
				radj += numbers[i]

		if letters[0] in CLIP_CHAR and numbers[0] >= MIN_SIZE_S:
			events.append((0, 'Sl', 0, numbers[0]))
		else:
			events.append((0, 'start', 0, 0))

		if len(letters) > 1 and letters[-1] in CLIP_CHAR and numbers[-1] >= MIN_SIZE_S:
			adj  = 0
			radj = 0
			for i in xrange(len(letters)-1):
				if letters[i] in REF_CHAR:
					adj += numbers[i]
				if letters[i] in READ_CHAR:
					radj += numbers[i]
			events.append((adj, 'Sr', radj, numbers[-1]))
		else:
			events.append((adj, 'end', radj, 0))

		# read name format:
		#
		# 0          1         2           3           4                         5
		# [readname]_[anchors]_[clip_type]_[clip_size]_[clip_starting_pos (ref)]_[clip_starting_pos (read)]

		# define "anchor" sequence length as the distance between noteworthy events, presumably comprised mostly of matches
		out_dat = []
		events  = sorted(events)
		if events[0][1] == 'Sl':
			my_anchor_len = events[1][0] - events[0][0]
			a_r = ref + ':' + str(pos+events[0][0]) + '-' + str(pos+events[0][0] + my_anchor_len) + ':' + orientation
			myRead = rdat[events[0][2]:events[0][2]+events[0][3]]
			if orientation == 'REV':
				myRead = RC(myRead)
			out_dat.append((rnm, a_r, events[0][1], events[0][3], events[0][0]+pos, events[0][2], myRead))

		if events[-1][1] == 'Sr':
			my_anchor_len = events[-1][0] - events[-2][0]
			a_l = ref + ':' + str(pos+events[-1][0] - my_anchor_len) + '-' + str(pos+events[-1][0]) + ':' + orientation
			myRead = rdat[events[-1][2]:events[-1][2]+events[-1][3]]
			if orientation == 'REV':
				myRead = RC(myRead)
			out_dat.append((rnm, a_l, events[-1][1], events[-1][3], events[-1][0]+pos, events[-1][2], myRead))

		for i in xrange(len(events)):
			if events[i][1] == 'I' or events[i][1] == 'X':
				my_anchor_left  = events[i][0] - events[i-1][0]
				my_anchor_right = events[i+1][0] - events[i][0]
				a_l = ref + ':' + str(pos+events[i][0] - my_anchor_left) + '-' + str(pos+events[i][0]) + ':' + orientation
				a_r = ref + ':' + str(pos+events[i][0] + 1) + '-' + str(pos+events[i][0] + my_anchor_right) + ':' + orientation
				myRead = rdat[events[i][2]:events[i][2]+events[i][3]]
				if orientation == 'REV':
					myRead = RC(myRead)
				out_dat.append((rnm, a_l+'+'+a_r, events[i][1], events[i][3], events[i][0]+pos, events[i][2], myRead))

		#
		# no more events found in read, interpret event
		#
		if not len(out_dat):
			if depth == 0:
				nSkip_boring += 1
				continue
			else:
				my_anchor_len = events[1][0] - events[0][0]

				if isValidCoord(EXCL_BED, ref, pos) and isValidCoord(EXCL_BED, ref, pos+my_anchor_len):

					a_i = ref + ':' + str(pos) + '-' + str(pos+my_anchor_len) + ':' + orientation
	
					(sl, sr, mm, xm) = (0, 0, 0, 0)
					if letters[0] == 'S': sl = numbers[0]
					if len(letters) > 1 and letters[-1] == 'S': sr = numbers[-1]
					for i in xrange(len(letters)):
						if letters[i] in MATCH_CHAR: mm += numbers[i]
						if letters[i] in XMTCH_CHAR: xm += numbers[i]
					myIdentity = float(mm)/float(sl+sr+mm+xm)

					if myIdentity >= CLIP_FULLALN_IDENTITY_THRESH:
						if ref not in leaf_svs:
							leaf_svs[ref] = []
						leaf_svs[ref].append(get_sv_traceback(myAnchors, myType, a_i, rname_orig))

					nSkip_done += 1
	
		# write output reads
		anySkip = False
		for n in out_dat:
			if isValidCoord(EXCL_BED, ref, n[4]):
				reads_written += 1
				#OUTPUT_READS_FILE_HANDLE.write('>'+n[0]+dchar+'_'.join([str(m) for m in n[1:6]])+'\n')
				OUTPUT_READS_FILE_HANDLE.write('>'+rname_orig+'_'+str(reads_written)+'\n')
				OUTPUT_READS_FILE_HANDLE.write(n[6]+'\n')

				NEXT_READNAME_DICT[rname_orig+'_'+str(reads_written)] = n[0]+dchar+'_'.join([str(m) for m in n[1:6]])

				if n[3] != len(n[6]):
					print n[3], '!=', len(n[6])
					print n
			else:
				anySkip = True
		if anySkip:
			nSkip_bed += 1

		# print progress
		if ref != last_printed[0]:
			last_printed = (ref,pos)
			print last_printed[0], ':', last_printed[1]
		elif pos - last_printed[1] >= print_every:
			last_printed = (ref,pos)
			print last_printed[0], ':', last_printed[1]

OUTPUT_READS_FILE_HANDLE.close()

# write out read data dict
f = open(OUTPUT_READS+'.readDat.json','w')
json.dump(NEXT_READNAME_DICT,f,sort_keys=True,indent=4)
f.close()

print ''
print 'reads skipped:'
print nSkip_boring+nSkip_done, 'no clipping/ins'
print nSkip_unmap, 'unmapped'
print nSkip_suppl, 'supplementary'
print nSkip_bed, 'bed-exclude regions'
print nLowQual, 'low mapq (treated as unmapped)'
print ''

#
#
#
#
#

if depth == 0:
	print 'done (depth 0).'
	exit(0)

MAX_SV_DIST    = 100
MIN_READ_COUNT = 1

TDUP_BREAKPOINT_TOLERANCE = 10

fo2 = open(OUTPUT_SV_DATA,'w')

for k_chr in sorted(leaf_svs.keys()):
	#if k_chr == '*':
	#	continue

	N = len(leaf_svs[k_chr])
	my_leaf_svs = leaf_svs[k_chr]
	print 'computing distance between', N, 'SVs (' + k_chr + ')...'
	adj_matrix = np.zeros((N,N),dtype='uint8')
	for i in xrange(N):
		#print i
		for j in xrange(i+1,N):
			if len(my_leaf_svs[i]) == len(my_leaf_svs[j]):
				d = [event_dist(my_leaf_svs[i][k], my_leaf_svs[j][k]) for k in xrange(len(my_leaf_svs[i]))]
				allBelowThresh = True
				for n in d:
					if n == None or n > MAX_SV_DIST:
						allBelowThresh = False
				if allBelowThresh:
					adj_matrix[i][j] = 1
					adj_matrix[j][i] = 1

	print 'clustering...'
	clusters = get_connected_subgraphs(adj_matrix)

	for i in xrange(len(clusters)):
		if len(clusters[i]) >= MIN_READ_COUNT:
		
			event_type = '+'.join([n[0] for n in my_leaf_svs[clusters[i][0]]])

			final_annotation = 'NO_ANNOTATION'

			if event_type == 'I':
				if my_leaf_svs[clusters[i][0]][0][1][0][0] == '*' or my_leaf_svs[clusters[i][0]][0][1][1][0] == '*':
					final_annotation = 'NOVEL_INS'
				else:
					bp_i_from = np.median([my_leaf_svs[j][0][1][0][2] for j in clusters[i]])
					bp_i_to   = np.median([my_leaf_svs[j][0][1][1][1] for j in clusters[i]])
					if my_leaf_svs[clusters[i][0]][0][1][0][0] == my_leaf_svs[clusters[i][0]][0][1][1][0]:
						if abs(bp_i_from - bp_i_to) <= TDUP_BREAKPOINT_TOLERANCE:
							final_annotation = 'TDUP'
					else:
						final_annotation = 'DDUP'

			elif event_type == 'Sr' or event_type == 'Sl':
				if my_leaf_svs[clusters[i][0]][0][1][0][0] != my_leaf_svs[clusters[i][0]][0][1][1][0]:
					final_annotation = 'TRA'

			# output format
			#
			# ###  total_event_type    event_depth                          num_reads                             annotation
			# ##   subevent_type       chr : left_anchor - median_bp_pos    chr : median_bp_pos - right_anchor    chr : inserted_sequence - inserted_sequence
			# ##   ...
			# #    readname_1
			# #    ...

			outStr = [['###', event_type, str(event_type.count('+')+1), str(len(clusters[i])), final_annotation]]

			# collapse all reads in cluster into consensus breakpoints
			for j2 in xrange(len(my_leaf_svs[clusters[i][0]])):
				subevents = [my_leaf_svs[j][j2][0] for j in clusters[i]]
				chr_1 = [my_leaf_svs[j][j2][1][0][0] for j in clusters[i]][0]
				chr_2 = [my_leaf_svs[j][j2][1][1][0] for j in clusters[i]][0]
				orr_1 = [my_leaf_svs[j][j2][1][0][3] for j in clusters[i]]
				orr_2 = [my_leaf_svs[j][j2][1][1][3] for j in clusters[i]]
				al_1  = [my_leaf_svs[j][j2][1][0][1] for j in clusters[i]]
				bp_1  = [my_leaf_svs[j][j2][1][0][2] for j in clusters[i]]		# from
				bp_2  = [my_leaf_svs[j][j2][1][1][1] for j in clusters[i]]		# to (start)
				al_2  = [my_leaf_svs[j][j2][1][1][2] for j in clusters[i]]

				if subevents[0][0] == 'S':

					#print chr_1, chr_2, subevents

					if chr_1 == '*' or chr_2 == '*':
						orr_state = 'FWD'
					else:
						if all([orr_1[n][:3] == orr_2[n][:3] for n in xrange(len(orr_1))]):
							orr_state = 'FWD'
						elif all([orr_1[n][:3] != orr_2[n][:3] for n in xrange(len(orr_1))]):
							orr_state = 'REV'
						else:
							print 'weird cluster'
							exit(1)

					# anchor length tends to depend on clipping direction, lets choose the longest on either side
					if 'Sr' in subevents:
						al1_sr = int(np.median([al_1[n] for n in xrange(len(al_1)) if subevents[n] == 'Sr']))
						al2_sr = int(np.median([al_2[n] for n in xrange(len(al_2)) if subevents[n] == 'Sr']))
					else:
						al1_sr = LARGE_NUMBER
						al2_sr = 0
					if 'Sl' in subevents:
						al1_sl = int(np.median([al_1[n] for n in xrange(len(al_1)) if subevents[n] == 'Sl']))
						al2_sl = int(np.median([al_2[n] for n in xrange(len(al_2)) if subevents[n] == 'Sl']))
					else:
						al1_sl = LARGE_NUMBER
						al2_sl = 0

					a_1 = chr_1 + ':' + str(min([al1_sr,al1_sl])) + '-' + str(int(np.median(bp_1))) + ':' + 'FWD'
					a_2 = chr_2 + ':' + str(int(np.median(bp_2))) + '-' + str(max([al2_sr,al2_sl])) + ':' + orr_state

					outStr.append(['##', 'S', a_1, a_2])
				
				if subevents[0] == 'I':
					chr_3 = [my_leaf_svs[j][j2][1][2][0] for j in clusters[i]][0]
					orr_3 = [my_leaf_svs[j][j2][1][2][3] for j in clusters[i]]
					bp_3  = [my_leaf_svs[j][j2][1][2][1] for j in clusters[i]]
					al_3  = [my_leaf_svs[j][j2][1][2][2] for j in clusters[i]]

					if chr_2 == '*':
						orr_state = 'FWD'
					else:
						if all([orr_1[n][:3] == orr_2[n][:3] for n in xrange(len(orr_1))]):
							orr_state = 'FWD'
						elif all([orr_1[n][:3] != orr_2[n][:3] for n in xrange(len(orr_1))]):
							orr_state = 'REV'
						else:
							print 'weird cluster'
							exit(1)

					a_1 = chr_1 + ':' + str(int(np.median(al_1))) + '-' + str(int(np.median(bp_1))) + ':' + 'FWD'
					a_i = chr_2 + ':' + str(int(np.median(bp_2))) + '-' + str(int(np.median(al_2))) + ':' + orr_state
					a_2 = chr_3 + ':' + str(int(np.median(bp_3))) + '-' + str(int(np.median(al_3))) + ':' + 'FWD'

					outStr.append(['##', 'I', a_1, a_2, a_i])

			#print i, clusters[i]
			#for j in clusters[i]:
			#	print '---', my_leaf_svs[j]

			# get readnames
			clust_rn = set()
			[clust_rn.add(my_leaf_svs[j][0][2]) for j in clusters[i]]
			for n in sorted(list(clust_rn)):
				outStr.append(['#',n])

			# redundant conversions between data structures because I can't make up my mind on output format
			meta_dat = ''
			bp_dat   = []
			rn_dat   = '"reads": "' + ', '.join(sorted(list(clust_rn))) + '"'
			for n in outStr:
				if n[0] == '###':
					meta_dat = '"metadata": "' + ', '.join(n[1:]) + '"'
				if n[0] == '##':
					bp_dat.append('"event '+str(len(bp_dat)+1)+'": "' + ', '.join(n[1:]) + '"')

			json_out = '{' + meta_dat + ', ' + ', '.join(bp_dat) + ', ' + rn_dat + '}'
			####print json_out
			fo2.write(json_out+'\n')

fo2.close()

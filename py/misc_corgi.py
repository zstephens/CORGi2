import os
import gzip
import numpy as np

#from scipy.sparse import csr_matrix
#from scipy.sparse.csgraph import reverse_cuthill_mckee

from mappability_corgi import MappabilityTrack
from virus_corgi import getViralName

# absolute path to resources
RESOURCE_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/resources/'

# bed dir
EXCL_BED = [MappabilityTrack(RESOURCE_PATH + 'hg38_centromere.bed.gz',      bed_buffer=10000),
            MappabilityTrack(RESOURCE_PATH + 'hg38_gap.bed.gz',             bed_buffer=1000),
            MappabilityTrack(RESOURCE_PATH + 'hg38_simpleRepeats.bed.gz',   bed_buffer=50),
            MappabilityTrack(RESOURCE_PATH + 'hg38_microsatellites.bed.gz', bed_buffer=10)]

HUMAN_CHR  = [str(n) for n in range(1,22+1)] + ['X', 'Y', 'M']
HUMAN_CHR += ['chr'+n for n in HUMAN_CHR]
HUMAN_CHR  = {n:True for n in HUMAN_CHR}

LEXICO_2_IND = {'chr1':1, 'chr2':2, 'chr3':3, 'chr10':10, 'chr11':11, 'chr12':12, 'chr19':19, 'chr20':20,
				'chr4':4, 'chr5':5, 'chr6':6, 'chr13':13, 'chr14':14, 'chr15':15, 'chr21':21, 'chr22':22,
				'chr7':7, 'chr8':8, 'chr9':9, 'chr16':16, 'chr17':17, 'chr18':18, 'chrX' :23, 'chrY' :24, 'chrM' :25}

def isValidCoord(mappability_track_list, myChr, myPos):
	if any([n.query(myChr,myPos) for n in mappability_track_list]):
		return False
	return True

RC_DICT = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def RC(s):
	return ''.join([RC_DICT[n] for n in s[::-1]])

# track for finding nearest-gene
TRANSCRIPT_TRACK = {}
f = gzip.open(RESOURCE_PATH + 'transcripts_hg38.bed.gz', 'r')
for line in f:
	splt = line.strip().split('\t')
	if splt[0] not in TRANSCRIPT_TRACK:
		TRANSCRIPT_TRACK[splt[0]] = []
	TRANSCRIPT_TRACK[splt[0]].append([int(splt[1]), int(splt[2]), splt[3], splt[4]])
f.close()

def get_nearest_transcript(myChr, pos, max_dist=20000):
	if myChr in TRANSCRIPT_TRACK:
		# lazy and slow, but it gets the job done!
		closest_dist = 99999999999
		closest_meta = ''
		for n in TRANSCRIPT_TRACK[myChr]:
			if pos >= n[0] and pos <= n[1]:
				closest_dist = 0
				closest_meta = [n[2], n[3]]
				break
			my_dist = min([abs(pos-n[0]), abs(pos-n[1])])
			if my_dist < closest_dist:
				closest_dist = my_dist
				closest_meta = [n[2], n[3]]
		if closest_dist <= max_dist:
			return (closest_dist, closest_meta)
	return None

# common SVs
COMMON_SVS = {}
f = gzip.open(RESOURCE_PATH + 'common_svs_hg38.bed.gz', 'r')
for line in f:
	splt = line.strip().split('\t')
	if splt[0] not in COMMON_SVS:
		COMMON_SVS[splt[0]] = []
	COMMON_SVS[splt[0]].append([int(splt[1]), int(splt[2]), splt[3]])
f.close()

SV_TYPE_MAP = {'deletion':                 ['DEL'],
               'herv deletion':            ['DEL'],
               'alu deletion':             ['DEL'],
               'line1 deletion':           ['DEL'],
               'sva deletion':             ['DEL'],
               'duplication':              ['DUP'],
               'copy number loss':         ['DEL'],
               'copy number gain':         ['DUP'],
               'copy number variation':    ['DUP','DEL'],
               'insertion':                ['INS'],
               'line1 insertion':          ['INS'],
               'alu insertion':            ['INS'],
               'sva insertion':            ['INS'],
               'mobile element insertion': ['INS']}

COLLAPSE_DUP_AND_INS = True

if COLLAPSE_DUP_AND_INS:
	for k in SV_TYPE_MAP.keys():
		if 'DUP' in SV_TYPE_MAP[k] and 'INS' not in SV_TYPE_MAP[k]:
			SV_TYPE_MAP[k].append('INS')
		if 'INS' in SV_TYPE_MAP[k] and 'DUP' not in SV_TYPE_MAP[k]:
			SV_TYPE_MAP[k].append('DUP')

def is_common_sv(myChr, pos_start, pos_end, my_type, max_endpoint_dist=50):
	if myChr in COMMON_SVS:
		c2 = sorted([pos_start, pos_end])
		for sv in COMMON_SVS[myChr]:
			if my_type in SV_TYPE_MAP[sv[2]]:	# same type
				c1 = sorted([sv[0], sv[1]])
				d1 = abs(c1[0]-c2[0])
				d2 = abs(c1[1]-c2[1])
				if d1 <= max_endpoint_dist and d2 <= max_endpoint_dist:	# beakpoints are same
					return True
	return False

#
UNMAPPED_ANCHOR = '*:0-0:FWD'
#
def get_sv_traceback(myAnchors, myType, a_i, rname_orig):
	out_sv = []
	for i in xrange(len(myAnchors)-1,-1,-1):
		current_anchor = myAnchors[i]
		current_type   = myType[i]

		if current_type == 'I':
			[a1,a2] = current_anchor.split('+')
			r1 = a1.split(':')[0:1] + [int(n) for n in a1.split(':')[1].split('-')] + a1.split(':')[2:3]
			r2 = a_i.split(':')[0:1] + [int(n) for n in a_i.split(':')[1].split('-')] + a_i.split(':')[2:3]
			r3 = a2.split(':')[0:1] + [int(n) for n in a2.split(':')[1].split('-')] + a2.split(':')[2:3]

			out_sv.append(['I', [tuple(r1), tuple(r2), tuple(r3)], rname_orig])

		elif current_type == 'Sl' or current_type == 'Sr':
			r1 = current_anchor.split(':')[0:1] + [int(n) for n in current_anchor.split(':')[1].split('-')] + current_anchor.split(':')[2:3]
			r2 = a_i.split(':')[0:1] + [int(n) for n in a_i.split(':')[1].split('-')] + a_i.split(':')[2:3]
			isSameOrientation = (r1[3] == r2[3])
			if not isSameOrientation:
				if current_type == 'Sl':
					r2 = [r2[0], r2[1], r2[2], r2[3]+'*']
				elif current_type == 'Sr':
					r2 = [r2[0], r2[2], r2[1], r2[3]+'*']

			if current_type == 'Sr':
				out_sv.append(['Sr', [tuple(r1), tuple(r2)], rname_orig])
			elif current_type == 'Sl':
				out_sv.append(['Sl', [tuple(r2), tuple(r1)], rname_orig])
	return out_sv

# [type, [regions], rname]
#            |
#            +--> [chr, p1, p2, orientation]
LARGE_NUMBER = (2**31) - 1
def event_dist(sv1, sv2):
	# event size
	if sv1[0][0] != sv2[0][0] or len(sv1[1]) != len(sv2[1]):
		return None
	# references
	for i in xrange(len(sv1[1])):
		r1 = getViralName(sv1[1][i][0])
		r2 = getViralName(sv2[1][i][0])
		if r1 != r2:
			return None

	if sv1[0] == 'I':
		d1 = abs(sv1[1][0][2] - sv2[1][0][2])
		d2 = abs(sv1[1][1][1] - sv2[1][1][1]) + abs(sv1[1][1][2] - sv2[1][1][2])
		d3 = abs(sv1[1][2][1] - sv2[1][2][1])
		return d1 + d2 + d3

	elif sv1[0][0] == 'S':

		isFwd = (sv1[1][0][3] == sv1[1][1][3] and sv2[1][0][3] == sv2[1][1][3])
		if True:
			d1 = abs(sv1[1][0][2] - sv2[1][0][2])
			d2 = abs(sv1[1][1][1] - sv2[1][1][1])
			return d1 + d2
		return None

		d1 = abs(sv1[1][0][1] - sv2[1][0][1])
		d2 = abs(sv1[1][1][2] - sv2[1][1][2])
		d1 = abs(sv1[1][sci[0]][sci[1]] - sv2[1][sci[4]][sci[5]])
		d2 = abs(sv1[1][sci[2]][sci[3]] - sv2[1][sci[6]][sci[7]])
		return d1 + d2
	#elif sv1[0] == 'Sr':
	#	d1 = abs(sv1[1][0][2] - sv2[1][0][2])
	#	d2 = abs(sv1[1][1][1] - sv2[1][1][1])
	#	return d1 + d2
	return None

# return all visited nodes given an adjaceny matrix and a starting index
def dfs_adj(A,start):
	stack   = [start]
	visited = set()
	while stack:
		vertex = stack.pop()
		visited.add(vertex)
		for i in xrange(len(A[vertex])):
			if i not in visited and A[vertex][i]:
				stack.append(i)
	return sorted(list(visited))

# return a list of all connected subgraphs
def get_connected_subgraphs(A):
	clusters_out = []
	not_visited = {i:True for i in xrange(len(A))}
	for k in not_visited.keys():
		if not_visited[k]:
			newCluster = dfs_adj(A,k)
			clusters_out.append([n for n in newCluster])
			for n in newCluster:
				not_visited[n] = False
	return clusters_out

# better, faster! (??)
def get_connected_subgraphs_rcm(A):
	A_sparse = csr_matrix(A)
	A_perm   = reverse_cuthill_mckee(A_sparse)
	out_clusters  = [[A_perm[0]]]
	current_start = 0
	for i in xrange(1,len(A_perm)):
		pass

CHR_WHITELIST  = [str(n) for n in xrange(1,22+1)] + ['X','Y','*']
CHR_WHITELIST += ['chr'+n for n in CHR_WHITELIST]
CHR_BLACKLIST  = ['M', 'chrM']
#BREAK_CHR      = ['MT','chrMT','M','chrM','2','chr2']
BREAK_CHR      = []

# adj clustering test
if __name__ == '__main__':

	test_A = [[0,1,0,0,0],
	          [1,0,0,0,0],
	          [0,0,0,1,1],
	          [0,0,1,0,1],
	          [0,0,0,1,0]]

	test_B = [[0,0,0,1,0,0,0],
	          [0,0,1,0,0,0,0],
	          [0,1,0,0,0,0,0],
	          [1,0,0,0,0,0,0],
	          [0,0,0,0,0,0,1],
	          [0,0,0,0,0,0,0],
	          [0,0,0,0,1,0,0]]

	print get_connected_subgraphs(test_A)
	print get_connected_subgraphs_rcm(np.array(test_A,dtype='uint8'))

	print get_connected_subgraphs(test_B)
	print get_connected_subgraphs_rcm(np.array(test_B,dtype='uint8'))

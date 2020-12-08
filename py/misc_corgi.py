import os
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

def isValidCoord(mappability_track_list, myChr, myPos):
	if any([n.query(myChr,myPos) for n in mappability_track_list]):
		return False
	return True

RC_DICT = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def RC(s):
	return ''.join([RC_DICT[n] for n in s[::-1]])

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

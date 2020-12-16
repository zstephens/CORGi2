import os
import sys
import json
import argparse

# absolute path to this script
SIM_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
sys.path.append(SIM_PATH+'/py/')

from misc_corgi import *

# parse input args
parser = argparse.ArgumentParser(description='process_svs.py')
parser.add_argument('-i',  type=str, required=True,  metavar='<str>', help="* input_vcf_dir/")
parser.add_argument('-o',  type=str, required=True,  metavar='<str>', help="* out_dir/")
parser.add_argument('-r',  type=str, required=True,  metavar='<str>', help="* ref.fa")
parser.add_argument('-t',  type=str, required=True,  metavar='<str>', help="* /path/to/truvari")
parser.add_argument('-b',  type=str, required=False, metavar='<str>', help="only_these_regions.bed", default=None)
args = parser.parse_args()

VCF_DIR = args.i
if VCF_DIR[-1] != '/':
	VCF_DIR += '/'
WORKING_DIR = args.o
if WORKING_DIR[-1] != '/':
	WORKING_DIR += '/'
if args.b != None:
	INCL_BED = ' --includebed ' + args.b
else:
	INCL_BED = ''
EXCL_BED = []
REF      = args.r
TRUVARI  = args.t

# feature vectors annotations:
BED_DIR  = '/Users/zach/Desktop/PacBioMetrics/annotation_bed_hg38/'


################################################################
######## raw SV caller output.
######## -- DIRECTLY FROM THE SOFTWARE, NO POST-PROCESSING YET
################################################################
####
####ALIGNER = 'ngmlr'
####
####RAW_DIR = 'ALL_DATA_COMPILATION/'
####
##### pbsv discover "${BAM}" "${PBSV_TEMP}"
##### pbsv call "${REF}" "${PBSV_TEMP}" "${VCF}" --ccs -t INS,DEL
######## using PBSV 0.12.0
####PBSV_VCF = RAW_DIR + 'hg002_'+ALIGNER+'_hg38_pbsv.vcf'
####
##### sniffles -s 3 --skip_parameter_estimation -m "${BAM_CALMD}" -v "${VCF}"
######## using the latest version of Sniffles as of 06/26/2019
####SNIFFLES_VCF = RAW_DIR + 'hg002_'+ALIGNER+'_hg38_sniffles.vcf'
####
##### svim alignment --cluster_max_distance 1.4 "${SVIM_TEMP_DIR}" "${BAM}"
######## SVIM 0.4.3
####SVIM_VCF = RAW_DIR + 'hg002_'+ALIGNER+'_hg38_svim.vcf'
####
##### NIST SVs VCF (crossmapped to hg38, sorted)
####RAW_GIAB = RAW_DIR + 'HG002_SVs_Tier1_v0.6_crossMap38_all.vcf'


############################################################
#### workflow helper functions
############################################################

# label collapsing one-liners, taken from: https://github.com/PacificBiosciences/sv-benchmark
def collapse_pbsv(fn_in, fn_out):
	cmd = 'bgzip -c ' + fn_in + ' > ' + fn_out
	os.system(cmd)
	cmd = 'tabix ' + fn_out
	os.system(cmd)

def collapse_sniffles(fn_in, fn_out):
	cmd  = 'cat ' + fn_in + ' | grep "^#" > ' + fn_out+'.temp_header'
	cmd += '; cat ' + fn_in + ' | grep -vE "^#" | grep \'DUP\|INS\|DEL\' | sed \'s/DUP/INS/g\' | sort -k1,1 -k2,2g > ' + fn_out+'.temp_data'
	os.system(cmd)
	cmd = 'cat ' + fn_out+'.temp_header' + ' ' + fn_out+'.temp_data' + ' | bgzip -c > ' + fn_out
	os.system(cmd)
	cmd = 'tabix ' + fn_out
	os.system(cmd)
	cmd = 'rm ' + fn_out+'.temp_header' + '; rm ' + fn_out+'.temp_data'
	os.system(cmd)

def collapse_svim(fn_in, fn_out, lenFilter=True):
	if lenFilter:
		awk_cond = 'if(($5=="<DEL>" || $5=="<INS>") && $6>40)'
	else:
		awk_cond = 'if($5=="<DEL>" || $5=="<INS>")'
	cmd = 'cat ' + fn_in + ' | sed \'s/INS:NOVEL/INS/g\' | sed \'s/DUP:INT/INS/g\' | sed \'s/DUP:TANDEM/INS/g\' | awk \'{ if($1 ~ /^#/) { print $0 } else { ' + awk_cond + ' { print $0 } } }\' > ' + fn_out+'.temp'
	os.system(cmd)
	cmd = 'bgzip -c ' + fn_out+'.temp' + ' > ' + fn_out
	os.system(cmd)
	cmd = 'tabix ' + fn_out
	os.system(cmd)
	cmd = 'rm ' + fn_out+'.temp'
	os.system(cmd)

def truvari(vcf_base, vcf_call, out_dir, out_log=None):
	if out_dir[-1] == '/':
		out_dir = out_dir[:-1]
	cmd = TRUVARI + ' -f ' + REF + ' -b ' + vcf_base + ' -c ' + vcf_call + ' -o ' + out_dir + ' -r 1000 -p 0.00 --passonly' + INCL_BED
	if out_log != None:
		cmd += ' > ' + out_log + ' 2>&1'
	os.system(cmd)

def bgzip_tabix(fn_in, fn_out):
	cmd = 'bgzip -c ' + fn_in + ' > ' + fn_out
	os.system(cmd)
	cmd = 'tabix ' + fn_out
	os.system(cmd)

def exists_and_is_nonZero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

def sort_vcf(fn_in, fn_out):
	header = ''
	dat_sort = []
	f = open(fn_in, 'r')
	for line in f:
		if line[0] == '#':
			header += line
		else:
			splt = line.split('\t')
			if splt[0] in LEXICO_2_IND:
				dat_sort.append([LEXICO_2_IND[splt[0]], int(splt[1]), line])
	dat_sort = sorted(dat_sort)
	f.close()
	f = open(fn_out, 'w')
	f.write(header)
	for n in dat_sort:
		f.write(n[2])
	f.close()

############################################################
#### the actual workflow itself
############################################################

INPUT_VCF   = sorted([VCF_DIR + n for n in os.listdir(VCF_DIR) if n[-4:] == '.vcf'])
INPUT_NAMES = [n.split('/')[-1].split('_')[0] for n in INPUT_VCF]

WDIR_0 = WORKING_DIR + '00_input_vcf/'
WDIR_1 = WORKING_DIR + '01_truvari_pairwise/'
WDIR_2 = WORKING_DIR + '02_unique_svs/'

TRUVARI_LOG = WORKING_DIR + 'truvari.log'

cmd = 'mkdir -p ' + '; mkdir -p '.join([WORKING_DIR, WDIR_0, WDIR_1, WDIR_2])
os.system(cmd)

#CLASS_LABELS  = WORKING_DIR + 'sv_labels.tsv'
#ANNOTATED_SVS = WORKING_DIR + 'sv_annotated.tsv'
#CLASS_COUNTS  = WORKING_DIR + 'sv_counts.json'

# prepare for truvari
VCF_FN_00 = [WDIR_0 + n + '.vcf.gz' for n in INPUT_NAMES]
print 'preparing vcf for truvari...'
for i in xrange(len(INPUT_VCF)):
	if exists_and_is_nonZero(VCF_FN_00[i]) == False:
		collapse_pbsv(INPUT_VCF[i], VCF_FN_00[i])

# truvari
print 'running truvari...'
VCF_FN_02 = [WDIR_2 + n + '_unique.vcf' for n in INPUT_NAMES]
for i in xrange(len(INPUT_VCF)):
	if exists_and_is_nonZero(VCF_FN_02[i]) == False:
		current_vcf = VCF_FN_00[i]
		for j in xrange(len(INPUT_VCF)):
			if i != j:
				print (i,j)
				tv_dir = WDIR_1 + INPUT_NAMES[i] + '_' + INPUT_NAMES[j]
				vcf_2  = VCF_FN_00[j]
				truvari(current_vcf, vcf_2, tv_dir, TRUVARI_LOG)
				current_vcf      = tv_dir + '/fn.vcf'
				current_vcf_sort = tv_dir + '/fn_sort.vcf'
				current_vcf_gz   = tv_dir + '/fn_sort.vcf.gz'
				if exists_and_is_nonZero(current_vcf):
					sort_vcf(current_vcf, current_vcf_sort)
					bgzip_tabix(current_vcf_sort, current_vcf_gz)
					current_vcf = current_vcf_gz
				else:
					print 'where did current_vcf go?', current_vcf
					exit(1)
		os.system('mv ' + current_vcf[:-3] + ' ' + VCF_FN_02[i])

# annotation: gene hits?
for i in xrange(len(INPUT_VCF)):
	f = open(VCF_FN_02[i], 'r')
	col_info = None
	col_form = None
	col_samp = None
	for line in f:
		if line[0] == '#':
			if line[1] != '#':
				splt = line[1:].strip().split('\t')
				col_info = splt.index('INFO')
				col_form = splt.index('FORMAT')
				col_samp = len(splt) - 1
		else:
			splt = line.strip().split('\t')
			#
			# are we supported by more than 1 read?
			#
			splt2  = splt[col_form].split(':')
			splt3  = splt[col_samp].split(':')
			col_ad = splt2.index('AD')
			readcounts = [int(n) for n in splt3[col_ad].split(',')]
			if len(readcounts) > 2:
				print 'multiple alleles, skipping this for now...'
				continue
			elif len(readcounts) == 2:
				if readcounts[1] <= 1:
					continue
			#
			# are we a common SV?
			#
			splt4 = splt[col_info].split(';')
			my_type = None
			my_len  = None
			my_end  = None
			for n in splt4:
				if n[:7] == 'SVTYPE=':
					my_type = n[7:]
				if n[:6] == 'SVLEN=':
					my_len = int(n[6:])
				if n[:4] == 'END=':
					my_end = int(n[4:])
			if my_type == None:
				print 'SVTYPE not found in info? skipping...'
				continue
			if my_type == 'INS':
				my_end = int(splt[1])
			elif my_type == 'DEL' and my_len < 0:
				my_end = int(splt[1]) - my_len
			elif my_type == 'DUP' and my_len > 0:
				my_end = int(splt[1]) + my_len
			elif my_end == None:
				print 'what kind of weird SV are you?'
				print splt
				exit(1)
			if is_common_sv(splt[0], int(splt[1]), my_end, my_type):
				continue
			#
			# are we near a gene?
			#
			nearest_hit = get_nearest_transcript(splt[0], int(splt[1]))
			if nearest_hit != None:
				print splt
				print nearest_hit[0], nearest_hit[1][1]
	f.close()
	break

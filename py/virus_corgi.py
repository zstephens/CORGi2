import os
import json

# absolute path to resources
RESOURCE_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/resources/'

# get viral accession ids and whatnot
f = open(RESOURCE_PATH + 'taxid10239.nbr', 'r')
ID_TO_ACCESSION       = {}
ACCESSION_TO_TAXONOMY = {}
for line in f:
	if line[0] != '#':
		splt = line.strip().split('\t')
		ID_TO_ACCESSION[splt[1]] = splt[0]
		ACCESSION_TO_TAXONOMY[splt[0]] = splt[4]
f.close()

# read in viral short-hand that I used for the long read workflow
f = open(RESOURCE_PATH + 'HumanViral_Reference_12-12-2018_simpleNames.json', 'r')
VIRAL_NAME_DICT = json.load(f)
f.close()

def isVirus(refname):
	if 'virus' in refname:
		return True
	else:
		return False

def getViralName(refname):
	if refname in VIRAL_NAME_DICT:
		refname = VIRAL_NAME_DICT[refname].split(' ')[0]
	if refname in ID_TO_ACCESSION:
		refname = ID_TO_ACCESSION[refname]
	if refname in ACCESSION_TO_TAXONOMY:
		refname = ACCESSION_TO_TAXONOMY[refname]
	return refname

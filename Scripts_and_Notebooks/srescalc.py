#!/usr/bin/env python

"""
This code is part of the general resources library in the Yu Lab and was originally
written by Michael Meyer.
"""

""" MJM March 2015. Revamped June 2015. Calculate for single structure, only chains specified. If only two chains present in structure, no need to specify.
    srescalc.py: Calculates interface residues between specified chains in PDB files."""

import os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from tempfile import mkdtemp
from shutil import rmtree, copyfile
from sys import stdin
from re import split as regexsplit
from collections import defaultdict
from helper_functions import naccess, open_pdb, natural_keys
import tarfile

import atexit

# Added to make sure scratch directory is deleted even if
# the program ends prematurely (when there are errors in the
# naccess script, it exits automatically)
def exit_handler(scratchDir):
   try:
      rmtree(scratchDir)
   except:
      pass


#---------------------------------- PARSE COMMAND-LINE ARGUMENTS ----------------------------------
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, description='Calculates surface residues in all chains in isolation in a PDB-formatted file. OR specify which chain to calculate.')
parser.add_argument('structure', nargs='?', default=stdin, help='File name of input PDB file (.pdb or .ent), or gzipped pdb file (.gz), or 4-letter pdb ID. If no argument given, reads structure on STDIN.')
parser.add_argument('-c', help='Calculate surfaces residues only for specified chain in isolation.', required=False)
parser.add_argument('-complex', help='Calculate surface residues on intact complex (not with each chain in isolation). Flag with paramter ALL will calculate sres in all chains, or with argument for specific chains separated by forward slashes (i.e. A/B/C/D).', default='Isolation', required=False)
parser.add_argument('-uSASA', help='Minimum percent solvent accessibility of unbound residues to be considered a surface residue. Integer [0-100].', type=int, default=15)

parser.add_argument('-o', '--output', choices=['list', 'residue_stats'], help='Format output options.', default='list')


args = parser.parse_args()

#---------------------------------- OPEN STRUCTURE ----------------------------------

# EDIT 2021_03_31 by Shobhita Gupta
# Modified the pdbinfile to be a list of pdb input files to
# account for cryoEM structures.
# Parsed the chainmapping file and stored it in a dictionary
# chainmapping.
# Moved up the creation of the scratchDir so that we can
# extract the tarball there for cryoEM structures.

# ORIGINAL
# pdbinfile = open_pdb(args.structure)

# NEW
pdbinfiles = [open_pdb(args.structure, verbose=False)]
chainmapping = None

# CREATE SCRATCH SPACE FOR WRITING INTERMEDIATE FILES
scratchDir = mkdtemp()
atexit.register(exit_handler, scratchDir)
tempfilename = os.path.join(scratchDir, 'temp.txt')

# Check for large cryoEM structures
if pdbinfiles[0] == None:
	sys.exit("Error: PDB Structure not found.")
		
#---------------------------------- READ STRUCTURE ----------------------------------

# READ STRUCTURE FILE AND SPLIT INTO CHAINS

# EDIT 2020_02_19 by Shayne Wierbowski
# Modified the structure parsing to only use the last model
# as is done in irescalc.py. NACCESS can fail when there are
# many models because all of the atoms get stacked together
# in one model.

# ORIGINAL
#pdbatomdict = defaultdict(list)
#
#for line in pdbinfile:
#	if line is not None and line[:4] == "ATOM":
#		chain = line[21]
#		pdbatomdict[chain].append(line)
#pdbinfile.close()

pdbatomdict = {}
curchain = ''
mapped_chain = ''

# EDIT 2021_03_31 by Shobhita Gupta
# Modified to iterate over all pdb input files
# Subsequently, modified to store ATOM lines using the 
# mapped_chain as the key instead of the chain from the PDB file
# NEW
for pdbinfile in pdbinfiles:
	curfile = os.path.basename(pdbinfile.name)
	for line in pdbinfile:
		# Handles edge case for PDBs with multiple models (but only one chain)
		# Will default to use atoms from the last model (all earlier models get
		# their atomdict replaced when subsequent models come around)
		if line[:6] == "ENDMDL":
			curchain = ''
		if line is not None and line[:4] == "ATOM":
			if curchain != line[21]:
				curchain = line[21]
				if chainmapping != None:
					mapped_chain = chainmapping[curfile][curchain]
				else:
					mapped_chain = curchain
				pdbatomdict[mapped_chain] = [line]
			else:
				pdbatomdict[mapped_chain].append(line)
	pdbinfile.close()


# CHECK THAT SELECTED CHAINS EXIST IN GIVEN PDB FILE
user_chains = set()
if args.complex != 'Isolation':
	if args.complex == 'ALL':
		user_chains = set(pdbatomdict.keys())
	else:
		user_chains = set(args.complex.split('/'))
elif args.c != None:
	user_chains = set([args.c])

if len(user_chains) > 0 and user_chains.intersection(set(pdbatomdict.keys())) != user_chains:
	sys.exit('Error: No chains in file match filter. Remember, PDB chains are case-sensitive.\nAvailable Chains: %s\nSelected Chains: %s' %(', '.join(pdbatomdict.keys()), ', '.join(list(user_chains))))


#---------------------------------- PREPARE INDIVIDUAL CHAIN FILES FOR NACCESS ----------------------------------

pdbchainlist = sorted(pdbatomdict.keys())
temp_pdb_files = []

if args.complex == 'Isolation':
	for h in pdbchainlist:	
		if args.c != None and h != args.c:
			continue
		
		if h==' ':
			temp_file_name = '_.pdb'
		else:
			temp_file_name = h+'.pdb'
		
		temp_pdb_files.append(os.path.join(scratchDir,temp_file_name))
		tempoutfile = open(temp_pdb_files[-1], 'w')
		
		for line in pdbatomdict[h]:
			tempoutfile.write(line)
		tempoutfile.close()
	
else:
	if args.complex == 'ALL':
		complex_chains = pdbchainlist
	else:
		complex_chains = args.complex.split('/')
	
	temp_pdb_files.append(os.path.join(scratchDir,'complex.pdb'))
	tempoutfile = open(temp_pdb_files[-1], 'w')
	
	for h in complex_chains:
		for line in pdbatomdict[h]:
			tempoutfile.write(line)
	
	tempoutfile.close()


#---------------------------------- RUN NACCESS ----------------------------------

# EDIT 2021_03_31 by Shobhita Gupta
# Split naccess output by filename so we can recover the chains later.

# ORIGINAL
#naccess_output = []
#for pdb_file in temp_pdb_files:
#	naccess_output += naccess(pdb_file)

# NEW
naccess_output = {}
for pdb_file in temp_pdb_files:
	naccess_output[pdb_file] = naccess(pdb_file)

#---------------------------------- PARSE NACCESS OUTPUT ----------------------------------

# asa = accessible surface area
asadict = defaultdict(list)

# EDIT 2021_03_31 by Shobhita Gupta
# Iterate through all chain files and split out the chain id
# instead of relying on naccess for large cryoEM structures.
for filename in naccess_output:
	for line in naccess_output[filename]:
		if line[:3] != 'RES':
			continue
		
		aa = line[3:7].strip()
		chain = line[7:9].strip()
		residue_index = line[9:13].strip()
		relative_perc_accessible = float(line[22:28])
		
		if relative_perc_accessible < args.uSASA:    #percent solvent accessible >= 15% ???
			continue
		
		# NEW: if this is a large cryoEM structure, then we get the chain name from the file.
		if chainmapping != None:
			chain = filename.split(".")[0].split("/")[-1]
		
		asadict[chain].append((residue_index, aa, relative_perc_accessible))

#---------------------------------- FORMAT AND PRINT OUTPUT ----------------------------------

#Format and print output
for k in sorted(asadict.keys(), key=lambda x: (len(x), x.isdigit(), x)):
	
	indices = sorted([res[0] for res in asadict[k]], key=natural_keys)
	
	if args.c != None:
		if args.output == 'list':
			print ','.join(indices)
		elif args.output =='residue_stats':
			for res in asadict[k]:
				print '%s\t%s\t%s\t%s' %(k, res[0], res[1], res[2])

	else:
		if args.output == 'list':
			print k + '\t'+ ','.join(indices)
		elif args.output =='residue_stats':
			for res in asadict[k]:
				print '%s\t%s\t%s\t%s' %(k, res[0], res[1], res[2])

#cleanup and goodbye
rmtree(scratchDir)
exit()

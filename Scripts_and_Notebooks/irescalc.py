#!/usr/bin/env python

"""
This code is part of the general resources library in the Yu Lab and was originally
written by Michael Meyer.
"""


""" MJM edited version. Calculate for single structure, only chains specified. If only two chains present in structure, not need to specify.
	irescalc.py: Calculates interface residues between specified chains in PDB files.
	
	Update 8/28/2014: Input of structure can be file path, PDB ID (reads from resources), or a full structure on standard input"""

import argparse, os, sys
from tempfile import mkdtemp
from shutil import rmtree
from helper_functions import open_pdb, naccess

## I/O arguments at command line
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Calculates interface residues between two chains in a PDB-formatted file. If more than two chains are available in file, specify which chains to use.')
parser.add_argument('structure', nargs='?', default=sys.stdin, help='File name of input PDB file (.pdb or .ent), or gzipped pdb file (.gz), or 4-letter pdb ID. If no argument given, reads structure on STDIN.')
parser.add_argument('-c1', help='First of chain pair for which to calculate interface residues', required=False)
parser.add_argument('-c2', help='Second of chain pair for which to calculate interface residues', required=False)
parser.add_argument('-uSASA', help='Minimum percent solvent accessibility of an unbound residue to be considered a surface residue. Integer [0-100]. Eliminate criteria with value of 0.', type=int, default=15)
parser.add_argument('-dSASA', help='Minimum change in solvent accessible surface area (in Angstroms squared) for a surface residue to be considered an interface residue. Float [0-inf]', type=float, default=1.0)


parser.add_argument('-o', '--output', choices=['list', 'residue_stats'], help='Format output options.', default='list')

args = parser.parse_args()

chain1 = args.c1
chain2 = args.c2

#---------------------------------- OPEN STRUCTURE ----------------------------------

pdbinfile = open_pdb(args.structure)

#------------------------------ PARSE PDB FILE ---------------------------------

pdbatomdict = {}
curchain = ''

for line in pdbinfile:
	if line is not None and line[:4] == "ATOM":
		if curchain != line[21]:
			curchain = line[21]
			pdbatomdict[curchain] = [line]
		else:
			pdbatomdict[curchain].append(line)
pdbinfile.close()


# CREATE SCRATCH SPACE FOR WRITING INTERMEDIATE FILES
scratchDir = mkdtemp()
tempfilename = os.path.join(scratchDir, 'temp')

## writing PDB files, and then generating RSA, ASA and LOG files by using Naccess to calculate solvent accessibilities, and finally using RSA files to calculate interface residues
asadict = {}
pdbchainlist = sorted(pdbatomdict.keys())


# EDIT: Made By Shayne (July 12, 2018)
# I have moved this block of code to identify the chains being tested prior to calculating
# the unbound ASA for each chain. I have also modified the calculation of unbound ASA
# to only calculate for the desired two chains rather than for every chain in the structure.
# This should provide a performance boost particularly for excessively large complexes.
#
# NOTE: In the future I would like to include unbound ASA as a database that can be pre-computed
# to avoid needed to calculate it for individual chains over and over again. 

#Pick pair of chains to calc ires based on user input or default to only two chains in structure.
if chain1 == None or chain2 == None:  #chains not specified
	if len(pdbchainlist) == 2:
		i = 0; j = 1
	else:
		sys.exit('No chains specified for structure with > 2 available chains.')
		
else: #specific chains specified
	try:
		i = pdbchainlist.index(chain1)
		j = pdbchainlist.index(chain2)
	except:
		sys.exit('One or both chains %s and %s not found in file' %(chain1, chain2))
		

for chain in pdbchainlist:
	# EDIT: Made By Shayne (July 12, 2018)
	# I have modified this loop to only calculate ASA for the two input chains
	# See above.
	if(not pdbchainlist.index(chain) in [i, j]):
		continue
	
	asadict[chain] = {}
	
	temp_pdb_file = tempfilename+'.pdb'
	
	#write unbound form of chain for naccess calculations
	tempoutfile = open(temp_pdb_file,'w')
	for line in pdbatomdict[chain]:
		tempoutfile.write(line)
	tempoutfile.close()
	
	#---------------------------------- RUN NACCESS ----------------------------------
	
	naccess_output = naccess(temp_pdb_file, parsed=True)
	
	for res in naccess_output:
		if res['all_atoms_rel'] >= args.uSASA:      # relative percent SASA >= 15% indicates surface residue
			residue_key = (res['res'], res['chain'], res['res_num'])
			asadict[chain][residue_key] = res['all_atoms_abs']   #save absolute SASA to compare with bound form



#-----------------------------------------------
	
#Calculalate interface residues for the given chain pair	
asadict[pdbchainlist[i]+pdbchainlist[j]] = {}

#write temporary file with only those two chains for processing with naccess
temp_pdb_file = tempfilename+'.pdb'
	
tempoutfile = open(temp_pdb_file,'w')
for line in pdbatomdict[pdbchainlist[i]]:
	tempoutfile.write(line)
for line in pdbatomdict[pdbchainlist[j]]:
	tempoutfile.write(line)
tempoutfile.close()

#---------------------------------- RUN NACCESS ----------------------------------
naccess_output = naccess(temp_pdb_file, parsed=True)

intreslist = {}
intreslist[pdbchainlist[i]] = []
intreslist[pdbchainlist[j]] = []

for res in naccess_output:	
	residue_key = (res['res'], res['chain'], res['res_num'])
	chain = res['chain']
	
	if residue_key not in asadict[chain]:
		continue
	
	res_dSASA = abs(asadict[chain][residue_key] - res['all_atoms_abs'])
	
	if res_dSASA >= args.dSASA:    #change in solvent accessible surface area >= 1 A squared
		intreslist[chain].append((chain, res['res_num'], res['res'], res_dSASA))

#---------------------------------- PRINT OUTPUT ----------------------------------


if args.output == 'list':
	print ','.join([q[1] for q in intreslist[pdbchainlist[i]]])
	print ','.join([q[1] for q in intreslist[pdbchainlist[j]]])
elif args.output == 'residue_stats':
	for res in intreslist[pdbchainlist[i]]:
		print '%s\t%s\t%s\t%s' %tuple(res)
	for res in intreslist[pdbchainlist[j]]:
		print '%s\t%s\t%s\t%s' %tuple(res)

#cleanup
rmtree(scratchDir)

# Sub Helper Libraries
import colors

from Bio import Entrez
import re
import subprocess as sp
import numpy as np
import pandas as pd
import math
import os
import sys
import glob
import requests
import time
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import random
import string
import inspect
import threading
import copy
from multiprocessing import Process, Pool, Queue, Value, Lock
from collections import defaultdict
import MySQLdb
import html2text
import mjm_tools
import traceback
import argparse
import getpass
import datetime
from difflib import SequenceMatcher
import hashlib
import primer3
import itertools

import scipy.stats as st


import jfb_tools

from bitarray import bitarray

MAXTHREADS = 20

user = getpass.getuser()
try:
	import authentication
	gris = open("[REDACTED_PATH]/bin/authentication.txt", "r")
	foo = gris.readline().strip()
	bar = eval(gris.readline().strip())
	gris.close()
except:
	foo = "nice"
	bar = "try"

startingDir = "/home/{0}/.tmp/".format(user)
def setStartingDir(s):
	startingDir = s
	if(not os.path.exists(startingDir)):
		print "Automatically creating tmp folder at {0} to store intermediate files.".format(startingDir)
		print "It is strongly recommended that you set up a recurring script in your crontab to clear out this directory periodically."
		print "You can edit the jobs that should be executed as part of you crontab using \"crontab -e\""
		print "For instance, my crontab contains the following two lines to handle this case..."
		print
		print "# Daily Cleanup of TMP Files Older than 7 Days"
		print "0 0 * * * find ~/.tmp/ -type f -mtime +7 -name '*' -print0 | xargs -r0 rm --"
		os.system("mkdir -p {0}".format(startingDir))
# FUNCTION END
setStartingDir(startingDir)

codon_mapping = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CAA': 'Q', 'AGT': 'S', 'CAG': 'Q', 'CCG': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TAG': '*', 'GGA': 'G', 'TAA': '*', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TGA': '*', 'TTG': 'L', 'CGT': 'R', 'TGG': 'W', 'CGC': 'R'}
codon_mapping_U = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CAA': 'Q', 'AGT': 'S', 'CAG': 'Q', 'CCG': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TAG': '*', 'GGA': 'G', 'TAA': '*', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TGA': 'U', 'TTG': 'L', 'CGT': 'R', 'TGG': 'W', 'CGC': 'R'}
# - Convert DNA Seq to AA Seq
# - Accepts as parameters...
#   @seq - Seqnence to convert
#   @codong_mapping - Dictionary for mapping
#   @include_stop_codon - Boolean indicating whether stop codon should be reported
#   @return_first - Boolean indicating whether to return first hit or longest hit
#   @return_full - Whether or not the full sequence (not necessarilly starting from ATG) should be returned
#   @include_U - Whether or not intermediate stop_codons should be converted to U (selenocysteine)
# - Returns the longest continuous transcript in any frame of the given sequence
# - Written for hw2 of BTRY 4381 03/01/17 by Shayne Wierbowski
def DNA2AA(seq, codon_mapping=codon_mapping, include_stop_codon=False, return_first=False, return_start=False, return_full=False, include_U=False):
	seq = seq.upper()
	start_codon = [x for x in codon_mapping if codon_mapping[x] == "M"][0]
	if(return_full):
		start_pos = [0]
	elif(return_first):
		start_pos = [seq.find(start_codon)]
		if(start_pos[0] == -1):
			return ""
	else:
		start_pos = [m.start() for m in re.finditer(start_codon, seq)]
	
	
	if(include_stop_codon):
		if(include_U):
			transcripts = ["".join([codon_mapping_U[seq[i:i+3]] for i in range(start, len(seq)-2, 3)]) for start in start_pos]
			def do(x):
				if(x[-1] == "U"):
					return x[:-1] + "*"
				return x
			# FUNCTION END
			transcripts = [do(x) for x in transcripts]
		else:
			transcripts = ["".join([codon_mapping[seq[i:i+3]] for i in range(start, len(seq)-2, 3)]).split("*")[0] + "*" for start in start_pos]
	else:
		if(include_U):
			transcripts = ["".join([codon_mapping_U[seq[i:i+3]] for i in range(start, len(seq)-2, 3)]) for start in start_pos]
			def do(x):
				if(x[-1] == "U"):
					return x[:-1]
				elif(x[-1] == "*"):
					return x[:-1]
				return x
			# FUNCTION END
			transcripts = [do(x) for x in transcripts]
		else:
			transcripts = ["".join([codon_mapping[seq[i:i+3]] for i in range(start, len(seq)-2, 3)]).split("*")[0] for start in start_pos]
	
	
	if(len(transcripts) == 0):
		return ""
	
	#if(return_first):
	#	if(len(transcripts) != 0):
	#		return transcripts[0]
	#	else:
	#		return ""
	#else:
	#	if(return_start):
	#		return max(zip(transcripts, start_pos), key=lambda x: len(x[0]))
	#	return(max(transcripts, key=len))
	if(return_start):
		return max(zip(transcripts, start_pos), key=lambda x: len(x[0]))
	else:
		return(max(transcripts, key=len))
# FUNCTION END

# - Calculates the degree of a node in a PPI network
# - Accepts as parameters...
#   @network - A PPI network defined as a list of tuples where
#              each tuple contains the two interacting proteins
#   @node - The node for which the degree should be calculated
# - Returns the degree of the requested node in the given network
# - Written for hw3 of BTRY 4381 03/15/17 by Shayne Wierbowski
def degree(network, node):
	return len(set([i for edge in network for i in edge if node in edge])) - (not has_homodimer(network, node))
# FUNCTION END

# - Helper function for degree, checks whether a node exists as a homodimer
# - Accepts as parameters...
#   @network - A PPI network defined as a list of tuples where
#              each tuple contains the two interacting proteins
#   @node - The node in question
# - Returns a boolean indicating whether or not the requested node is
#   present in the given network as a homodimer
# - Written for hw3 of BTRY 4381 03/15/17 by Shayne Wierbowski
def has_homodimer(network, node):
	return len([edge for edge in network if len(set(edge).union([node])) == 1]) != 0
# FUNCTION END

# - Provides basic stats for a PPI network
# - Accepts as parameters...
#   @network - A PPI network defined as a list of tuples where
#              each tuple contains the two interacting proteins
# - Returns in order, # nodes, # edges, # homodimers, average degree
# - Written for hw3 of BTRY 4381 03/15/17 by Shayne Wierbowski
def net_stats(network):
	nodes = set([i for edge in network for i in edge])
	edges = set(["\t".join(sorted([str(edge[0]), str(edge[1])])) for edge in network])
	homodimers = set([edge for edge in edges if len(set(edge.split("\t"))) == 1])
	degrees = [degree(network, i) for i in nodes]
	return len(nodes), len(edges), len(homodimers), np.mean(degrees)
# FUNCTION END

# - Gets the orthologs for a specified protein
# - Accepts as parameters...
#   @orthologs - A dictionary representation of a set of orthologs
#                where each key is a protein in one species and the
#                values are a list of orthologous proteins in another
#   @node - The query protein
# - Returns the list of orthologs for the given protein (the empty list
#   if the protein is not present in orthologs)
# - Written for hw3 of BTRY 4381 03/15/17 by Shayne Wierbowski
def getOrth(orthologs, node):
	try:
		return orthologs[node]
	except KeyError:
		return []
# FUNCTION END

# List of file for cleanup
myCleanupRemove = []

# Semaphore to protect myCleanupRemove
myCleanupSema = threading.BoundedSemaphore()

# - Adds the specified file to the cleanup list with specified priority
# - Accepts as parameters...
#   @baseName - The base name for the temp file (default to "tmp")
#   @priority - The cleanup priority for the given file (where higher number
#               priority files are retained through less thorough cleanups)
# - Written for general utility purposes 03/15/17 by Shayne Wierbowski
def myAddRemove(file, priority=0):
	global myCleanupRemove
	
	myCleanupSema.acquire()
	myCleanupRemove.append((os.path.abspath(file), priority))
	myCleanupSema.release()
# FUNCTION END

# - Removes all files in myCleanupRemove that are less than or equal to
#   the given priority
# - Accepts as parameters...
#   @priority - The cleanup priority (all files at or below this priority)
#               will be deleted
#   @all - A boolean indicating whether or not all files should be deleted
# - Written for general utility purposes 03/15/17 by Shayne Wierbowski
def myCleanup(priority=0, all=False):
	if(all):
		priority = sys.maxint
	
	global myCleanupRemove
	
	myCleanupSema.acquire()
	cleanup = [i[0] for i in myCleanupRemove if i[1] <= priority]
	for f in cleanup:
		myCleanupFile(f)
	myCleanupSema.release()
# FUNCTION END

# - Removes the specified file
# - Accepts as parameters...
#   @file - The file to delete
# - Written for general utility purposes 03/15/17 by Shayne Wierbowski
def myCleanupFile(file):
	if(os.path.exists(file)):
		#os.system("rm {0}".format(file))
		call("rm {0}".format(file), makeTmp=False)
	if(file in myCleanupRemove):
		myCleanupRemove.remove(file)
# FUNCTION END

# Sub Name for all Temp Files produced by a specific python instance
tmpSubName = ""

# Semaphore to protect temp reservation
tmpSema = threading.BoundedSemaphore()

# - Reserves the first available temporary file name, tmp_X
# - Accepts as parameters...
#   @startingDir - Directory at which to reserve temp file (default to cwd)
#   @baseName - The base name for the temp file (default to "tmp")
#   @addRemove - A boolean indicating whether or not the newly created
#                temp file should be added to the list for automatic
#                removal on cleanup
#   @priority - The cleanup priority for the newly created temp file (where
#               higher number priority files are retained through less thorough
#               cleanups)
#   @N - Number of random characters to use in extra security measure
# - Returns the name of the reserved temporary
# - Written for general utility purposes 03/15/17 by Shayne Wierbowski
def reserveTemp(startingDir=startingDir, baseName="tmp", addRemove=True, priority=0, N=6):
	global tmpSubName
	global tmpSema
	
	if(startingDir != "[REDACTED_PATH]/.tmp/"):
		if(startingDir[-1] != "/"):
			startingDir += "/"
	
	''' DITCHED IN FAVOR OF JUST ADDING A RANDOM NUMBER TO FILENAMES
	# Wait for "Mutext" to unlock
	while(os.path.exists("{0}_{1}_inUse1".format(startingDir, basename))):
		while(os.path.exists("{0}_{1}_inUse2".format(startingDir, basename))):
			time.sleep(0.1)
		call("touch {0}_{1}_inUse2".format(startingDir, basename), makeTmp=False)
		time.sleep(0.1)
	call("touch {0}_{1}_inUse1".format(startingDir, basename), makeTmp=False)
	'''
	
	# Reserve a Temp Sub Name for the currently running python shell
	while(tmpSubName == ""):
		tmpSema.acquire()
		if(tmpSubName != ""):
			tmpSema.release()
			break
		i = 1
		while(len(ls("{0}tmp_{1}_*".format(startingDir, i), fromcwd=False)) > 0):
			# Cleanup Temp Sub Names that are at least a week old
			#if(os.path.exists("{0}_{1}_0".format(startingDir, i)) and time.time() - os.stat("{0}_{1}_0".format(startingDir, i)).st_mtime > 60*60*24*7):
			#	call("rm {0}_{1}_0".format(startingDir, i), makeTmp=False)
			try:
				if(i < max([int(float(x[x.rfind("tmp_"):].split("_")[1])) for x in ls("{0}tmp_*".format(startingDir), fromcwd=False)])):
					i = max([int(float(x[x.rfind("tmp_"):].split("_")[1])) for x in ls("{0}tmp_*".format(startingDir), fromcwd=False)]) + 1
				else:
					i += 1
			except:
				i += 1
		ranStr = ''.join([random.choice(string.ascii_uppercase + string.digits) for _ in range(N)])
		call("touch {0}_{1}_0".format(startingDir, i), makeTmp=False)
		tmpSubName = "{0}_{1}".format(str(i), ranStr)
		call("chmod 777 {0}_{1}_0".format(startingDir, i), makeTmp=False)
		tmpSema.release()
	
	tmpSema.acquire()
	i = 1
	while(len(ls("{0}{1}_{2}_{3}_*".format(startingDir, baseName, tmpSubName, i), fromcwd=False)) > 0):
		#while(os.path.exists("{0}{1}_{2}_{3}".format(startingDir, baseName, tmpSubName, i))):
		try:
			if(i < max([int(float(x[x.rfind(baseName):].split("_")[3])) for x in ls("{0}{1}_{2}_*".format(startingDir, baseName, tmpSubName), fromcwd=False)])):
				i = max([int(float(x[x.rfind(baseName):].split("_")[3])) for x in ls("{0}{1}_{2}_*".format(startingDir, baseName, tmpSubName), fromcwd=False)]) + 1
			else:
				i += 1
		except:
			i += 1
	ranStr = ''.join([random.choice(string.ascii_uppercase + string.digits) for _ in range(N)])
	call("touch {0}{1}_{2}_{3}_{4}".format(startingDir, baseName, tmpSubName, i, ranStr), makeTmp=False)
	tempName = "{0}{1}_{2}_{3}_{4}".format(startingDir, baseName, tmpSubName, i, ranStr)
	call("chmod 777 {0}".format(tempName), makeTmp=False)
	tmpSema.release()
	#os.system("touch {0}".format(tempName))
	
	#myCleanupFile("{0}_{1}_0".format(startingDir, tmpSubName))
	
	if(addRemove):
		myAddRemove(tempName, priority)
	return tempName
# FUNCTION END

# - Times a function or evaluable string
# - Accepts as parameters...
#   @fn - The function or expression to execute
#   @times - The number of times to repeate execution
#   @HR - Flag indicating whether to return HR time or raw time
# - Returns the total elapsed time to execute
# - Written for general utility purposes 07/21/17 by Shayne Wierbowski
def mytime(fn, times=1, HR=True):
	if(type(fn) == str):
		start = time.time()
		for i in range(times):
			eval(fn)
		end = time.time()
	else:
		start = time.time()
		for i in range(times):
			fn()
		end = time.time()
	if(HR):
		return HRTime(end - start)
	else:
		return end - start
# FUNCTION EDN

# - Reads the lines for s specified file
# - Accepts as parameters...
#   @fileName - The base name of the file
#   @startingDir - Directory at which the file is located (default to cwd)
#   @skip - Number of lines to skip
#   @comment - Delimeter for a comment line to be ignored
#   @split - Boolean flag for whether or not to plit output by lines
#   @sep - Delimeter to to separate fields within a line by
# - Returns a list of the lines in the specified file
# - Written for general utility purposes 03/15/17 by Shayne Wierbowski
def easyReadLines(fileName, startingDir="", skip=0, comment=None, split=True, sep=""):
	if(startingDir != ""):
		if(startingDir[-1] != "/"):
			startingDir += "/"
	file = open("{0}{1}".format(startingDir, fileName), "r")
	lines = file.read().splitlines()[skip:]
	try:
		lines = [x for x in lines if x[0:len(comment)] != comment]
	except:
		pass
	if(sep):
		lines = [x.split(sep) for x in lines]
	if(not split):
		lines = "\n".join(lines)
	file.close()
	return lines
# FUNCTION END

# - Writes lines to a specified file
# - Accepts as parameters...
#   @fileName - The name for the output file
#   @lines - The lines to write (string, list of strings, or matrix)
#   @file_mode - The mode to open the file as (w, w+, a, a+)
#   @line_sep - Line separator
#   @field_sep - Field Separator
# - Written for general utility purposes 04/17/17 by Shayne Wierbowski
def easyWriteLines(fileName, lines, file_mode="w+", line_sep="\n", field_sep="\t"):
	# Should potentially check that you are not overwritting anything first?
	
	out = open(fileName, file_mode)
	if(lines == []):
		out.close()
		return
	# Lines is a string
	if(type(lines) == type("")):
		out.write(lines)
	# Lines is a list of strings
	elif(type(lines[0]) == type("")):
		out.write(line_sep.join(lines))
	# Lines is a list of fields
	else:
		out.write(line_sep.join([field_sep.join([str(y) for y in x]) for x in lines]))
	out.close()
# FUNCTION END

# Returns a checksum on a file
def hash_file(filename):
	 return hashlib.md5("".join(easyReadLines(filename))).hexdigest()
# FUNCTION END

# - Indicates whether or not an object is iterable
# - Accepts as parameters...
#   @x - The object
#   @exclude - Types to exclude from "iterable"
# - Returns a boolean indicating whether or not x is iterable
# - Written for general utility purposes 05/26/17 by Shayne Wierbowski
def is_iterable(x, exclude=[]):
	if(type(x) in exclude):
		return False
	try:
		_ = iter(x)
		return True
	except KeyboardInterrupt, SystemExit:
		raise
	except:
		return False
# FUNCTION END

# - Flattens a list
# - Accepts as parameters...
#   @l - The list to flatten
#   @exclude - Types to exclude flattening
# - Returns a flattened version of l
# - Written for general utility purposes 05/26/17 by Shayne Wierbowski
def flatten(l, exclude=[str, unicode, np.string_]):
	try:
		while(sum([is_iterable(x, exclude) for x in l]) > 0):
			l = [flatten(item, exclude) if is_iterable(item, exclude) else item for sub2 in [sub if is_iterable(sub, exclude) else [sub] for sub in l] for item in sub2]
		return l
	except KeyboardInterrupt, SystemExit:
		raise
	except:
		return l
# FUNCTION END

# Print Verbosity
verbosity = 0

# Semaphore to Protect vPrint Rights
vPrintSema = threading.BoundedSemaphore()

# - Set new print verbosity threshold
# - Accepts as parameters...
#   @v - The ner verbosity
# - Written for general utility purposes 03/20/17 by Shayne Wierbowski
def setVerbosity(v):
	global verbosity
	global vPrintSema
	
	vPrintSema.acquire()
	verbosity = v
	vPrintSema.release()
# FUNCTION END

# - Verbose Print (for controlled toggleable printing durring debugging)
# - Accepts as parameters...
#   @s - The string to print
#   @priority - Priority of printing (will only print if verbosity >= priority)
# - Written for general utility purposes 03/20/17 by Shayne Wierbowski
def vprint(s, priority=1):
	global verbosity
	global vPrintSema
	
	vPrintSema.acquire()
	if(verbosity >= priority):
		print(s)
	vPrintSema.release()
# FUNCTION END

# Current Working Directory for particular python process
my_cwd = "."
my_last_cwd = my_cwd
my_home = "[REDACTED_PATH]"

# Semaphore to protect CWD (THIS IS EXTREMELY NOT THREAD SAFE NEVER USE MY_CWD IN THREADS)
myCWDSema = threading.BoundedSemaphore()

# - Virtually changes directory
# - Accepts as parameters...
#   @dir - The virtual directory to change to
#   @fromcwd - Boolean flag for whether or not to use virtual cwd as basedir
# - Written for general utility purposes 03/24/17 by Shayne Wierbowski
# - Updated 04/12/17 to add - and ~ functionality
def chdir(dir, fromcwd=True):
	global my_cwd
	global my_last_cwd
	global my_home
	
	myCWDSema.acquire()
	dir.replace("~", my_home)
	if(dir == "-"):
		tmp = my_cwd
		my_cwd = os.path.abspath(my_last_cwd)
		my_last_cwd = tmp
	elif(fromcwd):
		my_last_cwd = my_cwd
		my_cwd = os.path.abspath(my_cwd + "/" + dir)
	else:
		my_last_cwd = my_cwd
		my_cwd = os.path.abspath(dir)
	myCWDSema.release()
# FUNCTION END

# Return absolute path from virtual cwd (unless otherwise specified)
# - Accepts as parameters...
#   @path - The relative path to the file
#   @fromcwd - Boolean flag for whether or not to use virtual cwd as basedir
# - Returns the absolute path for a specified path from virtual cwd
# - Written for general utility purposes 03/24/17 by Shayne Wierbowski
def path(path="", fromcwd=True):
	global my_cwd
	if(fromcwd):
		return os.path.abspath(my_cwd + "/" + path)
	else:
		return os.path.abspath(path)
# FUNCTION END

# Return base filename (blatantly copied from os, but I will remember easier)
# - Accepts as parameters...
#   @path - The path to get the base filename for
#   @Extension - Boolean flag indicating whether or not to retain the extension
# - Returns the basename of a path
# - Written for general utility purposes 03/24/17 by Shayne Wierbowski
def filename(path, Extension=True):
	if(not Extension and "." in path):
		path = path[:path.rfind(".")]
	return path[path.rfind("/") + 1 :]
# FUNCTION END

# List contents of directory (either os.listdir() or glob.glob() re match)
# - Accepts as parameters...
#   @path - The path to get the directory to ls
#   @fromcwd - Boolean flag for whether or not to use virtual cwd as basedir
# - Returns ls results for specified path
# - Written for general utility purposes 03/24/17 by Shayne Wierbowski
def ls(path="", fromcwd=True):
	global my_cwd
	if("*" in path):
		if(fromcwd):
			return [os.path.basename(x) for x in glob.glob(my_cwd + "/" + path)]
		else:
			return [os.path.basename(x) for x in glob.glob(path)]
	else:
		if(fromcwd):
			return os.listdir(my_cwd + "/" + path)
		else:
			return os.listdir(path)
# FUNCTION END

# Performs specified command and returns output as string
# - Accepts as parameters...
#   @cmd - The command
#   @readFrom - Output file to read from
#   @makeTmp - Boolean flag for whether or not output should be written to temp
#   @keepTmp - Boolean flag for whether or not the temp output should be deleted
#   @split - Boolean flag for whether or not the output should be split by line
#   @stdin - Value for the subprocess parameter
#   @stdout - Value for the subprocess parameter
#   @stderr - Value for the subprocess parameter
# - Returns output from the command
# - Written for general utility purposes 03/24/17 by Shayne Wierbowski
def call(cmd, readFrom="", makeTmp=True, keepTmp=False, split=True, stdin=None, stdout=sp.PIPE, stderr=None):
	tmp = ""
	if(makeTmp):
		# NOTE FOR LARGE OUTPUTS YOU MUST SEPCIFY AN OUTPUT FILE
		# "| tee tmp" waits for the command to finish before writting
		# the output --> subprocess freezes
		# This also might mean you can't run any chain of commands
		# involving large intermediary outputs
		tmp = reserveTemp(addRemove=not keepTmp)
		# This should fix the aforementioned problem
		tmp_cmd = cmd.split()
		if(tmp_cmd[len(tmp_cmd)-2] == ">"):
			cmd = cmd + " | tee " + tmp
		elif(tmp_cmd[len(tmp_cmd)-2] == "2>" and ">" in tmp_cmd):
			print "WARNING: I Suck at this so for some reason if you try to redirect stderr after stdout, the contents of stdout don't get redirected properly"
			raise RuntimeError("Shayne Sucks (swap stderr and stdout redirects)")
		else:
			cmd = cmd + " >& " + tmp
		#cmd = (cmd + " | tee " + tmp).split()
	#else:
		#cmd = cmd.split()
	#print cmd
	sp.call(cmd, stdin=stdin, stdout=stdout, stderr=stderr, shell=True)
	
	if(readFrom == ""):
		readFrom = tmp
	if(readFrom != ""):
		lines = easyReadLines(readFrom, split=split)
		if(keepTmp):
			return lines, tmp
		myCleanupFile(tmp)
		return lines
# FUNCTION END

# Obtain the FASTA sequence for an NCBI ID using Entrez
# - Accepts as parameters...
#   @ID - The NCBI ID for the querry FASTA
#   @dest - Destination to write output file
#   @name - Name for output file
#   @email - The Entrez email to use
#   @db - Value for the Entrez.efetch parameter
#   @retmode - Value for the Entrez.efetch parameter
# - Returns the FASTA sequence for the given NCBI ID
# - Written for general utility purposes 03/23/17 by Shayne Wierbowski
def get_Fasta(ID, dest="", name="", email="[REDACTED]", db="protein", retmode="xml"):
	if(name == ""):
		name = ID
	
	html = "http://www.uniprot.org/uniprot/{0}.fasta".format(ID)
	while(True):
		try:
			s = requests.Session()
			text = s.get(html).text
			seq = "".join(text.split("\n")[1:])
			break
		except requests.ConnectionError:
			print "ERROR: Retrying"
	return seq
	#tmp = reserveTemp()
	#call("wget -O {0} http://www.uniprot.org/uniprot/{1}.fasta".format(tmp, ID))
	#seq = "".join(easyReadLines(tmp)[1:])
	#myCleanupFile(tmp)
	#return seq
	
	'''
	Entrez.email = email
	try:
		handle = Entrez.efetch(db=db, id=ID, retmod=retmode)
	except KeyboardInterrupt:
		raise
	except:
		html = "http://www.uniprot.org/uniprot/{0}.fasta".format(ID)
		
		s = requests.Session()
		text = s.get(html).text
		
		seq = str("".join([x for x in text.split("\n") if not ">" in x]))
		if(set(list(seq)).union(set(list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))) == set(list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))):
			return seq
		else:
			print "ERROR: Could not identify FASTA for {0} | {1}".format(ID, set(list(seq)))
			return ""
		
	text = handle.read()
	
	pattern = re.compile("ncbieaa \"([A-Z]+?)\"")
	matches = pattern.finditer(text.replace("\n", ""))
	seqs = []
	for m in matches:
		seqs.append(m.group(1))
	if(len(seqs) > 1):
		print "ERROR: Multiple Sequences in " + ID
	
	if(dest != ""):
		call("mkdir -p {0}".format(dest))
		
		out = open(dest + "/" + name + ".fasta", "w+")
		out.write("> {0} : {1}\n".format(name, ID))
		out.write(seqs[0])
		out.close()
	
	seq = seqs[0].replace("\n", "")
	if(set(list(seq)).union(set(list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))) == set(list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))):
		return seq
	else:
		print "ERROR: Could not identify FASTA for {0} | {1}".format(ID, set(list(seq)))
		return ""
	'''
# FUNCTION END

# Runs a BLAST against PDB Structures and returns relevant structures
# - Accepts as parameters...
#   @fasta - The query sequence
#   @e_cutoff - Cutoff e value
#   @matrix - The scoring matrix
#   @positives_thresh - Threshold for positive percentage
#   @identifies_thresh - Threshold for identify percentage
#   @coCrystal - Flag to select only non-co-cyrstals, only co-crystals, or both (False, None, True)
#   @ligand - Flag to select only non-ligand containing, only ligand-containing, or both (False, None, True)
#   @resolution_thresh - Resolution threshol cutoff
#   @checkFn - A custom filtering function to be applied
#   @sortFn - A custom sorting function to be applied
# - Returns a list of tuples containing (pdbID:chainID, queryStart, queryEnd, positive, identity)
# - Written for general utility purposes 05/26/17 by Shayne Wierbowski
def pdbBLAST(fasta, e_cutoff=10, matrix="BLOSUM62", positives_thresh=0.8, identities_thresh=0.8, coCrystal=None, ligand=None, resolution_thresh=None, checkFn=None, sortFn=None):
	if(positives_thresh > 1):
		positives_thresh /= 100.0
	if(identities_thresh > 1):
		identities_thresh /= 100.0
	
	if(len(set(list(fasta.upper())) - set(list("ACTGN"))) == 0):
		# HANDLE NUCLEOTIDE BLAST
		pass
	
	html = "http://www.rcsb.org/pdb/rest/getBlastPDB1?sequence={0}&eCutOff={1}&matrix={2}".format(fasta.upper(), e_cutoff, matrix)
	
	s = requests.Session()
	text = s.get(html).text
	text = html2text.html2text(text)
	
	pattern = re.compile(".*\|pdbid\|entity\|chain\(s\)\|sequence.*")
	entries = [x for x in pattern.findall(text) if not ">" in x]
	
	def parse(s):
		ID = s[:s.find("\n")]
		values = [tuple(re.split(" = | *", y.replace("  ", " "))[:2]) for y in  [x.strip() for x in re.split(",|\n", s[s.find("\n"):s.find("Query")]) if "=" in x]]
		values = dict([(x[0], float(x[1].split("/")[0])/float(x[1].split("/")[1])) if "/" in x[1] else x for x in values])
		
		QueryPos = [int(y.strip()) for z in [x.split(" ") for x in s[s.find("Query"):].split("\n") if "Query" in x] for y in z if len(set(list(y.strip())) - set(list("1234567890"))) == 0 and y.strip() != ""]
		start = min(QueryPos)
		end = max(QueryPos)
		
		#return (ID, int(values["Length"]), float(values["Score"]), float(values["Expect"]), float(values["Identities"]), float(values["Positives"]), int(start), int(end))
		return (ID, int(end) - int(start) + 1, float(values["Score"]), float(values["Expect"]), float(values["Identities"]), float(values["Positives"]), int(start), int(end))
	# FUNCTION END
	
	# ID Length Score E-value %Ident %Pos Start End
	BLASTs = [parse(x.strip()) for x in text[text.find(entries[-1]) + len(entries[-1]):].split("\n    \n    \n    \n    >") if x.strip() != ""]
	
	BLASTs = [x for x in BLASTs if x[3] <= e_cutoff and x[4] >= identities_thresh and x[5] >= positives_thresh]
	# Sort arbitrarilly prioritized using Length, Positives^2, and Identities
	if(sortFn != None):
		BLASTs = sortFn(BLASTs)
	else:
		BLASTs = sorted(BLASTs, key = lambda x: x[1]*(x[4]**2)*x[5], reverse=True)
	
	def describePDB(pdb):
		s = requests.Session()
		
		resultsDict = dict()
		# PDB Header Description
		html = "http://www.rcsb.org/pdb/rest/describePDB?structureId={0}".format(pdb)
		try:
			text = s.get(html).text
			infoDict = [tuple(str(x.replace("\"", "").replace("\'", "")).split("=")) for x in text.split(" ") if x != "" and "=" in x]
			print "START"
			infoDict = dict([(x[0], flatten(x[1])) if len(x) == 2 else (x[0], flatten(x[1:])) for x in infoDict])
			print "DONE"
			keys = infoDict.keys()
			keepKeys = ["expMethod", "resolution", "deposition_date", "release_date", "last_modification_date", "pubmedId", "pubmedCentralId"]
			for k in keepKeys:
				if(k in keys):
					resultsDict[k] = infoDict[k]
			if(not "resolution" in keys):
				resultsDict["resolution"] = 10
		except:
			resultsDict["resolution"] = 10
		
		# PDB Description
		html = "http://www.rcsb.org/pdb/rest/describeMol?structureId={0}".format(pdb)
		text = s.get(html).text
		try:
			text = s.get(html).text
			pattern = re.compile("<chain.*?/>")
			
			chains = [x[x.find("\"")+1:x.rfind("\"")] for x in pattern.findall(text)]
			resultsDict["chains"] = chains
		except:
			resultsDict["chains"] = []
		
		# Ligand Info
		html = "http://www.rcsb.org/pdb/rest/ligandInfo?structureId={0}".format(pdb)
		text = s.get(html).text
		try:
			text = s.get(html).text
			infoDict = [tuple(str(x.replace("\"", "").replace("\'", "")).split("=")) for x in text.split(" ") if x != "" and "=" in x]
			resultsDict["chemicalID"] = [x[1:] for x in infoDict if x[0] == "chemicalID"]
		except:
			resultsDict["chemicalID"] = []
		
		return resultsDict
	# FUNCTION END
	
	keep = []
	covered = [False]*len(fasta)
	for b in BLASTs:
		# THIS SHOULD BE CALCUATED AHEAD OF TIME SO THAT RESOLUTION CAN BE FACTORED INTO THE SORTING
		b_desc = describePDB(b[0].split("|")[0].split(":")[0])
		print b_desc["chains"]
		print len(b_desc["chains"])
		print len(b_desc["chains"]) > 1
		print coCrystal
		print
		print
		
		if(checkFn != None and not checkFn(b_desc)):
			continue
		if(resolution_thresh != None and float(b_desc["resolution"]) > resolution_thresh):
			continue
		if(coCrystal != None and (len(b_desc["chains"]) > 1) != coCrystal):
			continue
		if(ligand != None and (len(b_desc["chemicalID"]) > 0) != ligand):
			continue
		
		if(sum(covered[b[6]:b[7]+1]) < 0.25*b[1]):
			keep.append([b[0].split("|")[0]] + list(b[1:]) + [b_desc])
			covered[b[6]:b[7]+1] = [True]*b[1]
		if(sum(covered) > 0.9*len(fasta)):
			break
	keep = sorted(keep, key = lambda x: x[6])
	return keep
# FUNCTION END

# Gets the surface residues (by aa position) in the specified PDB or UNIPROT ID
# - Accepts as parameters...
#   @ID - The PDB or UNIPROT ID of the query protein
#   @dest - Destination to write output file
#   @name - Name for output file
#   @chain - The desired chain
#   @type - Either PDB or UNIPROT
#   @sresDB - The database file containing surfae residues
#   @outIndexing - How the output should be indexed (either UNIPROT or PDB) defaults to type
# - Returns a list of amino acid positions that are surface residues or None if ID
#   could not be found
# - Written for general utility purposes 03/24/17 by Shayne Wierbowski
def getSurfRes(ID, dest="", name="", chain="[A-Za-z0-9]*?", type="PDB", sresDB = "[REDACTED_PATH]/ires/parsed_files/sres_perpdb_alltax.txt", outIndexing=""):
	if(name == ""):
		name = ID
	
	if(outIndexing == ""):
		outIndexing = type
	
	# GREP DOESNT SEEM TO RUN TO MANY MATCHES
	# I THINK THIS HAS TO DO WITH SP AND PIPE RUNNING OUT OF MEMORY?
	# LIMITED TO ONLY FIRST MATCH USING -m1
	if(type == "PDB"):
		lines = call("grep -Po -m1 '^{ID}[ \t]*?{chain}[ \t]*?[A-Za-z0-9]*?[ \t]*?[0-9]*?.*$' {DB}".format(ID=ID, chain=chain, DB=sresDB))
	elif(type == "UNIPROT"):
		lines = call("grep -Po -m1 '^[A-Za-z0-9]*?[ \t]*?{chain}[ \t]*?{ID}[ \t]*?[0-9]*?.*$' {DB}".format(ID=ID, chain=chain, DB=sresDB))
	else:
		return
	if(len(lines) > 1):
		print "ERROR: {0} {1}".format(ID, type)
	if(len(lines) == 0):
		return
	else:
		#print lines[0]
		if(outIndexing == "PDB"):
			intervals = [x.split("-") for x in lines[0].split("\t")[5].strip("[]").split(",")]
		elif(outIndexing ==  "UNIPROT"):
			intervals = [x.split("-") for x in lines[0].split("\t")[4].strip("[]").split(",")]
		
		# WEIRD EDGE CASE FOR SOME PDB ENTRIES THAT CONTAIN NEGATIVE NUMBERS AT THE END
		# MAYBE I SHOULD JUST TAKE FROM THE UNIPROT NUMBERING ISNTEAD?
		#print intervals
		intervals = [x if len(x) == 1 or (x[0] != "" and x[1] != "")  else ["-" + x[1]] for x in intervals]
		#print intervals
		numbers = [range(int(x[0]), int(x[1]) + 1) if len(x) > 1 else [int(x[0])] for x in intervals]
		#print numbers
		surfRes = [item for sublist in numbers for item in sublist]
		
		if(dest != ""):
			call("mkdir -p {0}".format(dest))
			
			out = open(dest + "/" + name + "_SurfRes.txt", "w+")
			for s in surfRes:
				out.write("{0}\n".format(s))
			out.close()
		
		return surfRes
# FUNCTION END

# Runs epitope prediction for given sequence
# - Accepts as parameters...
#   @f - File containing peptide sequence to run predictions for
#   @outDir - Directory for output files
#   @HLAs - List of HLA alleles to use
#   @method - Prediction method to use
#   @lens - List of lengths to predict for
#   @SB_thresh - Threshold for strong binders
#   @WB_thresh - Threshold for weak binders
#   @fresh - Boolean flag indicating whether output files (presumably from previous call) 
#            should be regenerated freshly or read / used without regenerating
#   @threads - An indicator for how many threads should be used. Defaults to 1 / calls a single
#              epitope prediction for all alleles. If set > 1, does each allele as a thread (process?)
# - Written for general utility purposes 03/27/17 by Shayne Wierbowski
# - Updated 04/13/17 to run only one netMHC call / parse output into files
#   NOTE - This is actually slightly slower (~7 seconds on a ~22 minute antigen prediction).
#          Should swap to old method and process / thread pool it for better speed. But this is
#          a tricky choice because you could alternatively want to parrallelize the predEpitopes calls
#          for when predicting epitopes on multiple FASTAs
def predEpitopes(f, outDir, HLAs=None, method="netMHC", lens=[8, 9, 10, 11], SB_thresh=0.5, WB_thresh=2.0, fresh=False, threads=1, header=["Pos", "HLA", "Peptide", "Core", "Offset", "I_pos", "I_len", "D_pos", "D_len", "iCore", "Identity", "1-log50k(add)", "Affinity (nM)", "%Rank", "Junk", "Bind_Level"]):
	if(HLAs == None):
		HLAs = call("netMHC -listMHC | grep \"^HLA\"")
	
	call("mkdir -p {0}".format(outDir))
	
	# SINCLE CALL APPROACH
	if(threads == 1):
		outs = []
		if(method == "netMHC"):
			outs.append("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_All.tsv"))
			if(fresh or not fileExists("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_All.tsv"))):
				call("netMHC -a {alleles} -f {input} -l {len} -rth {SB} -rlt {WB} > {outDir}/{input2}".format(alleles=",".join(HLAs), input=f, len=",".join([str(x) for x in lens]), SB=SB_thresh, WB=WB_thresh, outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_All.tsv"))
				lines = easyReadLines("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_All.tsv"))
				
				out = open("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_All.tsv"), "w+")
				out.write("\n".join(re.findall("\n *?([0-9]+.*)", "\n".join(lines))))
				out.close()
				
				# 2019_04_30 - Added better formatting
				fn = "{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_All.tsv")
				pd.read_csv(fn, delim_whitespace=True, names=header[:-1]).drop("Junk", axis=1).to_csv(fn, sep="\t", index=None)
				
			
			# DO GROUP BY IN PANDAS
			#DF = pd.read_csv("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_All.tsv"), delim_whitespace=True, header=None, names=range(16))
			DF = pd.read_csv("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_All.tsv"), sep="\t")
			grouped = DF.groupby("HLA")
			for name, group in grouped:
				outs.append("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + name + ".tsv"))
				group.to_csv("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + name + ".tsv"), sep="\t", index=False)
			
			# HORRENDOUSLY SLOW PARSING BY RE MATCH
			#lines = easyReadLines("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_All.tsv"))
			#for allele in HLAs:
			#	outs.append("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + allele + ".tsv"))
			#	if(fresh or not fileExists("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + allele + ".tsv"))):
			#		out = open("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + allele + ".tsv"), "w+")
			#		out.write("\n".join(re.findall(".*{0}.*".format(allele), "\n".join(lines))))
			#		out.close()
	# MULTIPLE CALL APPROACH
	else:
		outs = []
		if(method == "netMHC"):
			def singlePredEpitopes(f, outDir, allele, lens=[8, 9, 10, 11], SB_thresh=0.5, WB_thresh=2.0, fresh=False):
				outs.append("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + allele + ".tsv"))
				if(fresh or not fileExists("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + allele + ".tsv"))):
					#print("netMHC -a {alleles} -f {input} -l {len} -rth {SB} -rlt {WB} > {outDir}/{input2}".format(alleles=allele, input=f, len=",".join([str(x) for x in lens]), SB=SB_thresh, WB=WB_thresh, outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + allele + ".tsv"))
					call("netMHC -a {alleles} -f {input} -l {len} -rth {SB} -rlt {WB} > {outDir}/{input2}".format(alleles=allele, input=f, len=",".join([str(x) for x in lens]), SB=SB_thresh, WB=WB_thresh, outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + allele + ".tsv"))
					lines = easyReadLines("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + allele + ".tsv"))
					out = open("{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + allele + ".tsv"), "w+")
					out.write("\n".join(re.findall("\n *?([0-9]+.*)", "\n".join(lines))))
					out.close()
					
					# 2019_04_30 - Added better formatting
					fn = "{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + allele + ".tsv")
					pd.read_csv(fn, delim_whitespace=True, names=header[:-1]).drop("Junk", axis=1).to_csv(fn, sep="\t", index=None)
				
				return "{outDir}/{input2}".format(outDir=outDir, input2=filename(f).replace(".fasta", "") + "_prdct_on_" + allele + ".tsv")
			# FUNCTION END
			if(len(HLAs) == len(lens)):
				lens = [lens]*len(HLAs)
			outs = multiThreadFunction(func=singlePredEpitopes, times=len(HLAs), max_threads=threads, as_pool=False, f=f, outDir=outDir, allele=HLAs, lens=lens, SB_thresh=SB_thresh, WB_thresh=WB_thresh, fresh=fresh)
	return outs
# FUNCTION END

# Finds hotspots from a epitope predictions in a given directory
# - Accepts as parameters...
#   @inputDir - Directory containing input epitope prediction files
#   @len - Length of peptide
#   @dest - Destination to write output file
#   @name - Name for output file
#   @WB_thresh - Threshold for weak binders
# - Returns a list of hotspot counts for WB and SB
# - Written for general utility purposes 03/27/17 by Shayne Wierbowski
def findHotspots(inputDir, length, dest="", name=""):
	if(name == ""):
		name = filename(inputDir)
	
	startingDir = path(my_cwd, fromcwd=False)
	chdir(inputDir, fromcwd=False)
	files = [x for x in ls() if os.path.isfile(path(x))]
	
	SBcounts = np.array([0]*length)
	WBcounts = np.array([0]*length)
	for f in files:
		lines = easyReadLines(path(f), split=False)
		SBhits = re.findall(".*<= SB.*", lines)
		WBhits = re.findall(".*<= WB.*", lines)
		for h in WBhits:
			h = h.split()
			WBcounts[int(h[0]):int(h[0])+len(h[2])] += 1
		for h in SBhits:
			h = h.split()
			SBcounts[int(h[0]):int(h[0])+len(h[2])] += 1
	if(dest != ""):
		call("mkdir -p {0}".format(dest))
		
		out = open(dest + "/" + name + "_hotspots.csv", "w+")
		out.write(" ".join([str(x) for x in WBcounts]) + "\n")
		out.write(" ".join([str(x) for x in SBcounts]) + "\n")
		out.close()
	
	chdir(startingDir, fromcwd=False)
	return [list(WBcounts), list(SBcounts)]
# FUNCTION END


# Gets UNIPROT information based on gene name or NCBI ID
# - Accepts as parameters...
#   @gene - Gene name to search using
#   @ncbi - NCBI ID to search using
#   @extra - Boolean flag indicating whether or not to return data for entire entry
#            rather than just the ID
# - Returns information for the top hit in UNIPROT Query
# - Written for general utility purposes 03/27/17 by Shayne Wierbowski
def getUniprotID(gene="", ncbi="", extra=False):
	ncbi = ncbi.replace("_", "+").replace(" ", "+")
	html = "http://www.uniprot.org/uniprot/?query="
	if(gene != "" and ncbi != ""):
		html += "\"{0}\"+AND+gene:{1}&sort=score".format(ncbi, gene)
	elif(gene != ""):
		html += "gene:{0}&sort=score".format(gene)
	elif(ncbi != ""):
		html += "\"{0}\"&sort=score".format(ncbi)
	else:
		return
	
	s = requests.Session()
	text = s.get(html).text
	
	text = re.sub("\n|\t", " ", text)
	
	#pattern = re.compile("</script></td></tr></thead><tbody><tr id=\"(.*?)\" class=\".*?\">")
	#matches = pattern.findall(text)

	pattern = re.compile("</script></td></tr></thead><tbody><tr .*?</tr>")
	matches = pattern.findall(text)
	
	if(len(matches) == 0):
		return
	
	matches[0] = re.sub("</td>", "\n", matches[0])
	matches[0] = re.sub("<script>.*?</script>", "\t", matches[0])
	matches[0] = re.sub("<.*?>", "\t", matches[0])
	
	matches[0] = [[str(z) for z in x.strip().split("\t") if z != ""] for x in matches[0].strip().split("\n") if x != ""]
	
	if(extra):
		return matches[0][0][0], matches[0][1:]
	else:
		return matches[0][0][0]
# FUNCTION END

def fileExists(f):
	return os.path.exists(f) and os.stat(f).st_size != 0
# FUNCTION END

IEDB_Dir = "[REDACTED_PATH]/IEDB_DATA"

# Reads in an Epitope csv as a pandas DF
# - Accepts as parameters...
#   @f - The file for the data
#   @absPath - A boolean indicating whether or not an absolute path is give (if not, default to search in IEDB_Dir)
# - Returns the pandas DF from the given file
# - Written for general utility purposes 04/03/17 by Shayne Wierbowski
def readEpitopeDF(f, absPath=False):
	if(absPath):
		return pd.read_csv(f)
	else:
		return pd.read_csv(IEDB_Dir + "/" + f)
# FUNCTION END

# Writes an Epitope pandas DF to csv
# - Accepts as parameters...
#   @df - The dataframe to write
#   @out - The filename for the output file
#   @absPath - A boolean indicating whether or not an absolute path is give (if not, default to write in IEDB_Dir)
#   @overwrite - A boolean indicating wheterh or not existing files should be overwritten
# - Written for general utility purposes 04/03/17 by Shayne Wierbowski
def writeEpitopeDF(df, out, absPath=False, overwrite=False):
	if(absPath and (overwrite or not fileExists(out))):
		return df.to_csv(out, index=False)
	elif(overwrite or not fileExists(IEBD_Dir + "/" + out)):
		return df.to_csv(IEBD_Dir + "/" + out, index=False)
	else:
		return -1
# FUNCTION END

# Search a pandas DF for entries with a specific value for some attribute
# - Accepts as parameters...
#   @df - The DF to read
#   @col - The column to check value of
#   @value - The value to filter by
#   @outCol - Output Column (if empty, return all DF columns)
# - Returns entries in a pandas DF from a specific column for which another column takes on a specific value
# - Written for general utility purposes 04/03/17 by Shayne Wierbowski
def searchDF(df, col="label", value=True, outCol="Description", asSet=True):
	if(outCol != ""):
		return df[df[col] == value][outCol]
	else:
		return df[df[col] == value]
# FUNCTION END

# Generates a ROC Curve for a given set of labels and scores
# - Accepts as parameters...
#   @labels - The label values
#   @scores - The predicted scores
#   @forcePosAUC - A boolean flag indicating whether scores should be inverted to maximize AUC
#                  (i.e. disallow an AUC < 0.50)
# - Returns the figure object, auc value, false positive rate, true positive rate, and
#   a boolean flad indicating whether or not scores were inverted
# - Written for general utility purposes 04/03/17 by Shayne Wierbowski
def plotROC(labels, scores, forcePosAUC=False, gen_plot=True, figsize=(5, 5), fontsize=12):
	from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
	
	fpr, tpr, _ = roc_curve(labels, scores)
	roc_auc = auc(fpr, tpr)
	
	inverted = False
	if(forcePosAUC and roc_auc < 0.5):
		fpr, tpr, _ = roc_curve(labels, -1*np.array(scores))
		roc_auc = auc(fpr, tpr)
		inverted = True
	
	if(gen_plot):
		fig = plt.figure(figsize=figsize)
		lw = 2
		plt.plot(fpr, tpr, color="darkorange",
		lw=lw, label="ROC curve (area = %0.2f)" % roc_auc)
		plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
		plt.xlim([0.0, 1.0])
		plt.ylim([0.0, 1.05])
		plt.xlabel("False Positive Rate", fontsize=fontsize)
		plt.ylabel("True Positive Rate", fontsize=fontsize)
		if(inverted):
			plt.title("Receiver Operating Characteristic (Inversed Score)", fontsize=fontsize)
		else:
			plt.title("Receiver Operating Characteristic", fontsize=fontsize)
		plt.legend(loc="lower right")
	else:
		fig = None
	
	if(forcePosAUC):
		return fig, roc_auc, fpr, tpr, inverted
	else:
		return fig, roc_auc, fpr, tpr
# FUNCTION END

# Generates a Precision Recall Curve for a given set of labels and scores
# - Accepts as parameters...
#   @labels - The label values
#   @scores - The predicted scores
# - Returns the figure object, auc value, precision, and recall
# - Written for general utility purposes 04/07/17 by Shayne Wierbowski
def plotPRC(labels, scores, gen_plot=True, figsize=(5, 5), fontsize=12):
	from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
	
	# Compute Precision-Recall
	precision, recall, _ = precision_recall_curve(labels, scores)
	average_precision = average_precision_score(labels, scores)
	
	if(gen_plot):
		# Plot Precision-Recall curve
		fig = plt.figure(figsize=figsize)
		lw = 2
		plt.plot(recall, precision, lw=lw, color='navy', label='Precision-Recall curve')
		plt.xlabel('Recall', fontsize=fontsize)
		plt.ylabel('Precision', fontsize=fontsize)
		plt.ylim([0.0, 1.05])
		plt.xlim([0.0, 1.0])
		plt.title('Precision-Recall: AUC={0:0.2f}'.format(average_precision), fontsize=fontsize)
		plt.legend(loc="lower left")
	else:
		fig = None
	
	return fig, average_precision, precision, recall
# FUNCTION END

# Performs a function on a list of arguments using multiple threads
# - Accepts as parameters...
#   @func - The function to call
#   @times - The number of times to call the function (necessary to distinguish a list argument from a list of arguments)
#   @max_threads - The maximum number of threads to run at once (kind of useless since it needs to be explicitly
#                  given anyway (I can't figure out the best way to pass in arguments to this funciton because
#                  you can't give any keyword arguments after non-keyword arguments / if you leave out max_threads
#                  you cant include any unamed parameters for func)
#   @as_pool - Boolean flag indicating whether or not each thread should call its function as a process pool
#   @*args - The arguments to the function
#   @**kwargs - The keywarded arguments to the function
#   Returns a list of outputs from each function call in order
# - Written for general utility purposes 04/04/17 by Shayne Wierbowski
def multiThreadFunction(func, times, max_threads=MAXTHREADS, as_pool=False, *args, **kwargs):
	# WARNING THIS FUNCTION IS COMPLETELY USELESS IF FUNC DOES NOT HAVE AN OS CALL
	# BECAUSE THE GLOBAL INTERPRETER LOCK STOPS PYTHON FROM RUNNING IN PARALLEL
	
	# IF YOU RUN INTO ANY PROBLEMS HERE DONT BE A DUMBASS, JUST USE THREADPOOL INSTEAD
	# JUAN SPECIFICALLY TOLD YOU NOT TO WRITE YOUR OWN FUNCTION, SO YOU HAVE BEEN ADEQUATELY WARNED
	
	# From multiprocessing.pool import ThreadPool
	# pool = ThreadPool(max_threads)
	# output = pool.map(func, args) # args must be a list of tuples eg. [(arg1_c1, arg2_c1, arg3_c1), (arg1_c2, arg2_c2, arg3_c2), ...]
	
	# Get arguments for the function
	func_args = inspect.getargspec(func)
	
	# Create Dictionary for argumnets to the function
	func_params = dict([(x, None) for x in func_args[0]])
	
	# Fill in Default Values
	try:
		for i in range(len(func_args[3])):
			func_params[func_args[0][len(func_args[0]) - len(func_args[3]) + i]] = [func_args[3][i]]*times
	except TypeError:
		pass
	# Fill in Provided Values
	func_kwargs = dict()
	for key, value in kwargs.iteritems():
		# If Provided Param is an argument to the function
		# Add it to the dictionary
		if(key in func_args[0]):
			func_params[key] = value
		# If Provided Param is not an argument to the function, but function has **kwargs
		# Add it to the func_kwargs dictionary
		else:
			func_kwargs[key] = value
	
	# If the function being called does not accept **kwargs, set it to None (should probably throw an error)
	if(func_args[2] == None):
		func_kwargs = dict()
	
	# If the function being called does not accept *args, set it to None, else set it to *args
	if(func_args[1] == None):
		func_args = []
	else:
		func_args = list(args)
	
	# Create Thread Class to handle function
	class myThread (threading.Thread):
		def __init__(self, i, func, as_pool, params, args, kwargs):
			threading.Thread.__init__(self)
			self.i = i
			self.func = func
			self.as_pool = as_pool
			self.params = params
			self.args = args
			self.kwargs = kwargs
		# FUNCTION END
		
		def run(self):
			self.args = [self.params[key] for key in inspect.getargspec(self.func)[0]] + self.args
			if(self.as_pool):
				pool = Pool(1)
				# DONT KNOW WHAT TO DO WITH KWARGS HERE
				outputs[self.i] = pool.map(self.func, self.args)[0]
				pool.close()
				pool.join()
			else:
				outputs[self.i] = self.func(*self.args, **self.kwargs)
			'''
			if(self.args == None):
				if(selg.kwargs == None):
					output[self.i] = func(**params)
				else:
					output[self.i] = func(**params, **kwargs)
			'''
			sema.release()
		# FUNCTION END
	# CLASS END
	
	# Call Function as individual threads
	sema = threading.BoundedSemaphore(max_threads)
	threads = [None]*times
	outputs = [None]*times
	for i in range(times):
		# Prepare Thead params
		thread_params = copy.deepcopy(func_params)
		for key, value in thread_params.iteritems():
			try:
				if(len(value) == times):
					thread_params[key] = value[i]
				else:
					thread_params[key] = value
			except TypeError:
				thread_params[key] = value
		# Prepare Thead args
		thread_args = copy.deepcopy(func_args)
		for j in range(len(thread_args)):
			try:
				if(len(thread_args[j]) == times):
					thread_args[j] = thread_args[j][i]
				else:
					thread_args[j] = thread_args[j]
			except TypeError:
				thread_args[j] = thread_args[j]
		# Prepare Thead kwargs
		thread_kwargs = copy.deepcopy(func_kwargs)
		for key, value in thread_kwargs.iteritems():
			try:
				if(len(value) == times):
					thread_kwargs[key] = value[i]
				else:
					thread_kwargs[key] = value
			except TypeError:
				thread_kwargs[key] = value
		
		# Acquire Semaphore
		sema.acquire()
		threads[i] = myThread(i, func, as_pool, thread_params, thread_args, thread_kwargs)
		threads[i].start()
	
	# Join the threads and get any output they may have returned
	for t in threads:
		t.join()
	
	# Return outputs
	return outputs
# FUNCTION END

# Generates a plot for a given antigen comparing known vs predicted epitope hotspots (also runs predictions)
# - Accepts as parameters...
#   @antigen_ID - The ID for the antigen (read from provided antigen map)
#   @outDir - Location to store outputs
#   @antigen_map - Antigen Map file to read information from (known epitopes / antigen sequence)
#   @method - The parameter for predEpitopes
#   @lens - The parameter for predEpitopes
#   @SB_thresh - The parameter for predEpitopes
#   @WB_thresh - The parameter for predEpitopes
#   @fresh - The parameter for predEpitopes
#   @threads - The parameter for predEpitopes
#   Returns an epitope plot for an antigen along with hotspot values (WB and SB predictions) and actual epitope values
# - Written for general utility purposes 04/13/17 by Shayne Wierbowski
def plotAntigenPred(antigen_ID, outDir, antigen_map="[REDACTED_PATH]/Surface_NetMHC_Perf/AntigenMap_All", method="netMHC", lens=[8, 9, 10, 11], SB_thresh=0.5, WB_thresh=2.0, fresh=False, threads=1):
	call("mkdir -p {0}".format(outDir))
	
	# Handle list of antigens
	# Should actually modify this to run predictions on the whole list of antigens
	# and parse the outputs rather than calling recursively
	if(type(antigen_ID) != type("")):
		try:
			return [plotAntigenPred(a, outDir + "/" + a, antigen_map, method, lens, SB_thresh, WB_thresh, fresh, threads) for a in antigen_ID]
		except:
			return
	
	# Get List of Actual Epitopes / Antigen Sequence (from Antigen Map)
	AntMapDB = pd.read_csv(antigen_map, header=None, sep="\t")
	
	Epitopes = AntMapDB[AntMapDB[0] == antigen_ID].reset_index()
	
	antigen_Fasta = Epitopes[3].unique()[0]
	out = open("{0}/{1}.fasta".format(outDir, antigen_ID), "w+")
	out.write("> {0}\n".format(antigen_ID))
	out.write(antigen_Fasta)
	out.close()
	
	Epitopes["starts"] = [[m.start() for m in re.finditer(Epitopes[1][i], Epitopes[3][i])] for i in range(len(Epitopes))]
	Epitopes["ends"] = [[j + len(Epitopes[1][i]) for j in Epitopes["starts"][i]] for i in range(len(Epitopes))]
	
	Ecounts = np.array([0]*len(antigen_Fasta))
	for i in range(len(Epitopes)):
		for j in range(len(Epitopes["starts"][i])):
			Ecounts[Epitopes["starts"][i][j]:Epitopes["ends"][i][j]] += 1
	
	# Run Predictions
	outs = predEpitopes(f="{0}.fasta".format(antigen_ID), outDir=outDir + "/Predictions", method=method, lens=lens, SB_thresh=SB_thresh, WB_thresh=WB_thresh, fresh=fresh, threads=threads)
	
	# Analyze Predictions / Find Hotspots
	if(fresh or not fileExists("{outDir}/{ID}_hotspots.csv".format(outDir=outDir, ID=antigen_ID))):
		chdir(outDir + "/Predictions")
		files = [x for x in ls() if os.path.isfile(path(x))]
		SBcounts = np.array([0]*len(antigen_Fasta))
		WBcounts = np.array([0]*len(antigen_Fasta))
		total = 0
		for f in files:
			#print path(f)
			lines = easyReadLines(path(f), split=False)
			SBhits = re.findall(".*<= SB.*", lines)
			WBhits = re.findall(".*<= WB.*", lines)
			for h in WBhits:
				h = h.split()
				WBcounts[int(h[0]):int(h[0])+len(h[2])] += 1
			for h in SBhits:
				h = h.split()
				SBcounts[int(h[0]):int(h[0])+len(h[2])] += 1
		out = open(path("../{0}_hotspots.csv".format(antigen_ID)), "w+")
		out.write(" ".join([str(x) for x in WBcounts]) + "\n")
		out.write(" ".join([str(x) for x in SBcounts]) + "\n")
		out.close()
		chdir("-")
	else:
		lines = easyReadLines("{outDir}/{ID}_hotspots.csv".format(outDir=outDir, ID=antigen_ID))
		WBcounts = np.array(lines[0].split(), int)
		SBcounts = np.array(lines[1].split(), int)
	
	# Get Surface Residues
	surfRes = getSurfRes(antigen_ID, dest=outDir, type="UNIPROT")
	
	# Handle Plotting
	fig = plt.figure()
	
	# ADD NORMALIZATION?
	
	WBcounts = WBcounts.astype("float") / max([max(WBcounts), max(SBcounts)])
	SBcounts = SBcounts.astype("float") / max([max(WBcounts), max(SBcounts)])
	Ecounts = Ecounts.astype("float") / max(Ecounts)
	
	plt.plot(range(1, len(WBcounts) + 1), WBcounts)
	plt.plot(range(1, len(SBcounts) + 1), SBcounts)
	plt.plot(range(1, len(Ecounts) + 1), Ecounts)
	
	if(surfRes):
		for s in surfRes:
			plt.axvspan(s, s+1, facecolor="g", alpha=0.25)
	return fig, WBcounts, SBcounts, Ecounts, surfRes, Epitopes
# FUNCTION END

# Converts raw second time to Hour:Min:Second time
# - Accepts as parameters...
#   @time - time in raw seconds
#   Returns a human readable time string
# - Written for general utility purposes 04/17/17 by Shayne Wierbowski
def HRTime(time):
	return "{0:02.0f}:{1:02.0f}:{2:02g}".format(int(math.floor(float(time) / 60 / 60)), int(math.floor(float(time) / 60 % 60)), float(time)%60)
# FUNCTION END

start_T = 0
end_T = 0
def Tstart():
	global start_T
	
	start_T = time.time()
# FUNCTION END

def Tend(message="", reset=True, ret=False):
	global start_T
	global end_T
	
	end_T = time.time()
	
	net_T = end_T - start_T
	r = HRTime(net_T)
	
	if(reset):
		start_T = end_T
	if(ret):
		return net_T
	print message, r
# FUNCTION END

# Converts Hour:Min:Second time to raw second time
# - Accepts as parameters...
#   @time - human readable time string
#   Returns time in raw seconds
# - Written for general utility purposes 04/17/17 by Shayne Wierbowski
def rawSecTime(time):
	time = time.split(":")
	return 60*60*float(time[0]) + 60*float(time[1]) + float(time[2])
# FUNCTION END

# Times a function
# - Accepts as parameters...
#   @func - The function to call
#   @HR - Boolean flag for returning human readable times or raw seconds
#   @*args - The arguments to the function
#   @**kwargs - The keywarded arguments to the function
#   Returns a list of outputs from each function call in order
# - Written for general utility purposes 04/04/17 by Shayne Wierbowski
def funcTime(func, HR=True, *args, **kwargs):
	start = time.time()
	func(*args, **kwargs)
	end = time.time()
	if(HR):
		return HRTime(end - start)
	else:
		return end - start
# FUNCTION END

# Reads a Network from a file and returns it as a list of tuples
# - Written for midterm of BTRY 4381 03/27/17 by Shayne Wierbowski
def readNetwork(file, skip=0, comment=None, cols=[0, 1], directed=False):
	network = set(["\t".join([x.split()[cols[0]], x.split()[cols[1]]]) if directed else "\t".join(sorted([x.split()[cols[0]], x.split()[cols[1]]])) for x in [y for y in open(file, "r").read().splitlines()[skip:] if y[0] != comment]])
	return sorted([(x.split()[0], x.split()[1]) for x in network])
# FUNCTION END

# Converts a network to a dictionary with keys value pairs corresponding to each node and its connected nodes
# - Written for midterm of BTRY 4381 03/27/17 by Shayne Wierbowski
def net2Dict(network, directed=False):
	if(directed):
		netDict = {}
		[netDict.setdefault(x, []) for x in set([y[0] for y in network])]
		[netDict[x[0]].append(x[1]) for x in network]
	else:
		netDict = {}
		[netDict.setdefault(x, []) for x in set([y[0] for y in network]).union(set(z[1] for z in network))]
		[(netDict[x[0]].append(x[1]), netDict[x[1]].append(x[0])) for x in network]
	return dict([(key, list(set(netDict[key]))) for key in netDict])
# FUNCTION END

# Converts a fasta file to a fastq file
# - Written for midterm of BTRY 4381 03/27/17 by Shayne Wierbowski
def fastq2fasta(fasta):
	lines = open(fasta, "r").read().splitlines()
	out = open(fasta.replace(".fastq", ".fasta"), "w+")
	out.write("\n".join([lines[i].replace("@", ">") for i in range(len(lines)) if i % 4 == 0 or i % 4 == 1]))
	out.close()
	return fasta.replace(".fastq", ".fasta")
# FUNCTION END

# Performs a BLAST
# - Written for midterm of BTRY 4381 03/27/17 by Shayne Wierbowski
def runBlast(query, db, dbtype="nucl", outfmt=6, max_target_seqs=1000000000, fresh=False, makedbouts=[".nhr", ".nin", ".nsq"]):
	if(type(query) != type([])):
		query = [query]
	if(type(db) != type([])):
		db = [db]
	outs = []
	for d in db:
		if(fresh or False in [fileExists(d + x) for x in makedbouts]):
			call("makeblastdb -in {0} -dbtype {1}".format(d, dbtype))
		for q in query:
			nextout = "blast_{0}_to_{1}".format(filename(q)[:filename(q).find(".")], filename(d)[:filename(d).find(".")])
			outs.append(nextout)
			if(fresh or not fileExists(nextout)):
				call("blast{0} -query {1} -db {2} -out {3} -outfmt {4} -max_target_seqs {5}".format(dbtype[0], q, d, nextout, outfmt, max_target_seqs))
	return outs
# FUNCTION END

# Reads a BLAST output as pandas df
# Written for general utility purposes 01/4/18 by Shayne Wierbowski
def blast2df(file, header=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]):
	return pd.read_csv(file, sep="\t", header=None, names=header)
# FUNCTION END

# Reads a SAM output as pandas df
# Written for general utility purposes 01/4/18 by Shayne Wierbowski
def sam2df(file):
	'''
	lines = easyReadLines(file)
	lines = [x.split("\t") for x in lines if not x[0] == "@"]
	df = pd.DataFrame(lines)
	df = df.rename(index=str, columns={x[1]:x[0] for x in zip(["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"], range(11))})
	if(len(df) == 0):
		return pd.DataFrame([], columns=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"])
	return df
	'''
	with open(file, "r") as f:
		lines = []
		for l in f:
			if(l[0] == "@"):
				continue
			lines.append(l.strip().split("\t"))
		
		df = pd.DataFrame(lines)
		df.rename(index=str, inplace=True, copy=False, columns={x[1]:x[0] for x in zip(["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"], range(11))})
		if(len(df) == 0):
			return pd.DataFrame([], columns=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"])
		df["FLAG"] = df["FLAG"].astype(int)
		df["POS"] = df["POS"].astype(int)
		df["MAPQ"] = df["MAPQ"].astype(int)
		
		#while(len(lines)):
		#	lines.pop()
		del lines
		
		return df
# FUNCTION END

# Performs a BWA Alignment
# - Written for midterm of BTRY 4381 03/27/17 by Shayne Wierbowski
def runBWA(query, db, threads=12, fresh=False, indexouts=[".amb", ".ann", ".bwt", ".pac", ".sa"]):
	if(type(query) != type([])):
		query = [query]
	if(type(db) != type([])):
		db = [db]
	outs = []
	for d in db:
		#if(fresh or False in [fileExists(d + x) for x in indexouts]):
		if(False in [fileExists(d + x) for x in indexouts]):
			call("bwa index {0}".format(d))
		for q in query:
			nextout = "BWA_align_{0}_to_{1}".format(filename(q)[:filename(q).find(".")], filename(d)[:filename(d).find(".")])
			outs.append(nextout)
			if(fresh or not fileExists(nextout)):
				call("bwa mem -a -t {0} {1} {2} > {3}".format(threads, d, q, nextout))
				#call("bwa bwasw -t {0} {1} {2} > {3}".format(threads, d, q, nextout))
	return outs
# FUNCTION END

# Does things
# - Accepts as parameters...
#   @login - login
#   Returns a stuff
# - Written for general utility purposes 05/03/17 by Shayne Wierbowski
def retrieve(login):
	if(authentication.Master_Logged):
		return authentication.Master_recall_login(login)
	else:
		authentication.Master_Login(foo, bar)
		return authentication.Master_recall_login(login)
# FUNCTION END

# Logs into MYSQL
# - Accepts as parameters...
#   @host - the host to connect to
#   @db - the database to connect to
#   Returns a MYQSLdb cursor
# - Written for general utility purposes 05/03/17 by Shayne Wierbowski
#def MYSQL(host="localhost", db=""):
def MYSQL(host="[REDACTED]", db=""):
	if(not foo == "nice"):
		c = retrieve("MYSQL")
	else:
		c = [user, getpass.getpass("Enter MySQL Password:")]
	mydb = MySQLdb.connect(host=host, user=c[0], passwd=c[1], db=db)
	#mydb = MySQLdb.connect(host=host, user=c[0], passwd=c[1], db=db)
	return mydb.cursor()
# FUNCTION END

# Converts a MYSQL selection into a file
# - Accepts as parameters...
#   @db - the MYSQL database to use
#   @filename - the output filename for the selection
#   @selection - the selection to apply (by default, select *)
#   @table - the table to select from (only necessary for default selection)
# - Written for general utility purposes 09/28/17 by Shayne Wierbowski
def SQLTable2File(db, filename, selection = "select * from {0};", table=""):
	c = MYSQL(db=db)
	c.execute(selection.format(table))
	header = [x[0] for x in c.description]
	lines = c.fetchall()
	easyWriteLines(filename, [header] + list(lines))
	c.close()
# FUNCTION END

# Runs a specified .sql script
# - Accepts as parameters...
#   @script - the .sql script
#   @db - the database in which to run the .sql script
#   @asPandasDF - boolean flag indicating whether or not to try to return output
#                 as a pandas DataFrame (probably buggy)
#   Returns the output from the .sql script
# - Written for general utility purposes 06/19/17 by Shayne Wierbowski
def runMYSQLScript(script, db="", asPandasDF=False):
	if(not foo == "nice"):
		c = retrieve("MYSQL")
	else:
		c = [foo, bar]
	lines = call("mysql -u {0} --password={1} {2} < {3}".format(c[0], c[1], db, script))[1:]
	if(asPandasDF):
		try:
			return pd.DataFrame([x.split("\t") for x in lines[1:]], columns=lines[0].split("\t"))
		except:
			return lines
	else:
		return lines
# FUNCTION END

# Converts a text file to a MySQL table
# - Accepts as parameters...
#   @filename - the file to convert
#   @db - the database in which to run create the table
#   @columnType - a list of MySQL types for each column (NOTE can include text marking
#                 e.g. keys, index, not null)
#   @indices - a list of any indices to be created from column names
#   @primary_key - the column(s) to use as the primary key
#   @foreign_keys - the columns to use as foreign keys
#   @extraneousCreate - any extraneous text to include in the end of CREATE TABLE command
#                       e.g. INDEX, PRIMARY KEY
#   field_sep - the field separator
#   line_sep - the line separator
#   header - boolean flag indicating whether or not the file contains a header
#   skip - a specified number of lines to skip from the beginning (applied before header,
#          after comments removed)
#   comment - a symbol indicating that a line is a comment and should be ignored
#   names - a list of names to use for the columns
#   update - a boolean indicating whether or not the table should be updated (note: this just adds
#            all of the new enties into the existing table (WILL add duplicates, NOT a true update)
#   replace - a boolean indicating whether or not the table should be dropped and replaced
#   none2Null - a boolean indicating whether or not plain text "None" should be changed to "NULL"
#   foreign_key_on_delete - indicates what to do when the data a foreign key refers to is dropped (options, "NO ACTION", "Set NULL", "Set DEFAULT", "CASCADE")
# - Written for general utility purposes 06/19/17 by Shayne Wierbowski
def txt2MYSQL(fileName, db, columnTypes, indices=None, primary_key=None, foreign_keys=None, extraneousCreate="", field_sep="\t", line_sep="\n", header=True, skip=0, comment=None, names=None, update=False, replace=False, none2Null=True, foreign_key_on_delete="Set NULL"):
	columnTypes = [x.upper() for x in columnTypes]
	
	f = open("[REDACTED_PATH]/txt2MYSQL_Calls.txt", "a+")
	f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(fileName, db, columnTypes, indices, primary_key, foreign_keys))
	f.close()
	
	if(extraneousCreate != "" and extraneousCreate[:2] != ", "):
		extraneousCreate = ", " + extraneousCreate
	
	# Check if table exists
	c = MYSQL(db=db)
	c.execute("show tables;")
	table_exists = filename(fileName, Extension=False) in flatten(c.fetchall())
	
	if(table_exists and update == False and replace == False):
		print "ERROR: Table already exists; Exiting; To replace or update the table, set replace=Ture or update=True"
		return
	
	# Temporary .sql file
	tmp = reserveTemp()
	
	# Temporary file to write updated txt
	tmp2 = reserveTemp()
	
	# Read txt (ignoring commens when applicable)
	if(comment != None):
		flines = [[y.replace("\t", "\\t") for y in x.split(field_sep)] for x in "\n".join(easyReadLines(fileName)).split(line_sep) if x[0:len(comment)] != comment]
	else:
		flines = [[y.replace("\t", "\\t") for y in x.split(field_sep)] for x in "\n".join(easyReadLines(fileName)).split(line_sep)]
	
	# Skip lines
	flines = flines[skip:]
	
	# Remove header
	if(header):
		head = flines[0]
		flines = flines[1:]
		if(names == None):
			names = head
	
	# MANUALLY FIX BOOLEAN TYPES
	bools = [i for i in range(len(columnTypes)) if columnTypes[i] == "BOOL"]
	for b in bools:
		def do(x):
			try:
				return[str(int(int(x) == 1))]
			except KeyboardInterrupt:
				raise
			except:
				pass
			try:
				return [str(["FALSE", "TRUE"].index(x.upper()))]
			except KeyboardInterrupt:
				raise
			except:
				return [str(0)]
		# FUNCTION END
		flines = [x[:b] + do(x[b]) + x[b+1:] for x in flines]
	columnTypes = [x if x != "BOOL" else "BOOLEAN" for x in columnTypes]
	
	# Manually fix None (When I load NULL value in SQL into python and save it, it writes as "None".
	# I include a manual fix to change this back. This assumes that I never want a column to say "None")
	flines = [[x if not x.upper() in ["NONE", "", "NAN", "NULL"] else "\N" for x in y] for y in flines]
	# Write updated txt
	print tmp
	print tmp2
	print len(flines)
	print len(flines[0])
	print flines[0]
	easyWriteLines(tmp2, flines, line_sep=line_sep)
	
	def do(x):
		if(foreign_keys == None):
			return ""
		elif(True in [x in str(y[0]) for y in foreign_keys]):
			return [y[1] for y in foreign_keys if x in str(y[0])][0]
		return ""
	# FUNCTION END
	
	a = pd.DataFrame([[do(names[i]), names[i], columnTypes[i], flines[0][i]] for i in range(len(names))], columns=["Ref", "Field", "Type", "Example"])
	
	# Manually convert TEXT to VARCHAR(N)
	data = pd.read_csv(tmp2, sep="\t", names=names)
	def do(x, y):
		if(x != "TEXT"):
			return x
		else:
			lens = data[y].map(lambda a: len(str(a)))
			max_len = max(lens)
			avg_diff = np.mean([a - max_len for a in lens])
			print y, x
			print max_len, avg_diff
			if(max_len > sys.maxint):
				return x
			if(abs(avg_diff) <= 2):
				return "CHAR({0})".format(max_len)
			else:
				return "VARCHAR({0})".format(max_len)
	# FUNCTION END
	print len(names)
	print len(columnTypes)
	columnTypes = [do(columnTypes[i], names[i]) for i in range(len(columnTypes))]
	print columnTypes
	
	# Sanitize Names
	def sanitize(x):
		x = "`{0}`".format(x)
		
		# Sanitize Equality Signs
		#x = x.replace("<=", "Less Than_Equal")
		#x = x.replace("<", "Less Than")
		#x = x.replace(">=", "Greater Than_Equal")
		#x = x.replace(">", "Greater Than")
		#x = x.replace("=", "Equal")
		
		# Sanitize Quotes
		#x = x.replace("\"", "")
		#x = x.replace("\'", "")
		
		#x = x.replace("#", "Num")
		
		x = x.replace(" ", "_").replace("(", "").replace(")", "")
		x = x.replace(".", "_")
		return x
	# FUNCTION END
	names = [sanitize(x) for x in names]
	
	# Create .sql
	lines = "\n".join(easyReadLines("[REDACTED_PATH]/bin/MySQL_New_Table.template"))
	lines = lines.replace("^[TABLENAME]^", filename(fileName, Extension=False))
	lines = lines.replace("^[COLUMNDESCRIPTIONS]^", ", ".join([" ".join(list(x) + ["NULL"]) for x in zip([x for x in names], columnTypes)]))
	lines = lines.replace("^[EXTRANEOUS]^", extraneousCreate)
	lines = lines.replace("^[INPUTFILE]^", tmp2)
	lines = lines.replace("^[LINESEP]^", line_sep)
	
	constraints_dropped = []
	if(not table_exists):
		# Remove entire Drop line
		lines = re.sub("\^\[DROP_FLAG\]\^.*", "", lines)
		# Remove Create flag, keep line
		lines = re.sub("\^\[CREATE_FLAG\]\^", "", lines)
	elif(replace == True):
		# Drop all references to this table
		constraints_dropped = dropAllInReferences(filename(fileName, Extension=False), db)
		# Drop all references from this table
		dropAllOutReferences(filename(fileName, Extension=False), db)
		
		# Remove Drop flag, keep line
		lines = re.sub("\^\[DROP_FLAG\]\^", "", lines)
		# Remove Create flag, keep line
		lines = re.sub("\^\[CREATE_FLAG\]\^", "", lines)
	elif(update == True):
		# Remove entire Drop line
		lines = re.sub("\^\[DROP_FLAG\]\^.*", "", lines)
		# Remove entire Create line
		lines = re.sub("\^\[CREATE_FLAG\]\^.*", "", lines)
	else:
		print "ERROR; Existing"
		return
	
	# HANDLE INDICES
	if(not indices == None):
		# Prepare Index lines
		def auto_name(x):
			if(x[0] == ""):
				x[0] = "ind_" + "_".join([sanitize(y) for y in x[1]])
				x[0] = "`{0}`".format(x[0].replace("`", ""))
			return x
		# FUNCTION END
		indices = [auto_name(x) for x in indices]
		index_lines = "\n".join(["ALTER TABLE {0} ADD INDEX {1} ({2});".format(filename(fileName, Extension=False), x[0], ", ".join([sanitize(y) for y in x[1]])) for x in indices])
		
		# Update Index lines
		lines = re.sub("\^\[INDEX_FLAG\]\^", index_lines, lines)
	else:
		# Remove entire Index line
		lines = re.sub("\^\[INDEX_FLAG\]\^.*", "", lines)
	
	# HANDLE PRIMARY KEY
	if(not primary_key == None):
		# Prepare Primary Key lines
		pkey_lines = "\n".join(["ALTER TABLE {0} ADD PRIMARY KEY ({1});".format(filename(fileName, Extension=False), ", ".join([sanitize(y) for y in x])) for x in [primary_key]])
		
		# Update Primary Key lines
		lines = re.sub("\^\[PRIMARY_KEY_FLAG\]\^", pkey_lines, lines)
	else:
		# Remove entire Primary Key line
		lines = re.sub("\^\[PRIMARY_KEY_FLAG\]\^.*", "", lines)
	
	# HANDLE FOREIGN KEYS
	if(not foreign_keys == None):
		# Prepare Foreign Keys lines
		fkeys_lines = "\n".join(["ALTER TABLE {0} ADD FOREIGN KEY ({1}) REFERENCES {2}({3}) ON DELETE {4};".format(filename(fileName, Extension=False), ", ".join([sanitize(y) for y in x[0]]), x[1], ", ".join([sanitize(y) for y in x[2]]), foreign_key_on_delete) for x in foreign_keys])
		
		# Update Foreign Keys lines
		lines = re.sub("\^\[FOREIGN_KEYS_FLAG\]\^", fkeys_lines, lines)
	else:
		# Remove entire Foreign Keys line
		lines = re.sub("\^\[FOREIGN_KEYS_FLAG\]\^.*", "", lines)
	
	easyWriteLines(tmp, lines, line_sep=line_sep)
	
	r = runMYSQLScript(tmp, db)
	
	# Check number of rows in table
	c.execute("select * from {0};".format(filename(fileName, Extension=False)))
	table_len = len(c.fetchall())
	if(not table_len == len(flines)):
		print "WARNING: SOME ROWS HAVE BEEN DROPPED"
		print table_len
		print len(flines)
	else:
		print "SUCCESS!"
	
	# Recreate all references to this table if table had already existed
	for fk in constraints_dropped:
		print fk
		print("ALTER TABLE {0} ADD {1};".format(fk[0], fk[1].strip(",")))
		c.execute("ALTER TABLE {0} ADD {1};".format(fk[0], fk[1].strip(",")))
	
	#myCleanupFile(tmp)
	#myCleanupFile(tmp2)
	if(len(flines) < 30):
		b = pd.DataFrame(flines, columns=[x.replace("`", "") for x in names])
	else:
		b = pd.DataFrame(flines[:15] + flines[-15:], columns=[x.replace("`", "") for x in names])
	c.close()
	print "{0} Lines".format(len(flines))
	return r, a, b
# FUNCTION END

# Adds a foreign key to an existing table
# - Accepts as parameters...
#   @table - the table to modify
#   @db - the database the table is in
#   @foreign_key - the foreign key
#   foreign_key_on_delete - indicates what to do when the data a foreign key refers to is dropped (options, "NO ACTION", "Set NULL", "Set DEFAULT", "CASCADE")
# - Written for general utility purposes 04/03/18 by Shayne Wierbowski
def addForeignKey(table, db, foreign_key, foreign_key_on_delete="Set NULL"):
	# Check if table exists
	c = MYSQL(db=db)
	c.execute("show tables;")
	table_exists = table in flatten(c.fetchall())
	
	if(not table_exists):
		print "ERROR: Table does not exist"
		return
	
	# Check number of rows in table
	c.execute("select * from {0};".format(table))
	table_len_orig = len(c.fetchall())
	
	# Temporary .sql file
	tmp = reserveTemp()
	
	# Sanitize Names
	def sanitize(x):
		x = "`{0}`".format(x)
		
		# Sanitize Equality Signs
		#x = x.replace("<=", "Less Than_Equal")
		#x = x.replace("<", "Less Than")
		#x = x.replace(">=", "Greater Than_Equal")
		#x = x.replace(">", "Greater Than")
		#x = x.replace("=", "Equal")
		
		# Sanitize Quotes
		#x = x.replace("\"", "")
		#x = x.replace("\'", "")
		
		#x = x.replace("#", "Num")
		
		x = x.replace(" ", "_").replace("(", "").replace(")", "")
		return x
	# FUNCTION END
	
	# HANDLE FOREIGN KEYS
	# Prepare Foreign Keys lines
	lines = "\n".join(["ALTER TABLE {0} ADD FOREIGN KEY ({1}) REFERENCES {2}({3}) ON DELETE {4};".format(table, ", ".join([sanitize(y) for y in x[0]]), x[1], ", ".join([sanitize(y) for y in x[2]]), foreign_key_on_delete) for x in [foreign_key]])
	
	easyWriteLines(tmp, lines)
	print tmp
	
	r = runMYSQLScript(tmp, db)
	
	# Check number of rows in table
	c.execute("select * from {0};".format(table))
	table_len_after = len(c.fetchall())
	if(not table_len_orig == table_len_after):
		print "WARNING: SOME ROWS HAVE BEEN DROPPED"
	else:
		print "SUCCESS!"
	
	#myCleanupFile(tmp)
	c.close()
	return r
# FUNCTION END

# Drops all references to a specific MySQL Table in other tables (currently only in same database)
# - Accepts as parameters...
#   @table - the table to drop references to
#   @db - the database the table is in
# - Written for general utility purposes 04/04/18 by Shayne Wierbowski
def dropAllInReferences(table, db):
	# Get all Tables
	c = MYSQL(db=db)
	c.execute("SHOW TABLES;")
	tables = flatten(c.fetchall())
	
	fks_all = []
	constraints_dropped = []
	for t in tables:
		# Get CREATE TABLE Command
		c.execute("SHOW CREATE TABLE {0};".format(t))
		r = c.fetchall()
		
		# Get all Foreign Keys
		fks = [x for x in r[0][1].split("\n") if "FOREIGN KEY" in x]
		
		# Get all Foreign Keys Referencing Table
		fks = [x for x in fks if "REFERENCES `" + table in x]
		
		constraints_dropped += [(t, x) for x in fks]
		
		# Get names of all Foreign Keys
		fks = [x.split("FOREIGN KEY")[0].replace("CONSTRAINT", "").strip() for x in fks]
		
		fks_all += [(t, x) for x in fks]
	fks_all = list(set(fks_all))
	constraints_dropped = list(set(constraints_dropped))
	
	for fk in fks_all:
		#print fk
		c.execute("ALTER TABLE {0} DROP FOREIGN KEY {1};".format(fk[0], fk[1]))
	
	c.close()
	return constraints_dropped
# FUNCTION END

# Drops all references from a specific MySQL Table to other tables
# - Accepts as parameters...
#   @table - the table to drop references from
#   @db - the database the table is in
# - Written for general utility purposes 04/04/18 by Shayne Wierbowski
def dropAllOutReferences(table, db):
	# Get all Tables
	c = MYSQL(db=db)
	c.execute("SHOW TABLES;")
	tables = flatten(c.fetchall())
	
	if(not table in tables):
		print "ERROR: Table does not exist"
		return
	
	fks_all = []
	constraints_dropped = []
	for t in [table]:
		# Get CREATE TABLE Command
		c.execute("SHOW CREATE TABLE {0};".format(t))
		r = c.fetchall()
		
		# Get all Foreign Keys
		fks = [x for x in r[0][1].split("\n") if "FOREIGN KEY" in x]
		
		constraints_dropped += [(t, x) for x in fks]
		
		# Get names of all Foreign Keys
		fks = [x.split("FOREIGN KEY")[0].replace("CONSTRAINT", "").strip() for x in fks]
		
		fks_all += [(t, x) for x in fks]
	fks_all = list(set(fks_all))
	constraints_dropped = list(set(constraints_dropped))
	
	for fk in fks_all:
		#print fk
		c.execute("ALTER TABLE {0} DROP FOREIGN KEY {1};".format(fk[0], fk[1]))
	
	c.close()
	return constraints_dropped
# FUNCTION END

'''
def riceORFHomologs(orfs):
	if(type(orfs) != type([])):
		orfs = [orfs]
	orfsOrig = orfs
	orfs = [o.replace("LOC_", "") for o in orfs]
	orfs = [o[:o.rfind(".")].upper() if "." in o else o for o in orfs]
	
	orfs = list(set(orfs))
	
	htmls = ["http://www.ricechip.org/cgi-bin/searchrice7.pl?{0}".format(o) for o in orfs]
	
	#orfsMat = np.zeros([len(orfs), len(orfs)])
	orfsMat = defaultdict(lambda: defaultdict(int))
	
	for i in range(len(orfs)):
		s = requests.Session()
		text = s.get(htmls[i]).text
		
		if("Found 0 matches to:" in text):
			print "WARNING: no match to {0}".format(orfsOrig[i])
			orfsMat[orfsOrig[i]] = defaultdict(lambda: -1)
			for o in orfsOrig:
				orfsMat[o][orfsOrig[i]] = -1
			#orfsMat[i,:] = -1
			#orfsMat[:,i] = -1
			
		text = re.sub("\n|\t", " ", text)
		
		pattern = re.compile("http://www.ricechip.org/cgi-bin/.*?\">")
		matches = pattern.findall(text)
		
		matches = [str(m[m.rfind("?")+1:m.rfind(".")]).upper() for m in matches]
		
		for j in range(len(orfs)):
			for k in range(len(orfs)):
				orfsMat[orfsOrig[j]][orfsOrig[k]] += orfs[j] in matches and orfs[k] in matches
				#orfsMat[j][k] += orfs[j] in matches and orfs[k] in matches
	
	return orfsMat, [[orfsMat[o1][o2] for o2 in orfsMat.keys()] for o1 in orfsMat.keys()]
# FUNCTION END
'''

# Rice Homolog Variables - Prevent reading them in multiple times
RiceHomologDF = None
RiceHomologDict = None
RiceORFDict = None
RiceHomologIgnIso = None

# Checks if two [Rice] ORFs are homologs using existing BlastFile and HomologsFile
# - Accepts as parameters...
#   @orf1 - First ORF
#   @orf2 - Second ORF
#   @BlastFile - File containing BLAST results for known ORFS
#   @HomologsFile - File containing currently predicted Homologs for known ORFS (based on cutoff points in BLAST results)
#   @ignoreIsoforms - Flag indicating whether or not isoform extensions should be ignored
#   Returns a tuple containing aboolean indicating whether or not the orf pair is expected to be
#   homologs, and the top BLAST result for the pair
# - Written for general utility purposes 05/31/17 by Shayne Wierbowski
def checkHomologs(orf1, orf2, ORFFile="[REDACTED_PATH]/Plate_Seq/Homologs/RICE_ORFS_ALL.fasta", BlastFile="[REDACTED_PATH]/Plate_Seq/Homologs/blast_RICE_ORFS_ALL_to_RICE_ORFS_ALL", HomologsFile="[REDACTED_PATH]/Plate_Seq/Homologs/Rice_Homologs.txt", ignoreIsoforms=True):
	global RiceHomologDF
	global RiceHomologDict
	global RiceORFDict
	global RiceHomologIgnIso
	
	if(ignoreIsoforms):
		if(not fileExists(ORFFile.replace(".", "_Ignore_Isoforms."))):
			lines = easyReadLines(ORFFile)
			lines = [re.split("-|\.", l)[0] if ">" in l else l for l in lines]
			easyWriteLines(ORFFile.replace(".", "_Ignore_Isoforms."), lines)
		ORFFile = ORFFile.replace(".", "_Ignore_Isoforms.")
		
		if(not fileExists(BlastFile +  "_Ignore_Isoforms")):
			lines = [l.split("\t") for l in easyReadLines(BlastFile)]
			lines = [[re.split("-|\.", l[0])[0]] + [re.split("-|\.", l[1])[0]] + l[2:] for l in lines]
			easyWriteLines(BlastFile + "_Ignore_Isoforms", lines)
		BlastFile = BlastFile + "_Ignore_Isoforms"
		
		if(not fileExists(HomologsFile.replace(".", "_Ignore_Isoforms."))):
			lines = [l.split("\t") for l in easyReadLines(HomologsFile)]
			lines = [[re.split("-|\.", x)[0] for x in l] for l in lines]
			easyWriteLines(HomologsFile.replace(".", "_Ignore_Isoforms."), lines)
		HomologsFile = HomologsFile.replace(".", "_Ignore_Isoforms.")
	
	def simplifyORF(orf):
		return "_".join([str(int(x)) for x in re.split('OS|T|G|-|\.', orf.upper())[1:]])
	# FUNCTION END
	
	# Check BLAST Results
	if(RiceHomologIgnIso != ignoreIsoforms or type(RiceORFDict) == type(None)):
		orfs = easyReadLines(ORFFile)
		orfs = [[orfs[2*i].split(">")[1].strip(), orfs[2*i + 1]] for i in range(len(orfs)/2)]
		orfs = [[l[0], simplifyORF(l[0]), l[1]] for l in orfs]
		RiceORFDict = dict([(x[0], x[1:]) for x in orfs])
		RiceHomologIgnIso = ignoreIsoforms
	
	if(RiceHomologIgnIso != ignoreIsoforms or type(RiceHomologDF) == type(None)):
		RiceHomologDF = pd.read_csv(BlastFile, header=None, sep="\t", names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
		RiceHomologDF["qlength"] = [len(RiceORFDict[x][1]) for x in RiceHomologDF["qseqid"]]
		RiceHomologDF["slength"] = [len(RiceORFDict[x][1]) for x in RiceHomologDF["sseqid"]]
		RiceHomologIgnIso = ignoreIsoforms
	if(ignoreIsoforms):
		df = RiceHomologDF[[sorted([orf1, orf2]) == sorted(x) for x in zip(RiceHomologDF["qseqid"], RiceHomologDF["sseqid"])]].sort(columns=["length", "evalue", "pident"], ascending=[False, True, False])
	else:
		df = RiceHomologDF[[sorted([orf1, orf2])[0] in sorted(x)[0] and sorted([orf1, orf2])[1] in sorted(x)[1]for x in zip(RiceHomologDF["qseqid"], RiceHomologDF["sseqid"])]].sort(columns=["length", "evalue", "pident"], ascending=[False, True, False])
	if(len(df) > 0):
		row = df.iloc[0]
	else:
		row = None
	
	# Check Homolog List
	if(RiceHomologIgnIso != ignoreIsoforms or type(RiceHomologDict) == type(None)):
		RiceHomologDict = dict([(x.split("\t")[0], x.split("\t")[1:]) for x in easyReadLines(HomologsFile)])
		RiceHomologIgnIso = ignoreIsoforms
	if(orf1 == orf2):
		areHomologs = True
	else:
		try:
			if(ignoreIsoforms):
				areHomologs = orf2 in RiceHomologDict[orf1]
			else:
				keys = RiceHomologDict.keys()
				isoforms1 = [x for x in keys if orf1 in x]
				isoforms2 = [x for x in keys if orf2 in x]
				areHomologs = False
				for is1 in isoforms1:
					if(True in [is2 in RiceHomologDict[is1] for is2 in isoforms2]):
						areHomologs = True
						break
		except:
			areHomologs = False
	return (areHomologs, row)
# FUNCTION END

# Convert Plate Index Position to Human Readable Position (e.g. 1 to A01
def IndexPos2HRPos(index):
	row_ids = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
	nrows = 8
	ncols = 12
	index = int(index)
	i3 = (index - 1) / (nrows*ncols)
	i1 = ((index - 1) % (nrows*ncols)) / nrows + 1
	i2 = (index - 1) % nrows
	return "{0}{1:02d}".format(row_ids[nrows*i3 + i2], i1)
# FUNCTION END

# Convert Plate Human Readable Position to Index Position (e.g. A01 to 1
def HRPos2IndexPos(pos):
	row_ids = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
	nrows = 8
	ncols = 12
	i3 = row_ids.index(pos[0]) / 8
	i1 = int(pos[1:])
	i2 = row_ids.index(pos[0]) % 8 + 1
	return (i1-1)*nrows + i2 + i3*nrows*ncols
# FUNCTION END

# Obtains a specified pdb file (using mjm_tools open_pdb)
# - Accepts as parameters...
#   @pdb - The PDB ID for the requested structure
#   @name - The name for the output file (defaults to pdb.pdb
#   @dest - The destination directory for the output file (detault to cwd)
#   @fresh - Flag indicating whether or not existing files should be overwritten
#   Returns the name where the specified file is saved
# - Written for general utility purposes 05/31/17 by Shayne Wierbowski
def get_pdb(pdb, name="", dest=".", fresh=False):
	if(name == ""):
		name = pdb + ".pdb"
	
	if(fresh or not fileExists(dest + "/" + name)):
		try:
			easyWriteLines(dest + "/" + name, [x for x in mjm_tools.open_pdb(pdb)], line_sep="")
		except:
			return None
	return dest + "/" + name
# FUNCTION END

# Updates a local copy of a MYSQL table or selection
# - Accepts as parameters...
#   @db - The MYSQL Database to use
#   @filename - The local file to update
#   @table - The table to select from (only necessary for default selection)
#   @selection - The selection to use (defaults to *)
#   @force - Whether or not the update should be forced
#   @max_age - Maximum age at which the update will be ignored
#   Returns the name where the specified file is saved
# - Written for general utility purposes 09/28/17 by Shayne Wierbowski
def updateLocalMYSQL(db, filename, table="", selection="select * from {0}", force=True, max_age=60*60*24):
	if(force):
		SQLTable2File(db=db, filename=filename, table=table, selection=selection)
	elif(os.path.exists(filename) and time.time() - os.stat(filename) > max_age):
		SQLTable2File(db=db, filename=filename, table=table, selection=selection)
	return filename
# FUNCTION END

rap_msu_conversion = None
rice_all_ids = None

# Converts Rice RAP-DB ID to MSU ID
# - Accepts as parameters...
#   @id - The ID to convert
#   @ignoreIsoforms - A boolean flag indicating whether or not isoforms (.x) extensions should be stripped
#   Returns the converted ID
# - Written for general utility purposes 09/28/17 by Shayne Wierbowski
def rap2msu(id, ignoreIsoforms=True):
	global rap_msu_conversion
	if(rap_msu_conversion == None):
		f = updateLocalMYSQL(db="rice_interactome", filename="[REDACTED_PATH]/Plate_Seq/MSU_v7_ALL/rap2msu.txt", table="rap_msu_conversion")
		df = pd.read_csv(f, sep="\t")
		rap_msu_conversion = [dict(zip(df["ORF_Name"], df["LOCUS_NAME"])),
									 dict(zip(df["ORF_Name"], [x.split(".")[0] for x in df["LOCUS_NAME"]])),
									 dict(zip(df["LOCUS_NAME"], df["ORF_Name"])),
									 dict(zip([x.split(".")[0] for x in df["LOCUS_NAME"]], df["ORF_Name"]))]
	if(ignoreIsoforms):
		return rap_msu_conversion[1][id.split(".")[0]]
	else:
		return rap_msu_conversion[0][id]
# FUNCTION END

# Converts Rice MSU ID to RAP-DB
# - Accepts as parameters...
#   @id - The ID to convert
#   @ignoreIsoforms - A boolean flag indicating whether or not isoforms (.x) extensions should be stripped
#   Returns the converted ID
# - Written for general utility purposes 09/28/17 by Shayne Wierbowski
def msu2rap(id, ignoreIsoforms=True):
	global rap_msu_conversion
	if(rap_msu_conversion == None):
		f = updateLocalMYSQL(db="rice_interactome", filename="[REDACTED_PATH]/Plate_Seq/MSU_v7_ALL/rap2msu.txt", table="rap_msu_conversion")
		df = pd.read_csv(f, sep="\t")
		rap_msu_conversion = [dict(zip(df["ORF_Name"], df["LOCUS_NAME"])),
									 dict(zip(df["ORF_Name"], [x.split(".")[0] for x in df["LOCUS_NAME"]])),
									 dict(zip(df["LOCUS_NAME"], df["ORF_Name"])),
									 dict(zip([x.split(".")[0] for x in df["LOCUS_NAME"]], df["ORF_Name"]))]
	if(ignoreIsoforms):
		return rap_msu_conversion[3][id.split(".")[0]]
	else:
		return rap_msu_conversion[2][id]
# FUNCTION END

# Returns the merged 5.1 and 8.1 horfeome (ORF ID + Seq)
# - Accepts as parameters...
#   @as_dict - Boolean indicating whether to return a dictionary or a df
# - Written for general utility purposes 01/10/18 by Shayne Wierbowski
def horfeomeDict(as_df=False):
	horf81 = pd.read_csv(updateLocalMYSQL(db="human_orfeome", filename="[REDACTED_PATH]/bin/Local_MySQL_Copies/81_horfeome.txt", table="81_horfeome"), sep="\t")
	horf51 = pd.read_csv(updateLocalMYSQL(db="human_orfeome", filename="[REDACTED_PATH]/bin/Local_MySQL_Copies/51_horfeome_cDNA.txt", table="51_horfeome_cDNA"), sep="\t")
	
	keep = horf81[horf81["ORF"].isin(horf81["ORF"])]
	keep2 = horf51[horf51["ORF"].isin(set(horf51["ORF"]).difference(keep["ORF"].values))]
	
	horf_dict81 = keep.set_index('ORF').to_dict()['Sequence']
	horf_dict51 = keep2.set_index('ORF').to_dict()['CDS_SEQ']
	
	horf_dict = horf_dict81.copy()
	horf_dict.update(horf_dict51)
	
	if(not as_df):
		return horf_dict
	else:
		return pd.DataFrame([[x, horf_dict[x]] for x in horf_dict], columns=["ORF", "Seq"]).sort_values("ORF")
# FUNCTION END

# Converts (tries to convert) Rice MSU ID to UniProt information by manually searching UniProt
# - Accepts as parameters...
#   @id - The ID to convert
#   @use_rapdbBackup - A boolean flag indicating whether or not the alternative ID should be used as a backup
#   @ignoreIsoforms - A boolean flag indicating whether or not isoforms (.x) extensions should be stripped
#   Returns the converted ID
# - Written for general utility purposes 10/02/17 by Shayne Wierbowski
def MSU2UniProt(id, use_rapdbBackup=True, ignoreIsoforms=True):
	if(ignoreIsoforms):
		id = id.split(".")[0]
	
	html = "http://www.uniprot.org/uniprot/?query={0}+AND+organism:Oryza%20sativa&sort=score".format(id)
	
	# Start Session to Query UniProt
	s = requests.Session()
	text = s.get(html).text
	
	text = re.sub("\n|\t", " ", text)
	
	# Pull out entire Results Table
	pattern = re.compile("</script></td></tr></thead><tbody><tr .*?</tbody?")
	matches = pattern.findall(text)
	
	if(len(matches) == 0):
		return
	
	text = matches[0]
	
	# Pull out individual results
	pattern = re.compile("<tr .*?</tr>")
	matches = pattern.findall(text)
	
	if(len(matches) == 0):
		return
	
	#if(len(matches) > 1):
	#	print "FOUND GREATER THAN 1 MATCH FOR {0} MANUALLY CHECK THIS".format(id)
	
	# Parse individual results
	for i in range(len(matches)):
		# Remove HTML Blocks
		matches[i] = re.sub("</td>", "\n", matches[i])
		matches[i] = re.sub("<script>.*?</script>", "\t", matches[i])
		matches[i] = re.sub("<.*?>", "\t", matches[i])
		matches[i] = [[str(z) for z in x.strip().split("\t") if z != ""] for x in matches[i].strip().split("\n") if x != ""]
		
		# Make sure entry 3 and 4 split properly
		matches[i][3] = flatten([x.split(",") for x in matches[i][3]])
		matches[i][3] = [x.strip().replace("(", "").replace(")", "") for x in matches[i][3]]
		matches[i][4] = flatten([x.split(",") for x in matches[i][4]])
		matches[i][4] = [x.strip().replace("(", "").replace(")", "") for x in matches[i][4]]
		
		# Go to Uniprot Page
		html = "http://www.uniprot.org/uniprot/{0}".format(matches[i][0][0])
		
		text = s.get(html).text
		
		# Pull out all non-string hyperlinks
		pattern = re.compile("^(<a href=\"http://string-db.org.*?</a>)")
		#print pattern.findall(text)
		text2 = re.sub("<a href=\"http://string-db.org.*?</a>", "", text)
		
		# Skip this match (does not contain id outside of string interactions)
		if(not id.upper() in text2.upper()):
			matches[i].append(False)
		else:
			matches[i].append(True)
	
	matches = [x for x in matches if x[-1] == True]
	
	if(len(matches) == 0 and use_rapdbBackup):
		return RAPDB2UniProt(msu2rap(id, ignoreIsoforms=ignoreIsoforms), use_msuBackup=False)
	
	return matches
# FUNCTION END

# Converts (tries to convert) Rice MSU ID to UniProt information by manually searching UniProt
# - Accepts as parameters...
#   @id - The ID to convert
#   @use_rapdbBackup - A boolean flag indicating whether or not the alternative ID should be used as a backup
#   @ignoreIsoforms - A boolean flag indicating whether or not isoforms (.x) extensions should be stripped
#   Returns the converted ID
# - Written for general utility purposes 10/02/17 by Shayne Wierbowski
def RAPDB2UniProt(id, use_msuBackup=True, ignoreIsoforms=True):
	if(ignoreIsoforms):
		id = id.split(".")[0]
	
	html = "http://www.uniprot.org/uniprot/?query={0}+AND+organism:Oryza%20sativa&sort=score".format(id)
	
	# Start Session to Query UniProt
	s = requests.Session()
	text = s.get(html).text
	
	text = re.sub("\n|\t", " ", text)
	
	# Pull out entire Results Table
	pattern = re.compile("</script></td></tr></thead><tbody><tr .*?</tbody?")
	matches = pattern.findall(text)
	
	if(len(matches) == 0):
		return
	
	text = matches[0]
	
	# Pull out individual results
	pattern = re.compile("<tr .*?</tr>")
	matches = pattern.findall(text)
	
	if(len(matches) == 0):
		return
	
	#if(len(matches) > 1):
	#	print "FOUND GREATER THAN 1 MATCH FOR {0} MANUALLY CHECK THIS".format(id)
	
	# Parse individual results
	for i in range(len(matches)):
		# Remove HTML Blocks
		matches[i] = re.sub("</td>", "\n", matches[i])
		matches[i] = re.sub("<script>.*?</script>", "\t", matches[i])
		matches[i] = re.sub("<.*?>", "\t", matches[i])
		matches[i] = [[str(z) for z in x.strip().split("\t") if z != ""] for x in matches[i].strip().split("\n") if x != ""]
		
		# Make sure entry 3 and 4 split properly
		matches[i][3] = flatten([x.split(",") for x in matches[i][3]])
		matches[i][3] = [x.strip().replace("(", "").replace(")", "") for x in matches[i][3]]
		matches[i][4] = flatten([x.split(",") for x in matches[i][4]])
		matches[i][4] = [x.strip().replace("(", "").replace(")", "") for x in matches[i][4]]
		
		# Go to Uniprot Page
		html = "http://www.uniprot.org/uniprot/{0}".format(matches[i][0][0])
		
		text = s.get(html).text
		
		# Pull out all non-string hyperlinks
		pattern = re.compile("^(<a href=\"http://string-db.org.*?</a>)")
		#print pattern.findall(text)
		text2 = re.sub("<a href=\"http://string-db.org.*?</a>", "", text)
		
		# Skip this match (does not contain id outside of string interactions)
		if(not id.upper() in text2.upper()):
			matches[i].append(False)
		else:
			matches[i].append(True)
	
	matches = [x for x in matches if x[-1] == True]
	
	if(len(matches) == 0 and use_msuBackup):
		return MSU2UniProt(rap2msu(id), use_rapdbBackup=False)
	
	return matches
# FUNCTION END

# Runs analysis on Sanger Sequencing data
# - Accepts as parameters...
#   @startingDir - The base directory where the Sequencing Results are stored
#   @expected - Path to a file containing the expected results. File should include...
#               Source Plate, Source Well, Expected ORF, SEQ Plate *, SEQ Plate Position, + any extraneous
#   @blastdb - The path to the fasta to use in the BLAST
#   @out - Optional path to write output df
#   @seqFilePattern - Pattern for matching sequencing result files (defaults to like X.txt)
#   @platePattern - Pattern for matching sequencing plate number from file name(defaults to like plate[X])
#   @posPattern - Pattern for matching  sequencing position from file name (defaults to like X[A01])
#   Returns the DataFrame containing analyisis of the Sanger results
# - Written for general utility purposes 10/02/17 by Shayne Wierbowski
def SangerAnalyze(startingDir, expected, blastdb, out=None, seqFilePattern="*.txt", platePattern="Plate[0-9]", posPattern="[A-H][0-9][0-9]", fresh=True):
	# Find Sequencing Files
	files = []
	for root, dirnames, filenames in os.walk(startingDir):
		files += glob.glob("{0}/*.txt".format(root))
	
	# Load Sequencing Files DF
	df = pd.DataFrame(files, columns=["File"])
	df["SEQ Plate Position"] = df["File"].map(lambda x: re.findall(posPattern, x)[-1])
	df["SEQ Plate #"] = df["File"].map(lambda x: int(re.findall(platePattern, x)[-1].replace("Plate", "")))
	df["Sequence"] = df["File"].map(lambda x: "".join(call("cat {0} | grep -v \">\"".format(x))))
	
	# Join with Expectation DF
	df = pd.read_csv(expected, sep="\t").join(df.set_index(["SEQ Plate #", "SEQ Plate Position"]), on=["SEQ Plate #", "SEQ Plate Position"], how="left")
	
	# Merge into single tmp FASTA
	tmpFasta = reserveTemp()
	fastaEntries = [">{0}_{1}\n{2}".format(x[0], x[1], x[2]) for x in zip(df["Source Plate"], df["Source Well"], df["Sequence"])]
	easyWriteLines(tmpFasta, [x for x in fastaEntries if not "nan" in x])
	#call("cat {0} > {1}".format(" ".join(matches), tmpFasta))
	
	# RUN BLAST
	blastResults = runBlast(tmpFasta, db=blastdb, fresh=fresh)[0]
	
	# Read BLAST Results
	blastResults = pd.read_csv(blastResults, sep="\t", header=None, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
	blastResults["Source Plate"] = [int(x.split("_")[0]) for x in blastResults["qseqid"]]
	blastResults["Source Well"] = [x.split("_")[1] for x in blastResults["qseqid"]]
	#blastResults["Dest Plate"] = [int(re.findall(platePatten, x)[-1].replace("Plate", "")) for x in blastResults["qseqid"]]
	#blastResults["Dest Well"] = [re.findall(posPattern, x)[-1] for x in blastResults["qseqid"]]
	
	# Join BLAST Results
	df = df.join(blastResults.set_index(["Source Plate", "Source Well"]), on=["Source Plate", "Source Well"])
	df["Matches Expected"] = df["Expected ORF"] == df["sseqid"]
	
	# Sort and drop_duplicates
	df = df.sort_values(by=["Source Plate", "Source Well", "evalue"], ascending=True)
	
	if(not out == None):
		df.to_csv(out, sep="\t", index=None)
	
	return df
# FUNCTION END


# Generates the reverse complement of a DNA sequence
# - Accepts as parameters...
#   @seq - The sequence to reverse complement
#   Returns the reverse complement
# - Written for general utility purposes 1/03/18 by Shayne Wierbowski
def reverse_complement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}
	try:
		return "".join(complement.get(base, base) for base in reversed(seq))
	except IndexError:
		print "ERROR: BAD DNA SEQ"
		return None
# FUNCTION END

# Reads in a fasta format file as a ID --> Seq Dictionary
# - Accepts as parameters...
#   @filename - The fasta file to read
#   Returns the dictionary representation of the file
# - Written for general utility purposes 1/17/18 by Shayne Wierbowski
def fasta2dict(filename, only_len=False, key_transform=str, keep=None):
	#lines = easyReadLines(filename)
	if(keep != None):
		keep = set(keep)
	
	with open(filename, "r") as f:
		res = dict()
		#IDs = []
		#Seqs = []
		cur = None
		if(only_len):
			tmp = 0
		else:
			tmp = ""
		for l in f:
			if(">" in l):
				#IDs.append(l.split(">")[1].strip())
				#if(tmp != ""):
				#	Seqs.append(tmp)
				if(cur == None):
					cur = key_transform(l.split(">")[1].strip())
				else:
					if(keep == None or cur in keep):
						res[cur] = tmp
					cur = key_transform(l.split(">")[1].strip())
				if(only_len):
					tmp = 0
				else:
					tmp = ""
			else:
				if(only_len):
					tmp += len(l.strip())
				else:
					tmp += l.strip()
		#Seqs.append(tmp)
		if(keep == None or cur in keep):
			res[cur] = tmp
	#return dict(zip(IDs, Seqs))
	return res
# FUNCTION END

# Reads in a fastq format file as a ID --> Seq Dictionary
# - Accepts as parameters...
#   @filename - The fasta file to read
#   Returns the dictionary representation of the file
# - Written for general utility purposes 1/17/18 by Shayne Wierbowski
def fastq2dict(filename):
	l = easyReadLines(filename)
	r = dict(zip([x[1:].strip() for x in l[::4]], [dict(zip(["SEQ", "INFO", "QUAL"], x)) for x in zip(l[1::4], l[2::4], l[3::4])]))
	del l
	return r
# FUNCTION END

def dict2fasta(filename, dict):
	easyWriteLines(filename, [">" + str(x) + "\n" + str(dict[x]) for x in dict])
# FUNCTION END

def dict2fastq(filename, d, qual="I"):
	if(type(d.values()[0]) != dict):
		easyWriteLines(filename, ["@" + str(x) + "\n" + str(d[x]) + "\n+\n" + qual*len(str(d[x])) for x in d])
	else:
		lines = []
		for k, v in d.iteritems():
			lines.append("@" + str(k))
			lines.append(v["SEQ"])
			lines.append(v["INFO"])
			lines.append(v["QUAL"])
		easyWriteLines(filename, lines)
# FUNCTION END

# Download GRCh38 from ensemble ftp server (fasta + gtf annotations) and generates relevant local index (just change the release-XX and GRChXX for other versions)
# - Accepts as parameters...
#   @out_dir - The output directory to store the downloaded files in
#   @fresh - Boolean flag indicating whether or not pre-existing copies should be replaced
#   @index - Boolean flag indicating whether or not to separate out the files into an index
#   Returns the path to the final merged fasta
# - Written for general utility purposes 01/24/18 by Shayne Wierbowski
def download_HG38(out_dir, fresh=False, index=True):
	call("mkdir -p {0}".format(out_dir))
	starting_dir = os.getcwd()
	os.chdir("{0}".format(out_dir))
	if(fresh or not os.path.exists("GRCh38_r85.all.fa")):
		call("wget ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz")
		call("gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > GRCh38_r85.all.fa")
		call("rm H*.gz")
	if(fresh or not os.path.exists("Homo_sapiens.GRCh38.85.gtf")):
		call("wget ftp://ftp.ensembl.org/pub/release-85/gtf/homo_sapiens/Homo_sapiens.GRCh38.85.gtf.gz")
		call("gunzip Homo_sapiens.GRCh38.85.gtf.gz")
	os.chdir(starting_dir)
	if(index):
		index_HG(out_dir, "GRCh38", "GRCh38_r85.all.fa")
	return out_dir + "/" + "GRCh38_r85.all.fa"
# FUNCTION END

# Download GRCh37 (hg19) from ensemble ftp server (fasta + gtf annotations) and generates relevant local index (just change the release-XX and GRChXX for other versions)
# - Accepts as parameters...
#   @out_dir - The output directory to store the downloaded files in
#   @fresh - Boolean flag indicating whether or not pre-existing copies should be replaced
#   @index - Boolean flag indicating whether or not to separate out the files into an index
#   Returns the path to the final merged fasta
# - Written for general utility purposes 01/24/18 by Shayne Wierbowski
def download_HG37(out_dir, fresh=False, index=True):
	call("mkdir -p {0}".format(out_dir))
	starting_dir = os.getcwd()
	os.chdir("{0}".format(out_dir))
	if(fresh or not os.path.exists("GRCh37_r75.all.fa")):
		call("wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.{1..22}.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.X.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.Y.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.MT.fa.gz")
		call("gunzip -c Homo_sapiens.GRCh37.75.dna.chromosome.* > GRCh37_r75.all.fa")
		call("rm H*.gz")
	if(fresh or not os.path.exists("Homo_sapiens.GRCh37.75.gtf")):
		call("wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz")
		call("gunzip Homo_sapiens.GRCh37.75.gtf.gz")
	os.chdir(starting_dir)
	if(index):
		index_HG(out_dir, "GRCh37", "GRCh37_r75.all.fa")
	return out_dir + "/" + "GRCh37_r75.all.fa"
# FUNCTION END

# Modified the above to download HG37, Specifically release 64. For Junke / exact reference
# genome build used by TumorFusions / PRADA
# - Written for general utility purposes 05/20/19 by Shayne Wierbowski
def download_HG37_64(out_dir, fresh=False, index=True):
	call("mkdir -p {0}".format(out_dir))
	starting_dir = os.getcwd()
	os.chdir("{0}".format(out_dir))
	if(fresh or not os.path.exists("GRCh37_r64.all.fa")):
		call("wget ftp://ftp.ensembl.org/pub/release-64/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.64.dna.chromosome.{1..22}.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-64/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.64.dna.chromosome.X.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-64/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.64.dna.chromosome.Y.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-64/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.64.dna.chromosome.MT.fa.gz")
		call("gunzip -c Homo_sapiens.GRCh37.64.dna.chromosome.* > GRCh37_r64.all.fa")
		call("rm H*.gz")
	if(fresh or not os.path.exists("Homo_sapiens.GRCh37.64.gtf")):
		call("wget ftp://ftp.ensembl.org/pub/release-64/gtf/homo_sapiens/Homo_sapiens.GRCh37.64.gtf.gz")
		call("gunzip Homo_sapiens.GRCh37.64.gtf.gz")
	os.chdir(starting_dir)
	if(index):
		index_HG(out_dir, "GRCh37", "GRCh37_r64.all.fa")
	return out_dir + "/" + "GRCh37_r64.all.fa"
# FUNCTION END

# Download GRCh36 (hg18) from ensemble ftp server (fasta + gtf annotations) and generates relevant local index (just change the release-XX and GRChXX for other versions)
# - Accepts as parameters...
#   @out_dir - The output directory to store the downloaded files in
#   @fresh - Boolean flag indicating whether or not pre-existing copies should be replaced
#   @index - Boolean flag indicating whether or not to separate out the files into an index
#   Returns the path to the final merged fasta
# - Written for general utility purposes 01/24/18 by Shayne Wierbowski
def download_HG36(out_dir, fresh=False, index=True):
	call("mkdir -p {0}".format(out_dir))
	starting_dir = os.getcwd()
	os.chdir("{0}".format(out_dir))
	if(fresh or not os.path.exists("GRCh36_r54.all.fa")):
		call("wget ftp://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.54.dna.chromosome.{1..22}.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.54.dna.chromosome.X.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.54.dna.chromosome.Y.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.54.dna.chromosome.MT.fa.gz")
		call("gunzip -c Homo_sapiens.NCBI36.54.dna.chromosome.* > GRCh36_r54.all.fa")
		call("rm H*.gz")
	if(fresh or not os.path.exists("Homo_sapiens.NCBI36.54.gtf")):
		call("wget ftp://ftp.ensembl.org/pub/release-54/gtf/homo_sapiens/Homo_sapiens.NCBI36.54.gtf.gz")
		call("gunzip Homo_sapiens.NCBI36.54.gtf.gz")
	os.chdir(starting_dir)
	if(index):
		index_HG(out_dir, "GRCh36", "GRCh36_r54.all.fa")
	return out_dir + "/" + "GRCh36_r54.all.fa"
# FUNCTION END

# Download GRCm38 from ensemble ftp server (fasta + gtf annotations) and generates relevant local index (just change the release-XX and GRChXX for other versions)
# - Accepts as parameters...
#   @out_dir - The output directory to store the downloaded files in
#   @fresh - Boolean flag indicating whether or not pre-existing copies should be replaced
#   @index - Boolean flag indicating whether or not to separate out the files into an index
#   Returns the path to the final merged fasta
# - Written for general utility purposes 04/11/19 by Shayne Wierbowski
def download_MG38(out_dir, fresh=False, index=True):
	call("mkdir -p {0}".format(out_dir))
	starting_dir = os.getcwd()
	os.chdir("{0}".format(out_dir))
	if(fresh or not os.path.exists("GRCm38_r96.all.fa")):
		call("wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.{1..19}.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.X.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.Y.fa.gz")
		call("wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.MT.fa.gz")
		call("gunzip -c Mus_musculus.GRCm38.dna.chromosome.* > GRCm38_r96.all.fa")
		call("rm M*.gz")
	if(fresh or not os.path.exists("Mus_musculus.GRCm38.96.gtf")):
		call("wget ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.96.gtf.gz")
		call("gunzip Mus_musculus.GRCm38.96.gtf.gz")
	os.chdir(starting_dir)
	if(index):
		index_HG(out_dir, "GRCm38", "GRCm38_r96.all.fa")
	return out_dir + "/" + "GRCm38_r96.all.fa"
# FUNCTION END

# Creates a local index of a particular human genome build (separates out single folder per chromosome and breaks each chromosome into binned ranges)
# - Accepts as parameters...
#   @out_dir - Output directory to create the index in (i.e. where the genome build is stored)
#   @name - Name for the index (i.e. sub-directory the index will be organized in)
#   @fasta - Filename for the reference fasta containing the genome build
#   @bin_size - The bin size to break individual chromosome sequences into
#   Returns the final path to the index 
# - Written for general utility purposes 01/24/18 by Shayne Wierbowski
def index_HG(out_dir, name, fasta, bin_size=1000000):
	call("mkdir -p {0}".format(out_dir))
	starting_dir = os.getcwd()
	os.chdir("{0}".format(out_dir))
	name = "{0}_{1}_Index".format(name, bin_size)
	call("mkdir {0}".format(name))
	fasta = fasta2dict(fasta)
	for x in fasta:
		x2 = "Chr{0}".format(x.split()[0])
		call("mkdir {0}/{1}".format(name, x2))
		for bin_start in range(0, len(fasta[x]), bin_size):
			easyWriteLines("{0}/{1}/{2}".format(name, x2, bin_start), fasta[x][bin_start:bin_start + bin_size])
	os.chdir(starting_dir)
	return out_dir + "/" + name
# FUNCTION END

HG_fasta_dict = defaultdict(lambda: defaultdict(dict))
HG_files_dict = defaultdict(dict)

# Fetches a desired sub-sequence from a locally indexed human genome build
# - Accepts as parameters...
#   @chrom - The chromosome of interest (for current HG buiilds, should be "ChrX", but will automatically add the "Chr" part)
#   @start - The 1-indexed sequence start site (inclusive)
#   @end - The 1-indexed sequence end site (inclusive)
#   @version - The name of the locally indexed genome build (i.e. folder name in [REDACTED_PATH]/resources)
#   @index_source - Direct path to index directory (i.e. [REDACTED_PATH]/GRCh38/GRCh38_1000000_Index) (Optional)
#   @reverse - Boolean flag for returning information from the negative strand (i.e. reverse complement)
#   Returns the DNA sequence for the desired range
# - Written for general utility purposes 01/24/18 by Shayne Wierbowski
def fetch_HG_Seq(chrom, start, end, version="GRCh38", index_source=None, reverse=False):
	#start_t = time.time()
	# Update input parameters
	
	# Check for mis-ordered start / end (or negative strand queries)
	if(start > end):
		if(not reverse):
			print "ERROR: start cannot be greater than end"
			return ""
		else:
			tmp = start
			start = end
			end = tmp
	
	# Move 1-index coordinates to 0-index
	start = int(start - 1)
	end = int(end - 1)
	
	# Allow integer based chromosomes
	if(not "chr" in str(chrom).lower()):
		chrom = "Chr" + str(chrom)
	elif(not "Chr" in str(chrom)):
		#chrom = str(chrom).lower().replace("chr", "Chr")
		chrom = str(chrom).replace("chr", "Chr")
	
	# Allow integer based genome versions (only for 37 and 38, and 36)
	# NOTE: 2020_02_12
	# Why did I need to make this exception in the first place? Presumably because of the
	# mice genomes? To add to this, also need to modify read_HG_GTF.
	# This is a generally non-scalable solution since it could be misused if there were
	# ever a version X from a different organism included in my resources
	if(str(version) == "37" or str(version) == "38" or str(version) == "36"):
		version = "GRCh" + str(version)
	
	# Try to fetch local genome version indexes
	if(index_source == None):
		index_source = glob.glob("[REDACTED_PATH]/resources/{0}/*Index".format(version))[0]
	bin_size = int(index_source.split("_")[-2])
	
	# Move to index location
	starting_dir = os.getcwd()
	os.chdir(index_source)
	
	# Identify relevant index files
	#files = sorted([x for x in glob.glob("{0}/*".format(chrom)) if start - bin_size <= int(x.split("/")[-1]) <= end], key=lambda x: int(x.split("/")[-1]))
	global HG_files_dict
	try:
		files = HG_files_dict[version][chrom]
	except:
		files = [(f, int(f.split("/")[-1])) for f in glob.glob("{0}/*".format(chrom))]
		HG_files_dict[version][chrom] = files
	files = sorted([x[0] for x in files if start - bin_size <= x[1] <= end], key=lambda x: x[1])
	
	# Keeps everything that has been loaded previously in memory (saves on time re-accessing the same positions over and over est. >1000 fold)
	global HG_fasta_dict
	def do(x):
		try:
			return HG_fasta_dict[version][chrom][x]
		except KeyError:
			#HG_fasta_dict[version][chrom][x] = call("cat {0}".format(x))[0]
			HG_fasta_dict[version][chrom][x] = open(x).read().strip()
			return HG_fasta_dict[version][chrom][x]
	# FUNCTION END
	#seq = call("cat {0}".format(" ".join(files)))[0]
	
	#print chrom, start, end, version, files
	seq = "".join(do(x) for x in files)
	try:
		seq_start = int(files[0].split("/")[-1])
	except IndexError:
		print chrom, start, end, version, index_source, reverse
		print files
		raise
	
	os.chdir(starting_dir)
	
	#end_t = time.time()
	#print "fetch_HG", end_t - start_t
	# Return desired range from index files
	if(not reverse):
		return seq[start - seq_start:end - seq_start + 1]
	else:
		# Option for returning negative strand
		return reverse_complement(seq[start - seq_start:end - seq_start + 1])
# FUNCTION END

HG_GTF_dict = defaultdict(dict)

# Reads in a human genome GTF file and as a pandas df
# - Accepts as parameters...
#   @version - The name of the genome build (i.e. folder name in [REDACTED_PATH]/resources)
#   @chrom - The chromosome of interest (for current HG buiilds, should be "ChrX", but will automatically add the "Chr" part)
#   @source - Direct path to a specific GTF file (Optional)
#   Returns a pandas df representing the specified GTF file (or specific chromosome from it)
# - Written for general utility purposes 01/24/18 by Shayne Wierbowski
def read_HG_GTF(version="GRCh38", chrom=None, source=None):
	# For keeping data in memory / cutting down on re-access time
	global HG_GTF_dict
	
	# Parse Chromosome Input
	chrom = str(chrom)
	if("chr" in chrom.lower()):
		chrom = chrom[3:]
	
	# Allow integer based genome versions (only for 37 and 38, and 36)
	# NOTE: 2020_02_12
	# Why did I need to make this exception in the first place? Presumably because of the
	# mice genomes? To add to this, also need to modify read_HG_GTF.
	# This is a generally non-scalable solution since it could be misused if there were
	# ever a version X from a different organism included in my resources
	if(str(version) == "37" or str(version) == "38" or str(version) == "36"):
		version = "GRCh" + str(version)
	
	# Try to infer GTF file based on version / [REDACTED_PATH]/resources content
	if(source == None):
		source = glob.glob("[REDACTED_PATH]/resources/{0}/*gtf".format(version))[0]
	
	# Check if this part has already been loaded
	try:
		df = HG_GTF_dict[version][chrom]
		return df
	except KeyError:
		pass
	
	if(chrom == "None"):
		df = pd.DataFrame([x.split("\t") for x in call("grep -v \"#\" {0}".format(source))], columns=["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"])
		#df = pd.DataFrame([x.strip().split("\t") for x in open(source, "r") if x[0] != "#"], columns=["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"])
	else:
		#df = pd.DataFrame([x.split("\t") for x in call("grep -v \"#\" {0} | grep \"^{1}\"".format(source, chrom))], columns=["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"])
		df = pd.DataFrame([x.split("\t") for x in call("grep -v \"#\" {0} | grep \"^{1}\s\"".format(source, chrom))], columns=["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"])
	
	df["Attributes"] = df["Attributes"].map(lambda x: dict([(str(y.split()[0]), str(y.split()[1].strip("\""))) for y in x.split(";") if not y == ""]))
	for att in set(flatten(df["Attributes"].map(lambda x: list(x)))):
		def do(x):
			if(att in x):
				return x[att]
			else:
				return None
		# FUNCTION END
		df[att] = df["Attributes"].map(lambda x: do(x))
	
	df["Start"] = df["Start"].map(lambda x: int(x))
	df["End"] = df["End"].map(lambda x: int(x))
	
	HG_GTF_dict[version][chrom] = df[[x for x in list(df) if x != "Attributes"]]
	
	return HG_GTF_dict[version][chrom]
# FUNCTION END

# Converts a number of bits to a human readable size
# Code copied directly from https://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
def hr_size(num, suffix='B'):
	for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
		if abs(num) < 1024.0:
			return "%3.1f %s%s" % (num, unit, suffix)
		num /= 1024.0
	return "%.1f %s%s" % (num, 'Yi', suffix)
# FUNCTION END

# Returns the size of an object in python
# Code copied directly from https://goshippo.com/blog/measure-real-size-any-python-object/
# Added modification for human readable flag
def get_size(obj, seen=None, h=True):
	"""Recursively finds size of objects"""
	size = sys.getsizeof(obj)
	if seen is None:
		seen = set()
	obj_id = id(obj)
	if obj_id in seen:
		 return 0
	# Important mark as seen *before* entering recursion to gracefully handle
	# self-referential objects
	seen.add(obj_id)
	if isinstance(obj, dict):
		size += sum([get_size(v, seen, h=False) for v in obj.values()])
		size += sum([get_size(k, seen, h=False) for k in obj.keys()])
	elif hasattr(obj, '__dict__'):
		size += get_size(obj.__dict__, seen, h=False)
	elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
		size += sum([get_size(i, seen, h=False) for i in obj])
	if(h):
		return hr_size(size)
	else:
		return size
# FUNCTION END

# Class for handling my preferred typ of multiprocessing
# Basically a clone of Pool / Pool.map() but written by me
# Heavilly based on code available here...
# https://stackoverflow.com/questions/10415028/how-can-i-recover-the-return-value-of-a-function-passed-to-multiprocessing-proce
class Multiprocessor():
	
	def __init__(self, max_cores=10, niceness=0, progress_bar=list, random_seeds=True):#, raise_errors=False):
		self.__max_cores = max_cores
		self.__niceness = niceness
		self.__progress_bar = progress_bar
		self.__random_seeds = random_seeds
		#self.__raise_errors=raise_errors
		self.__processes = []
		self.__outqueue = Queue()
		self.__errqueue = Queue()
		self.__lock = Lock()
		self.__counter = Value("i", 0)
		self.__running = Value("i", 0)
		self.__rets = []
		self.__errs = []
		self.__int_rets = []
		self.__int_errs = []
		self.__last_wait = time.time()
	# FUNCTION END
	
	def get_outs(self):
		return self.__rets
	# FUNCTION END
	
	def get_errs(self):
		return self.__errs
	# FUNCTION END

	def get_i_outs(self):
		return self.__int_rets
	# FUNCTION END
	
	def get_i_errs(self):
		return self.__int_errs
	# FUNCTION END
	
	@staticmethod
	def _wrapper(func, outqueue, errqueue, lock, counter, running, seed, niceness, args, kwargs):
		# Handle jobs that end in error
		err = (counter, None)
		try:
			if(seed != None):
				np.random.seed(seed)
			if(niceness != 0):
				cur_nice = os.nice(0)
				os.nice(max([0, niceness - cur_nice]))
			ret = (counter, func(*args, **kwargs))
		except KeyboardInterrupt:
			outqueue.put((counter, None))#, block=True, timeout=None)
			outqueue.close()
			outqueue.join_thread()
			errqueue.put(counter, traceback.format_exc())#, block=True, timeout=None)
			errqueue.close()
			errqueue.join_thread()
			raise
		except:
			#if(raise_errors):
			#	raise
			ret = (counter, None)
			err = (counter, traceback.format_exc())
			print err[1]
		
		# Add output / err to queues
		#for i in range(1000):
		outqueue.put(ret)#, block=True, timeout=None)
		outqueue.close()
		outqueue.join_thread()
		errqueue.put(err)#, block=True, timeout=None)
		errqueue.close()
		errqueue.join_thread()
		
		# Decrement running
		lock.acquire()
		running.value -= 1
		lock.release()
		return
	# FUNCTION END
	
	def run(self, func, *args, **kwargs):
		# Wait until runnning < max_cores
		self.__lock.acquire()
		while(self.__running.value > self.__max_cores):
			self.__lock.release()
			if(len(self.__processes) == 0):
				running.value = 1
				self.__lock.acquire()
				print "IT DID NOT WORK"
				continue
			time.sleep(5)
			self.intermediate_wait()
			self.__lock.acquire()
		
		# Increment Counter / Running Valuess for process about to be spawned
		self.__running.value += 1
		self.__counter.value += 1
		#print "R:", self.__running.value
		#print "C:", self.__counter.value
		#if(self.__counter.value % self.__max_cores == 0):
		#	pass
		#	self.intermediate_wait()
		if(not self.__outqueue.empty() or not self.__errqueue.empty()):
			if(time.time() - self.__last_wait > 10):
				self.intermediate_wait()
		self.__lock.release()
		
		# Correct Random Seeds
		if(self.__random_seeds):
			seed = np.random.randint(0, 10000000)
		else:
			seed = None
		
		# Submit job
		#args2 = [func, self.__outqueue, self.__errqueue,  self.__lock, int(self.__counter.value), self.__running, self.__raise_errors, args, kwargs]
		args2 = [func, self.__outqueue, self.__errqueue,  self.__lock, int(self.__counter.value), self.__running, seed, self.__niceness, args, kwargs]
		#self.__processes.append(Process(target=self._wrapper, args=args2))
		#self.__processes[len(self.__processes) - 1].start()
		
		p = Process(target=self._wrapper, args=args2)
		p.daemon = True
		self.__processes.append(p)
		p.start()
	# FUNCTION END
	
	def wait(self):
		# Wait as usual
		try:
			return self.wait_inner()
		# Safely terminate children is something is interrupted
		except (KeyboardInterrupt, SystemExit):
			for p in self.__processes:
				p.terminate()
			raise
	# FUNCTION END
	
	def intermediate_wait(self):
		joined = True
		while(joined == True):
			if(not self.__outqueue.empty() or not self.__errqueue.empty()):
				ret = self.__outqueue.get(block=True, timeout=None)
				time.sleep(0.1)
				err = self.__errqueue.get(block=True, timeout=None)
				time.sleep(0.1)
				self.__int_rets.append(ret)
				self.__int_errs.append(err)
				continue
			else:
				break
			joined = False
			for p in self.__processes:
				if(p.exitcode == None):
					continue
				p.join(timeout=5)
				if(p.exitcode != None):
					joined = True
					self.__processes.remove(p)
					break
			if(not joined):
				#print "THIS SHOULDNT HAPPEN"
				continue
			'''
			if(self.__outqueue.empty()):
				print "ERROR: Outqueue is empty after process join"
				ret = None
			else:
				ret = self.__outqueue.get()
			if(self.__errqueue.empty()):
				print "ERROR: Errqueue is empty after process join"
				err = "ERROR: QUEUES WERE EMPTY AFTER PROCESS JOIN"
			else:
				err = self.__errqueue.get()
			'''
		self.__last_wait = time.time()
	# FUNCTION END
	
	def wait_inner(self):
		# Retrieve all values in queue
		#self.__rets = []
		#self.__errs = []
		rets = []
		errs = []
		while(len(self.__processes) > 0):
			if(not self.__outqueue.empty() or not self.__errqueue.empty()):
				ret = self.__outqueue.get(block=True, timeout=None)
				err = self.__errqueue.get(block=True, timeout=None)
				rets.append(ret)
				errs.append(err)
				continue
			else:
				pass
			#	#continue
			joined = False
			for p in self.__processes:
				if(p.exitcode == None):
					continue
				p.join(timeout=5)
				if(p.exitcode != None):
					joined = True
					self.__processes.remove(p)
					if(not self.__outqueue.empty() or not self.__errqueue.empty()):
						ret = self.__outqueue.get(block=True, timeout=None)
						err = self.__errqueue.get(block=True, timeout=None)
						rets.append(ret)
						errs.append(err)
					break
			if(not joined):
				#print "THIS SHOULDNT HAPPEN"
				continue
			'''
			if(self.__outqueue.empty()):
				print "ERROR: Outqueue is empty after process join"
				continue
				ret = None
			else:
				ret = self.__outqueue.get()
			if(self.__errqueue.empty()):
				print "ERROR: Errqueue is empty after process join"
				err = "ERROR: QUEUES WERE EMPTY AFTER PROCESS JOIN"
			else:
				err = self.__errqueue.get()
			'''
			#rets.append(ret)
			#errs.append(err)
		
		# Join all processes
		#for p in self.__processes:
		#	p.join()
		
		
		# Add in intermediate waiting results
		rets = rets + self.__int_rets
		errs = errs + self.__int_errs
		
		# Return outs, errs list
		self.__rets = [x[1] for x in sorted(rets, key=lambda x: x[0])]
		self.__errs = [x[1] for x in sorted(errs, key=lambda x: x[0])]
		return self.__rets, self.__errs
	# FUNCTION END
	
	def run_multiple(self, func, inputs):
		self.__processes = []
		self.__outqueue = Queue()
		self.__errqueue = Queue()
		self.__lock = Lock()
		self.__counter = Value("i", 0)
		self.__running = Value("i", 0)
		self.__rets = []
		self.__errs = []
		self.__int_rets = []
		self.__int_errs = []
		self.__last_wait = time.time()
		
		# Automatically submit multiple jobs / return outputs
		for i in self.__progress_bar(inputs):
			self.run(func, i)
		ret, err = self.wait()
		#self.__rets = ret
		#self.__errs = err
		if(not len(ret) == len(err) and len(ret) == len(inputs)):
			print "ERROR: OUTPUTS NOT SAME SIZE AS INPUTS"
		return ret, err
	# FUNCTION END
	
# CLASS END

# Creates arg parser for a pipeline by reading a file describing the parameters. Extremely dangerous...
# - Accepts as parameters...
#   @filename - The filename containing the setup for the desired parameters
#   @description - The description of the pipeline
#   @max_line - The maximum number of characters to be printed per line in the help text
#   @set_default_params - A boolean flag indicating whether default parameter values should be used automatically when command line parameter not included.
#                         I set this to False because I prefer to provide the option of providing a text file instead of command line parameters and then allowing command line parameters to overrule the contents of the text file
#   create_default_namespace - A boolean flag indicating whether global variables for the parameters should be created.
#                              If set to false, a dictionary containing variables will be created / returned instead. Much safer.
#                              Note: This is extremely dangerous since it allows python execution to be run based on the contents of a text file. Extremely succeptible to malicious insertions.
#   Returns a pandas df representing the parameters and the arg parser created (and a dictionary with all parameters if "safe" version used)
# - Written for general utility purposes 02/05/18 by Shayne Wierbowski
def gen_parser(filename, description, gs, ls, max_line=60, set_default_params=False, create_default_namespace=True):
	parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
	if(not os.path.exists(filename)):
		print "ERROR: Could not identify parameter description file. Someone has seriously goofed."
	args = pd.read_csv(filename, sep="\t", keep_default_na=False, na_values="", quotechar="\'")
	
	def do(x):
		if(x == ""):
			return np.nan
		if(pd.isnull(x)):
			return x
		try:
			return eval(str(x), gs, ls)
		except (ValueError, SyntaxError, NameError):
			return x
	# FUNCTION END
	args["action"] = args["action"].map(lambda x: do(x))
	args["nargs"] = args["nargs"].map(lambda x: do(x))
	args["const"] = args["const"].map(lambda x: do(x))
	args["default"] = args["default"].map(lambda x: do(x))
	args["type"] = args["type"].map(lambda x: do(x))
	args["choices"] = args["choices"].map(lambda x: do(x))
	
	done = []
	def do(x, parser):
		def a(x):
			if(type(x) == str):
				#print "HEE"
				return "\"" + x + "\""
			if(type(x) == type):
				return re.findall("\'(.*)\'", str(x))[0]
			#print "BOO", x""
			return x
		# FUNCTION END
		
		if(str(x) in done):
			return
		done.append(str(x))
		command = "parser.add_argument("
		if(not pd.isnull(x["name"])):
			command += "\"-{0}\", ".format(x["name"])
		
		if(not pd.isnull(x["name2"])):
			command += "\"--{0}\", ".format(x["name2"])
		
		if(not pd.isnull(x["nargs"]) and not x["nargs"] in [1, "1"]):
			command += "nargs={0}, ".format(a(x["nargs"]))
		
		if(not pd.isnull(x["const"])):
			command += "const={0}, ".format(a(x["const"]))
		
		if(set_default_params and not pd.isnull(x["default"])):
			command += "default={0}, ".format(a(x["default"]))
		
		if(not pd.isnull(x["type"])):
			def str2bool(v):
				if v.lower() in ('yes', 'true', 't', 'y', '1'):
					return True
				elif v.lower() in ('no', 'false', 'f', 'n', '0'):
					return False
				else:
					raise argparse.ArgumentTypeError('Boolean value expected.')
			# FUNCTION END
			if(a(x["type"]) == "bool"):
				command += "type=str2bool, "
			else:
				command += "type={0}, ".format(a(x["type"]))
		
		if(type(x["choices"]) == list or not pd.isnull(x["choices"])):
			command += "choices={0}, ".format(a(x["choices"]))
		
		if(not pd.isnull(x["required"])):
			command += "required={0}, ".format(a(x["required"]))
		
		if(not pd.isnull(x["metavar"])):
			command += "metavar=\"<{0}>\", ".format(x["metavar"])
		
		if(not pd.isnull(x["dest"])):
			command += "dest={0}, ".format(a(x["dest"]))
		
		if(not pd.isnull(x["help"])):
			h = x["help"].split("\\n")
			h2 = []
			for l in h:
				l = l.split()
				if(len(l) == 0):
					h2.append("")
				while(len(l) > 0):
					tmp = ""
					while(len(l) > 0 and (len(tmp) == 0 or len(tmp) + len(l[0]) < max_line)):
						tmp += l.pop(0) + " "
					h2.append(tmp)
			h2.append("")
			h2.append("")
			command += "help=\"{0}\"".format("\\n".join(h2).replace("\"", "\\\""))
		
		command = command.strip(", ")
		
		command += ")"
		
		exec(command) in globals(), locals()
	# FUNCTION END
	args[list(args)].apply(lambda x: do(x, parser), axis=1)
	
	if(create_default_namespace):
		def do(x, y):
			def a(x):
					if(type(x) == str):
						return "\"" + x + "\""
					if(type(x) == type):
						return re.findall("\'(.*)\'", str(x))[0]
					return x
			#print("global {0}\n{0} = {1}".format(x, a(y)))
			exec("global {0}\n{0} = {1}".format(x, a(y))) in gs, ls
		# FUNCTION END
		args[["dest", "default"]].apply(lambda (x, y): do(x, y), axis=1)
		return args, parser
	
	return args, parser, {x[0]:x[1] for x in args[["dest", "default"]].values}
# FUNCTION END


# Profiles a python script using rkern's line_profiler (https://github.com/rkern/line_profiler) and reads it into a pandas df
# - Accepts as parameters...
#   @filename - The file to profile
#   @param_string - Any command line parameters to include when running the script for profiling (as a string, e.g. "-i in.txt -o out.txt")
#   @out - Location to save raw profile output to
#   Returns a pandas df representing the profiler output
#   NOTE: Susceptible to breakage since this script manually inserts "@Profile\n" before every function declaration using regex
#   NOTE: Does not work with processes
# - Written for general utility purposes 02/08/18 by Shayne Wierbowski
def line_profile(filename, param_string="", out=None):
	if(out == None):
		out = reserveTemp()
	
	# Read in code
	lines = open(filename, "r").read()
	
	# Place @Profile at the beginning of each function
	p = re.compile("[^\n\s]*def .*:")
	
	for m in sorted([m for m in p.finditer(lines)], key=lambda x: -1*x.start()):
		# Figure out indentation
		offset = lines[lines[:m.start()].rfind("\n"):m.start()].replace("\n", "")
		#print len(offset)
		# Insert @Profile + newline + proper indentation
		lines = lines[:m.start()] + "@profile\n" + offset + lines[m.start():]
	
	lines = ["import os\nline_profiler_starting_dir = os.getcwd()\n@profile\ndef line_profile_main():"] + ["\t" + l for l in lines.split("\n")] + ["# FUNCTION END\nline_profile_main()\nos.chdir(line_profiler_starting_dir)"]
	lines = "\n".join(lines)
	
	print lines
	
	# Write Profile-ready script to tmp location
	tmp = os.path.basename(reserveTemp())
	print tmp
	easyWriteLines(tmp, lines)
	
	
	# Call Profiler
	print call("kernprof -l {0} {1}".format(tmp, param_string))
	print out
	print call("python -m line_profiler {0}.lprof > {1}".format(tmp, out))
	
	
	
	myCleanupFile(tmp)
	myCleanupFile(tmp + ".lprof")
	
	return parse_line_profile(out)
# FUNCTION END

# Parses the raw output from rkern's line_profiler (https://github.com/rkern/line_profiler) and into a pandas df
# - Accepts as parameters...
#	@filename - The file to parse
#	Returns a pandas df representing the profiler output
# - Written for general utility purposes 02/08/18 by Shayne Wierbowski
def parse_line_profile(filename):
	lines = easyReadLines(filename)
	time_unit = float(lines.pop(0).split("unit: ")[-1].split(" s")[0])
	
	keep = []
	l = lines.pop(0)
	while(len(lines) > 0):
		# Get to start of entry
		while(l != ""):
			l = lines.pop(0)
		
		# Read Entry Info
		total_time = float(lines.pop(0).split("time: ")[-1].split(" s")[0])
		filename = lines.pop(0).split("File: ")[-1]
		function = lines.pop(0).split("Function: ")[-1]
		l = lines.pop(0)
		while(l[:6] != "Line #"):
			l = lines.pop(0)
		
		# Parse Header
		orig = l
		tmp = ""
		while(tmp != l):
			tmp = l
			l = l.replace("   ", "  ")
		cols = l.split("  ")
		#cols[3] = cols[3] + " "
		starts = [0] + [orig.find(cols[i]) + len(cols[i]) for i in range(len(cols))][:-1]
		#cols[3] = cols[3].strip()
		
		# Skip "=====" row
		lines.pop(0)
		
		def do(i):
			try:
				sub = l[starts[i]:starts[i+1]]
				if(sub.strip() == ""):
					return 0
				return float(sub)
			except IndexError:
				return l[starts[i]:]
			except:
				tmp = l.strip()
				tmp2 = ""
				while(tmp != tmp2):
					tmp2 = tmp
					tmp = tmp.replace("  ", " ")
				print i
				print starts
				print l
				print orig
				#raise
				print tmp.split(" ")[i]
				return float(tmp.split(" ")[i])
				
		# FUNCTION END
		
		l = lines.pop(0)
		rows = []
		while(l != ""):
			rows.append([do(i) for i in range(len(starts))])
			l = lines.pop(0)
		tmp = pd.DataFrame(rows, columns=cols)
		tmp["Total Time"] = total_time
		tmp["File"] = filename
		tmp["Function"] = function
		keep.append(tmp)
	df = pd.concat(keep)
	df["Time"] = df["Time"]*time_unit
	df["Per Hit"] = df["Per Hit"]*time_unit
	return df.ix[:, ["File", "Function", "Total Time"] + cols]
# FUNCTION END

def longest_substring(s1, s2, return_starts=True):
	match = SequenceMatcher(a=s1, b=s2).find_longest_match(0, len(s1), 0, len(s2))
	
	if(return_starts):
		return match.a, match.b, s1[match.a: match.a + match.size]
	return s1[match.a: match.a + match.size]
# FUNCTION END

# Copied from https://matplotlib.org/examples/color/colormaps_reference.html
def show_color_maps():
	"""
	==================
	Colormap reference
	==================
	
	Reference for colormaps included with Matplotlib.

	This reference example shows all colormaps included with Matplotlib. Note that
	any colormap listed here can be reversed by appending "_r" (e.g., "pink_r").
	These colormaps are divided into the following categories:
	
	Sequential:
		These colormaps are approximately monochromatic colormaps varying smoothly
		between two color tones---usually from low saturation (e.g. white) to high
		saturation (e.g. a bright blue). Sequential colormaps are ideal for
		representing most scientific data since they show a clear progression from
		low-to-high values.
	
	Diverging:
		These colormaps have a median value (usually light in color) and vary
		smoothly to two different color tones at high and low values. Diverging
		colormaps are ideal when your data has a median value that is significant
		(e.g.  0, such that positive and negative values are represented by
		different colors of the colormap).
	
	Qualitative:
		These colormaps vary rapidly in color. Qualitative colormaps are useful for
		choosing a set of discrete colors. For example::
	
			color_list = plt.cm.Set3(np.linspace(0, 1, 12))
	
		gives a list of RGB colors that are good for plotting a series of lines on
		a dark background.
	
	Miscellaneous:
		Colormaps that don't fit into the categories above.
	
	"""
	
	# Have colormaps separated into categories:
	# http://matplotlib.org/examples/color/colormaps_reference.html
	cmaps = [('Perceptually Uniform Sequential', [
					'viridis', 'plasma', 'inferno', 'magma']),
			 ('Sequential', [
					'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
					'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
					'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']),
			 ('Sequential (2)', [
					'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
					'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
					'hot', 'afmhot', 'gist_heat', 'copper']),
			 ('Diverging', [
					'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
					'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']),
			 ('Qualitative', [
					'Pastel1', 'Pastel2', 'Paired', 'Accent',
					'Dark2', 'Set1', 'Set2', 'Set3',
					'tab10', 'tab20', 'tab20b', 'tab20c']),
			 ('Miscellaneous', [
					'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
					'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
					'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])]
	
	nrows = max(len(cmap_list) for cmap_category, cmap_list in cmaps)
	gradient = np.linspace(0, 1, 256)
	gradient = np.vstack((gradient, gradient))
	
	def plot_color_gradients(cmap_category, cmap_list, nrows):
		fig, axes = plt.subplots(nrows=nrows)
		fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
		axes[0].set_title(cmap_category + ' colormaps', fontsize=14)
		
		for ax, name in zip(axes, cmap_list):
			ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
			pos = list(ax.get_position().bounds)
			x_text = pos[0] - 0.01
			y_text = pos[1] + pos[3]/2.
			fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)
		
		# Turn off *all* ticks & spines, not just the ones with colormaps.
		for ax in axes:
			ax.set_axis_off()
		
	for cmap_category, cmap_list in cmaps:
		try:
			plot_color_gradients(cmap_category, cmap_list, nrows)
		except KeyboardInterrupt:
			raise
		except:
			pass
	
	plt.show()
# FUNCTION END

# Modified Juan's usageStats to include amount free, totals, GB of MEM, and %CPU as a 100% 
def usageStats(n_samples=10, sample_wait=1, ngigs=252, ncores=144):
	# Get Juan's usageStats / Add GB of MEM Usage
	df = jfb_tools.usageStats(n_samples, sample_wait)
	df[('GB MEM', 'max')] = df[('%MEM', 'max')].map(lambda x: x*ngigs/100.0)
	df[('GB MEM', 'min')] = df[('%MEM', 'min')].map(lambda x: x*ngigs/100.0)
	df[('GB MEM', 'mean')] = df[('%MEM', 'mean')].map(lambda x: x*ngigs/100.0)
	
	# Add Amount Free to df
	free = ["FREE:", 
			100 - sum(df[('%MEM', 'max')]), 100 - sum(df[('%MEM', 'min')]), 100 - sum(df[('%MEM', 'mean')]),
			ncores - sum(df[('#CPUs', 'max')]), ncores - sum(df[('#CPUs', 'min')]), ncores - sum(df[('#CPUs', 'mean')]),
			100*ncores - sum(df[('%CPU', 'max')]), 100*ncores - sum(df[('%CPU', 'min')]), 100*ncores - sum(df[('%CPU', 'mean')]),
			ngigs - sum(df[('GB MEM', 'max')]), ngigs - sum(df[('GB MEM', 'min')]), ngigs - sum(df[('GB MEM', 'mean')])] 
	free = pd.DataFrame(free, ["USER"] + list(df)).T.set_index("USER")
	df = pd.concat([df, free])
	
	# Convert %CPU to a 100 max
	df[('%CPU', 'max')] = df[('%CPU', 'max')].map(lambda x: 100*x/sum(df[("%CPU", "max")]))
	df[('%CPU', 'min')] = df[('%CPU', 'min')].map(lambda x: 100*x/sum(df[("%CPU", "min")]))
	df[('%CPU', 'mean')] = df[('%CPU', 'mean')].map(lambda x: 100*x/sum(df[("%CPU", "mean")]))
	
	# Set Totals in df
	total = ["TOTAL:", 
			sum(df[('%MEM', 'max')]), sum(df[('%MEM', 'min')]), sum(df[('%MEM', 'mean')]),
			sum(df[('#CPUs', 'max')]), sum(df[('#CPUs', 'min')]), sum(df[('#CPUs', 'mean')]),
			sum(df[('%CPU', 'max')]), sum(df[('%CPU', 'min')]), sum(df[('%CPU', 'mean')]),
			sum(df[('GB MEM', 'max')]), sum(df[('GB MEM', 'min')]), sum(df[('GB MEM', 'mean')])] 
	total = pd.DataFrame(total, ["USER"] + list(df)).T.set_index("USER")
	
	df = pd.concat([df, total])
	df = df.loc[:, [('GB MEM', 'max'),
					('GB MEM', 'min'),
					('GB MEM', 'mean'),
					('%MEM', 'max'),
					('%MEM', 'min'),
					('%MEM', 'mean'),
					('#CPUs', 'max'),
					('#CPUs', 'min'),
					('#CPUs', 'mean'),
					('%CPU', 'max'),
					('%CPU', 'min'),
					('%CPU', 'mean')]]
	return df
# FUNCTION END

# Stolen from https://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale
# Generates colormaps from a sequence of (c1, c2, transition_midpoint)
def make_colormap(seq):
	'''Return a LinearSegmentedColormap
	seq: a sequence of floats and RGB-tuples. The floats should be increasing
	and in the interval (0,1).
	'''
	seq = [(None,) * 4, 0.0] + list(seq) + [1.0, (None,) * 4]
	cdict = {'red': [], 'green': [], 'blue': [], "alpha": []}
	for i, item in enumerate(seq):
		if isinstance(item, float):
			if(len(seq[i - 1]) == 3):
				r1, g1, b1 = seq[i - 1]
				a1 = 1
			else:
				r1, g1, b1, a1 = seq[i - 1]
			if(len(seq[i + 1]) == 3):
				r2, g2, b2 = seq[i + 1]
				a2 = 1
			else:
				r2, g2, b2, a2= seq[i + 1]
			cdict['red'].append([item, r1, r2])
			cdict['green'].append([item, g1, g2])
			cdict['blue'].append([item, b1, b2])
			cdict['alpha'].append([item, a1, a2])
	return matplotlib.colors.LinearSegmentedColormap('CustomMap', cdict)
# FUNCTION END

AA_3to1d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
			  'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
			  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
			  'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
# Converts 3 character AA symbol to 1 character AA symbol
def AA_3to1(seq, raise_errors=True):
	 if len(seq) % 3 != 0: 
		  if(raise_errors):
				raise ValueError('Input length should be a multiple of three')
		  else:
				return "*"
	 
	 return "".join([AA_3to1d.get(seq[3*i:3*i+3], "*") for i in range(len(seq)/3)])
# FUNCTION END

AA_1to3d = {v: k for k, v in AA_3to1d.iteritems()}
# Converts 1 character AA symbol to 3 character AA symbol
def AA_1to3(seq):
	return "".join([AA_1to3d[x] for x in list(seq)])
# FUNCTION END

# Reads in a PDB File as a DataFrame (2018, 07, 10)
def pdb2df(filename, whitelist=["ATOM", "HETATM"]):
	# Try accessing file using mjm_tools open_pdb
	if(not os.path.exists(filename)):
		fh = mjm_tools.open_pdb(filename)
		
		lines = [x.strip() for x in fh.readlines()]
	else:
		# Handle unzipping gzipped files
		if(filename.split(".")[-1] == "gz"):
			tmp = reserveTemp()
			call("gunzip -c {0} > {1}".format(filename, tmp))
			filename = tmp
	
		# Read PDB File
		lines = easyReadLines(filename)
	 
	# Filter lines to only interesting records saving header and tailer
	tmp = [i for i in range(len(lines)) if lines[i][:6].strip() in whitelist]
	header = lines[:min(tmp)]
	tailer = lines[max(tmp) + 1:]
	
	model_ends = [i for i in range(len(lines)) if lines[i][:6] == "ENDMDL"]
	if(len(model_ends) != 0):
		lines = lines[:min(model_ends)]
	lines = [x for x in lines if x[:6].strip() in whitelist]
	 
	# Check that lines are formatted as expected
	if(len(set([len(x) for x in lines]).difference(set([81, 80, 78, 26]))) >= 1):
		#        |     1 | |    | |  2| | |    | |  3|        | 4      |   5    |     6|      |   7  |    |  | 8
		#    1234|5678901|2|3456|7|890|1|2|3456|7|890|12345678|90123456|78901234|567890|123456|789012|3456|78|90
		#    ATOM|    500| |CA  |A|LEU| |J| 700| |   | -13.123| -20.123| -30.123|     N|     0|      |0   | 0|\n
		print "WARNING: PDB lines do not match expectation (length 81, 80, 78, or 26 characters). Observed length", set([len(x) for x in lines]).difference(set([81, 78, 26]))
		print "I do not know why this is happening, maybe consult the PDB format summary again (http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html)"
	
	# Extract columns / correct type
	names = ["Data Type", "Atom ID", "Atom Name", "Alternate", "Residue Name", "Chain", "Residue ID", "Insertions?", "X", "Y", "Z", "Occupancy", "Temperature Factor", "Segment ID", "Element", "Charge"]
	types = [str, int, str, str, str, str, int, str, float, float, float, float, float, str, str, str, str]
	points = [0, 6, 12, 16, 17, 21, 22, 26, 30, 38, 46, 54, 60, 72, 76, 78, 81]
	def do(x, t):
		if(x == ""):
			return np.nan
		else:
			return t(x)
	lines = [[do(x[points[i]:points[i+1]].strip(), types[i]) for i in range(len(points) - 1)] for x in lines]
	
	# Create DataFrame
	df = pd.DataFrame(lines, columns=names)
	
	# Convert Resn to 1 Character
	df["Residue Sym"] = df["Residue Name"].map(lambda x: AA_3to1(x, raise_errors=False))
	
	# Set header / tailer attributes
	df.header = header
	df.tailer = tailer
	
	# Return
	return df
# FUNCTION END

# Writes a DF representation of a PDB to a PDB file (2019, 01, 15)
# SOMETHING IS VERY WRONG WITH THIS
def df2pdb(filename, df, columns=["Data Type", "Atom ID", "Atom Name", "Alternate", "Residue Name", "Chain", "Residue ID", "Insertions?", "X", "Y", "Z", "Occupancy", "Temperature Factor", "Segment ID", "Element"]):
	out = open(filename, "w+")
	indices = [list(df).index(c) for c in columns]
	s = ""
	for i in range(len(df)):
		atom, atmi, atmn, alt, resn, chain, resi, insertions, x, y, z, occupancy, temp, segid, esym, = df.iloc[i, indices].replace(np.nan, "")
		#out.write("{atom:>4}  {atmi:>5} {atmn:<4}{alt:<1}{resn:>3} {chain}{resi:>4.4}{insertions:>1}   {x:>8.8}{y:>8.8}{z:>8.8}{occupancy:>6}{temp:>6}     {segid:<4}{esym:>2}\n".format(atom=atom, atmi=atmi, atmn=atmn, alt=alt, resn=resn, chain=chain, resi=str(resi), insertions=insertions, x=str(x), y=str(y), z=str(z), occupancy=occupancy, temp=temp, segid=segid, esym=esym))
		
		# Changed spacing between ATMI and ATMN column 2019_03_15
		#out.write("{atom:<6}{atmi:>5} {atmn:<4}{alt:<1}{resn:>3} {chain}{resi:>4.4}{insertions:>1}   {x:>8.8}{y:>8.8}{z:>8.8}{occupancy:>6}{temp:>6}      {segid:<4}{esym:>2}\n".format(atom=atom, atmi=atmi, atmn=atmn, alt=alt, resn=resn, chain=chain, resi=str(resi), insertions=insertions, x=str(x), y=str(y), z=str(z), occupancy=occupancy, temp=temp, segid=segid, esym=esym))
		#out.write("{atom:<6}{atmi:>5}  {atmn:<3}{alt:<1}{resn:>3} {chain}{resi:>4.4}{insertions:>1}   {x:>8.8}{y:>8.8}{z:>8.8}{occupancy:>6}{temp:>6}      {segid:<4}{esym:>2}\n".format(atom=atom, atmi=atmi, atmn=atmn, alt=alt, resn=resn, chain=chain, resi=str(resi), insertions=insertions, x=str(x), y=str(y), z=str(z), occupancy=occupancy, temp=temp, segid=segid, esym=esym))
		# Swap to writting all output at once ~2020_02_20
		#s += "{atom:<6}{atmi:>5}  {atmn:<3}{alt:<1}{resn:>3} {chain}{resi:>4.4}{insertions:>1}   {x:>8.8}{y:>8.8}{z:>8.8}{occupancy:>6}{temp:>6}      {segid:<4}{esym:>2}\n".format(atom=atom, atmi=atmi, atmn=atmn, alt=alt, resn=resn, chain=chain, resi=str(resi), insertions=insertions, x=str(x), y=str(y), z=str(z), occupancy=occupancy, temp=temp, segid=segid, esym=esym)
		# Undid spacing change between ATMI and ATMN 2020_02_24 (Seems like it was NOT correct, but why did I make that change in the first place?)
		s += "{atom:<6}{atmi:>5} {atmn:<4}{alt:<1}{resn:>3} {chain}{resi:>4.4}{insertions:>1}    {x:>7.7} {y:>7.7} {z:>7.7}{occupancy:>6}{temp:>6}      {segid:<4}{esym:>2}\n".format(atom=atom, atmi=atmi, atmn=atmn, alt=alt, resn=resn, chain=chain, resi=str(resi), insertions=insertions, x=str(x)[:7], y=str(y)[:7], z=str(z)[:7], occupancy=occupancy, temp=temp, segid=segid, esym=esym)
	out.write(s)
	out.close()
# FUNCTION END

# Displays a pre-generated CloneSeq Plot (July 13, 2018)
def plotCloneSeq(project, orf, pos):
	# Find index of image containing the position
	i = 0
	while((i+1)*100 < pos):
		i += 1
	
	# Load image
	try:
		img = matplotlib.image.imread("[REDACTED_PATH]/NAS/Clone_Seq/{0}/Plots/{1}/{1}_{2}.png".format(project, orf, i))
	except:
		print "ERROR: Image {0} could not be found".format("[REDACTED_PATH]/NAS/Clone_Seq/{0}/Plots/{1}/{1}_{2}.png".format(project, orf, i))
		return None
	# Trim Excess Whitespace
	img = img[400:-400, 300:-300, :]
	
	# Prep figure to be correct size (See https://stackoverflow.com/questions/8056458/display-image-with-a-zoom-1-with-matplotlib-imshow-how-to)
	dpi = 80
	margin = 0.00 # (5% of the width/height of the figure...)
	xpixels, ypixels = img.shape[0], img.shape[1]

	# Make a figure big enough to accomodate an axis of xpixels by ypixels
	# as well as the ticklabels, etc...
	figsize = (1 + margin) * ypixels / dpi, (1 + margin) * xpixels / dpi

	fig = plt.figure(figsize=figsize, dpi=dpi)
	# Make the axis the right size...
	ax = fig.add_axes([margin, margin, 1 - 2*margin, 1 - 2*margin])
	
	# Turn off axis for the plot contianing the image
	ax.axis('off')
	
	ax.imshow(img, interpolation='none')
	
	return fig
# FUNCTION END

# Aligment Functions


from Bio import pairwise2
from Bio.SubsMat import MatrixInfo

# Alignment Helper Functions (Copied from Rice Eclair Predictions Notebook July 17, 2018)
def alignPident(align1, align2, useShorter=True):
	align1 = align1.upper()
	align2 = align2.upper()
	
	# Convert alignments to array, sum up identical AA
	identical = sum(np.array((list(align1)) == np.array(list(align2)))&(np.array(list(align1)) != "-"))
	# Figure out what to calculate percentage based off of
	# Shorter Sequence
	# Longer Sequence
	# Aligned Portion?
	# Total Alignment Length?
	if useShorter:
		totalLength = min([sum(np.array(list(align1)) != "-"), sum(np.array(list(align2)) != "-")])
	else:
		totalLength = max([sum(np.array(list(align1)) != "-"), sum(np.array(list(align2)) != "-")])

	# Return Percentage
	return float(identical)/totalLength
# FUNCTION END

# Alignment Helper Functions (Copied from Rice Eclair Predictions Notebook July 17, 2018)
def alignPositives(align1, align2, matrix=MatrixInfo.blosum62, useShorter=True):
	# Function for reading blosum score (since matrix uses non-symetrical (x, y) tuples as keys)
	def getBlosumScore(x, y):
		try:
			return matrix[(x, y)]
		except KeyError:
			return matrix[(y, x)]
	# FUNCTION END

	# Calculate Blosum score over alilgnment
	blosumScores = [getBlosumScore(align1[i], align2[i]) for i in range(len(align1)) if not (align1[i] == "-" or align2[i] == "-")]

	# Sum Positives
	positive = sum(np.array(blosumScores) > 0)

	# Figure out what to calculate percentage based off of
	# Shorter Sequence
	# Longer Sequence
	# Aligned Portion?
	# Total Alignment Length?
	if useShorter:
		totalLength = min([sum(np.array(list(align1)) != "-"), sum(np.array(list(align2)) != "-")])
	else:
		totalLength = max([sum(np.array(list(align1)) != "-"), sum(np.array(list(align2)) != "-")])

	# Return Percentage
	return float(positive)/totalLength
# FUNCTION END

# Alignment Helper Functions (Copied from Rice Eclair Predictions Notebook July 17, 2018)
def alignCoverage(align1, align2, matrix=MatrixInfo.blosum62, useShorter=True):
	# Sum up non-gapped regions of alignment
	non_gaps = len(align1) - sum([align1[i] == "-" or align2[i] == "-" for i in range(len(align1))])

	# Divide by sequence length of align1
	return float(non_gaps)/len(align1.replace("-", ""))
# FUNCTION END

# Alignment Helper Functions (Copied from Rice Eclair Predictions Notebook July 17, 2018)
def getBlosumScore(x, y, matrix=MatrixInfo.blosum62):
	if((x, y) in matrix):
		return matrix[(x, y)]
	else:
		return matrix[(y, x)]
# FUNCTION END

def getDNASubMatrix():
	df = pd.read_csv("[REDACTED_PATH]/resources/substitution_matrices/EDNAFULL.sub", sep="\t", index_col=0)
	header = list(df)
	matrix = dict()
	for i in df.index:
		for j in header:
				matrix[(i, j)] = df.loc[i, j]
				matrix[(i.lower(), j)] = df.loc[i, j]
				matrix[(i, j.lower())] = df.loc[i, j]
				matrix[(i.lower(), j.lower())] = df.loc[i, j]
	return matrix
# FUNCTION END

# Alignment Helper Functions (Copied from Rice Eclair Predictions Notebook July 17, 2018)
def NWSeqAlignment(s1, s2, matrix=None, useShorter=True, best=False, show_align=False, genPosMap=False):
	keep = None
	
	kind = "prot"
	if(matrix == None):
		matrix = MatrixInfo.blosum62
		if(len(set(list(s1.upper()) + list(s2.upper())).difference(set(["A", "C", "G", "T", "N"]))) == 0):
			kind = "nucl"
			matrix = getDNASubMatrix()
	
	# Added 2019_10_08 to resolve the gapped end of alignment problem.
	# Effectively prioritized fewest gap closes to break ties where possible.
	#
	# e.g.
	#
	#		AAAAAG--------G
	#		AAAAAGGCCCCCCCG
	#
	# instead of...
	#
	#		AAAAAGG--------
	#		AAAAAGGCCCCCCCG
	#
	if(kind == "prot"):
		align = pairwise2.align.globalds(s1, s2, matrix, -10, -0.5)
	else:
		align = pairwise2.align.globalds(s1, s2, matrix, -400, 0)
	align = sorted(align, key=lambda x: (-x[2], len(re.findall("-[^-]", x[0])) + len(re.findall("-[^-]", x[1]))))
	for a in align:
		a1, a2, score, begin, end = a
		r = dict([("Align1", a1), ("Align2", a2), ("Score", score), ("Begin", begin), ("End", end), ("Pident", alignPident(a1, a2, useShorter=useShorter)), ("Positives", alignPositives(a1, a2, matrix=matrix, useShorter=useShorter)), ("Coverage1", alignCoverage(a1, a2, matrix=matrix, useShorter=useShorter)), ("Coverage2", alignCoverage(a2, a1, matrix=matrix, useShorter=useShorter))])
		align_string = pd.DataFrame([list(r["Align1"]), list(r["Align2"])]).T
		def do(x, y):
			if(x.upper() == y.upper()):
				return "|"
			if(x == "-" or y == "-"):
				return " "
			else:
				return ["-", "+"][getBlosumScore(x, y, matrix=matrix) > 0]
		# FUNCTION END
		align_string["C"] = align_string[[0, 1]].apply(lambda (x, y): do(x, y), axis=1)
		align_string = "".join(align_string["C"].values)
		r["Alignment"] = align_string
		
		if(keep == None or r["Pident"] >= keep["Pident"]):
			keep = r
		if(not best):
			break
	
	if(show_align):
		alignPrint(keep)
	if(genPosMap):
		return keep, alignPosMap(keep)
	else:
		return keep
# FUNCTION END


# Gets the time since last modification to a file (July 19, 2018)
def time_since_mod(f):
	return time.time() - os.path.getmtime(f)
# FUNCTION END

# Generate Posiiton Map from alignment info (April 22, 2019)
def alignPosMap(align1, align2=None, alignment=None, name1="Align1", name2="Align2", zero_index=False):
	if(align2 == None and type(align1) == dict):
		aln = align1
		align1 = aln["Align1"]
		align2 = aln["Align2"]
		alignment = aln["Alignment"]
	elif align2 == None:
		raise ValueError

	pos_map = []
	
	qpos = 0
	rpos = 0
	
	if(zero_index):
		qpos -= 1
		rpos -= 1
	
	for i in range(len(align1)):
		if(align1[i] == "-" and (align2[i] == "-")):
			continue
		elif(align1[i] == "-" and (align2[i] != "-")):
			rpos += 1
			pos_map.append([-1, "", rpos, align2[i]])
		elif(align1[i] != "-" and (align2[i] == "-")):
			qpos += 1
			pos_map.append([qpos, align1[i], -1, ""])
		else:
			rpos += 1
			qpos += 1
			pos_map.append([qpos, align1[i], rpos, align2[i]])
	
	# 2019_08_05 - This used to have "_Seq2AA" as the final column.
	# I think that was a mistake, but could have been written into some of my old code?
	pos_map = pd.DataFrame(pos_map, columns=[name1 + "_Pos", name1 + "_AA", name2 + "_Pos", name2 + "_AA"])
	
	return pos_map
# FUNCTION END

# Pretty Prints an alignemnt (July 31, 2018)
def alignPrint(align1, align2=None, alignment=None, name1="Align1", name2="Align2", rows=60, spacing=10, display=True):
	if(align2 == None and type(align1) == dict):
		aln = align1
		align1 = aln["Align1"]
		align2 = aln["Align2"]
		alignment = aln["Alignment"]
	elif align2 == None:
		raise ValueError
	
	def space_string(string, length):
		return ' '.join(string[i:i+length] for i in range(0, len(string), length))
	# FUNCTION END
	
	l1Base = name1 + ": "
	l2Base = name2 + ": "
	if(len(l1Base) > len(l2Base)):
		l2Base += " "*(len(l1Base) - len(l2Base))
	else:
		l1Base += " "*(len(l2Base) - len(l1Base))
	l3Base = " "*len(l1Base)
	
	lines = []
	num_size = len(str(len(alignment))) + 2
	
	pos1 = 1
	pos2 = 1
	
	for i in range(0, len(alignment), rows):
		startNum1 = " "*(num_size - 1 - len(str(pos1))) + str(pos1) + " "
		pos1 += len(align1[i:i+rows].replace("-", ""))
		endNum1 = " " + str(pos1 - 1) + " "*(num_size - 1 - len(str(pos1 - 1)))
		
		startNum2 = " "*(num_size - 1 - len(str(pos2))) + str(pos2) + " "
		pos2 += len(align2[i:i+rows].replace("-", ""))
		endNum2 = " " + str(pos2 - 1) + " "*(num_size - 1 - len(str(pos2 - 1)))
		
		l1 = l1Base + startNum1 + space_string(align1[i:i+rows], spacing) + endNum1
		l2 = l2Base + startNum2 + space_string(align2[i:i+rows], spacing) + endNum2
		l3 = l3Base + " "*num_size + space_string(alignment[i:i+rows], spacing) + " "*num_size
		
		lines.append(l1)
		lines.append(l3)
		lines.append(l2)
		lines.append("")
	
	if(display):
		print("\n".join(lines))
	else:
		return "\n".join(lines)
# FUNCTION END

# Renames a column in a DF
def rename_column(df, orig, new):
	if(orig in df.columns):
		cols = list(df.columns)
		cols[cols.index(orig)] = new
		df.columns = cols
# FUNCTION END

def odds_ratio(exposure_mask, case_mask, CI_value=0.05, two_sided=True, log_odds=False, verbose=False, long_output=False, expose_label="Exposed", case_label="Case", error="CI"):
	# Catch all for various types formats I've provided the masks in before (first case is for direct binary arrays)
	if(type(exposure_mask) == bitarray):
		a = float((exposure_mask & case_mask).count())
		b = float((exposure_mask & ~case_mask).count())
		c = float((~exposure_mask & case_mask).count())
		d = float((~exposure_mask & ~case_mask).count())
	else:
		# Otherwise they are probably numpy arrays which could either be Booleans or Ints
		# If they are Ints I want to make sure they actually only contain 0 and 1 (otherwise the input is trying to represent
		# a concept not supported by the code (e.g. if you tried to report there are multiple mutations at one location)
		d1 = set(exposure_mask).difference(set([0, 1, 0.0, 1.0, True, False]))
		d2 = set(case_mask).difference(set([0, 1, 0.0, 1.0, True, False]))
		if(d1 or d2):
			print "WARNING: Masks are not formatted as expected (do not contain just 0 and 1)"
			print set(exposure_mask)
			print set(case_mask)
			print "END"
		a = float(sum(exposure_mask*case_mask))
		b = float(sum(exposure_mask*(1 - case_mask)))
		c = float(sum((1 - exposure_mask)*case_mask))
		d = float(sum((1 - exposure_mask)*(1 - case_mask)))
	
	# Print the contingency table if desired
	if(verbose):
		tmp = pd.DataFrame([[a, b], [c, d]], columns=[case_label, "Non-" + case_label], index=[expose_label, "Non-" + expose_label])
		print tmp
	
	# Cannot calculate a valid OR if there are any 0 boxes
	if(0 in [a, b, c, d]):
		if(long_output):
			return np.nan, np.nan, np.nan, np.nan, a, b, c, d
		else:
			return np.nan, np.nan, np.nan, np.nan
	
	# Calculate OR / p-value
	OR = (a/c)/(b/d)
	p_value = 1 - (st.norm.cdf(abs(np.log(OR)) / np.sqrt(1/a + 1/b + 1/c + 1/d)))
	
	# Make adjustment for one-sided vs. two-sided test. Practically speaking we almost always want to use two_sided unless
	# we have a strong prior expectation that the enrichment will go one way verses the other.
	if(two_sided):
		CI_value = CI_value / 2.0
		p_value = p_value*2
	else:
		pass
	z_score = st.norm.ppf((1 - CI_value), )
	
	if(error == "CI"):
		upperCI = np.exp(np.log(OR) + z_score*np.sqrt(1/a + 1/b + 1/c + 1/d))
		lowerCI = np.exp(np.log(OR) - z_score*np.sqrt(1/a + 1/b + 1/c + 1/d))
	elif(error == "SE"):
		upperCI = np.exp(np.log(OR) + np.sqrt(1/a + 1/b + 1/c + 1/d))
		lowerCI = np.exp(np.log(OR) - np.sqrt(1/a + 1/b + 1/c + 1/d))		
	
	# Decide how much output to return (i.e. do I want the raw a, b, c, d values back)
	# and whether to conver the OR to log2(OR) (log2 transofmration is desirable because it makes
	# the magnitudes for enrichment vs. depletion symetrical (e.g. 1, 0.5, and 2 in OR space) vs. (0, -1, and 1 in logOR space). 
	if(long_output):
		if(log_odds):
			return np.log2(OR), np.log2(upperCI), np.log2(lowerCI), p_value, a, b, c, d
		else:
			return OR, upperCI, lowerCI, p_value, a, b, c, d
	else:
		if(log_odds):
			return np.log2(OR), np.log2(upperCI), np.log2(lowerCI), p_value
		else:
			return OR, upperCI, lowerCI, p_value
# FUNCTION END

# For creating Shayne readable versions of Scipy Linkage Matrices (2019_04_36)
def linkage2readable(linkage):
	 # Convert to DF
	 df = pd.DataFrame(linkage, columns=["Source_A", "Source_B", "Height", "Node_Size"])
	 
	 # Figure Out Number of Leaves
	 # Does not work because not all leaves are merged with other leaves
	 #offset = df[df["Node_Size"] == 2][["Source_A", "Source_B"]].max().max()
	 
	 offset = max(df[["Source_A", "Source_B"]].values[-1]) + 2 - len(df)
	 
	 # Add ID of newly created node
	 df["Merged_ID"] = [offset + x for x in range(len(df))]
	 
	 node2children = defaultdict(set)
	 for ida, idb, dist, size, idm in df.values:
		  node2children[idm] = node2children[ida].union(node2children[idb]).union(set([ida, idb]))
	 df["Children"] = df["Merged_ID"].map(lambda x: node2children[x])
	 df["Leaf_Children"] = df["Merged_ID"].map(lambda x: set([a for a in node2children[x] if a < offset]))
	 
	 return df
# FUNCTION END


# Shayne friendly wrapper for shitty primer3-py IO
class PrimerDesign():
	
	# Parameters controling primer design
	global_args = {# Basic Primer Size / TM / GC params (Set based on defaults in http://primer.yulab.org/)
	               'PRIMER_MIN_SIZE': 16,
	               'PRIMER_OPT_SIZE': 20,
	               'PRIMER_MAX_SIZE': 28,
	               'PRIMER_MIN_TM': 55.0,
	               'PRIMER_OPT_TM': 60.0,
	               'PRIMER_MAX_TM': 65.0,
	               'PRIMER_MIN_GC': 30.0,
	               'PRIMER_MAX_GC': 80.0,
	               
	               # What kind / how many primers to pick
	               'PRIMER_PICK_LEFT_PRIMER': 1,
	               'PRIMER_PICK_RIGHT_PRIMER': 1,
	               'PRIMER_PICK_INTERNAL_OLIGO': 0,
	               'PRIMER_NUM_RETURN': 5,
	               
	               # Valid product ranges for primers (Set to primer3-py internal default. Should be modified on a per-use basis)
	               'PRIMER_PRODUCT_SIZE_RANGE': [[100,300]],
	               'PRIMER_PRODUCT_OPT_SIZE': 0,
	               'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': 0,
	               'PRIMER_PAIR_WT_PRODUCT_SIZE_GT': 0,
	               
	               # These were just in the example I started from (No clue what they do?)
	               'PRIMER_INTERNAL_MAX_SELF_END': 8,
	               'PRIMER_MAX_POLY_X': 100,
	               'PRIMER_INTERNAL_MAX_POLY_X': 100,
	               'PRIMER_SALT_MONOVALENT': 50.0,
	               'PRIMER_DNA_CONC': 50.0,
	               'PRIMER_MAX_NS_ACCEPTED': 0,
	               'PRIMER_MAX_SELF_ANY': 12,
	               'PRIMER_MAX_SELF_END': 8,
	               'PRIMER_PAIR_MAX_COMPL_ANY': 12,
	               'PRIMER_PAIR_MAX_COMPL_END': 8,
	               
	               # Plus there are WAY more options here (http://primer3.ut.ee/primer3web_help.htm)
	              }
	
	# Not sure what these are used for
	misprime_lib = None
	mishyb_lib = None
	
	def __init__(self, primer_size_min=16, primer_size_max=28, primer_size_opt=20, TM_min=55.0, TM_max=65.0, TM_opt=60.0, GC_min=30.0, GC_mac=80.0, n_ret=5, **kwargs):
		# Save the constructor parameters
		self.primer_size_min = primer_size_min
		self.primer_size_max = primer_size_max
		self.primer_size_opt = primer_size_opt
		self.TM_min = TM_min
		self.TM_max = TM_max
		self.TM_opt = TM_opt
		self.GC_min = GC_min
		self.GC_mac = GC_mac
		self.n_ret = n_ret
		
		# Store any parameters the user explicitly adds
		self.explicit_params = kwargs
		
		# Update the global_args dictionary
		self.__update_internal_params__()
	# FUNCTION END
	
	def __update_internal_params__(self):
		# Update the internal parameters provided in the constructor
		update_dict1 = {'PRIMER_MIN_SIZE': self.primer_size_min,
		                'PRIMER_OPT_SIZE': self.primer_size_opt,
		                'PRIMER_MAX_SIZE': self.primer_size_max,
		                'PRIMER_MIN_TM': self.TM_min,
		                'PRIMER_OPT_TM': self.TM_opt,
		                'PRIMER_MAX_TM': self.TM_max,
		                'PRIMER_MIN_GC': self.GC_min,
		                'PRIMER_MAX_GC': self.GC_mac,
		                'PRIMER_NUM_RETURN': self.n_ret,
		               }
		self.global_args.update(update_dict1)
		
		# Update any parameters the user adds to the constructor
		self.global_args.update(self.explicit_params)
		
		# Update these parameters in the actual primer3 module
		primer3.setP3Globals(self.global_args)
	# FUNCTION END
	
	def design_primers(self, seqs, labels=None, included_regions=None, excluded_regions=None, n_ret=None, header=['Sequence_ID', 'Rank', 'Fwd_Seq', 'Fwd_Start', 'Fwd_End', 'Fwd_Size', 'Fwd_TM', 'Fwd_GC', 'Fwd_Hairpin', 'Rev_Seq', 'Rev_Start', 'Rev_End', 'Rev_Size', 'Rev_TM', 'Rev_GC', 'Rev_Hairpin', 'Overall_Penalty', 'Product_Size', 'Product_Seq']):
		# Make sure global params are up to date
		self.__update_internal_params__()
		
		# Update formatting on all inputs to handle list or non-list
		if(not is_iterable(seqs, exclude=[str])):
			seqs = [seqs]
		
		if(labels == None):
			labels = range(len(seqs))
		if(not is_iterable(labels, exclude=[str])):
			labels = itertools.cycle([labels])
		
		if(not is_iterable(included_regions)):
			included_regions = itertools.cycle([included_regions])
		
		if(not is_iterable(excluded_regions)):
			excluded_regions = itertools.cycle([excluded_regions])
		
		
		# Update max number of primer returned
		if(n_ret != None):
			primer3.setP3Globals({"PRIMER_NUM_RETURN":n_ret})
		
		
		headers_verbose = ["Sequence_ID", "Rank", "Fwd_Seq", "Fwd_Start", "Fwd_End", "Fwd_Size", "Fwd_TM", "Fwd_GC", "Fwd_Hairpin", "Fwd_End_Stability", "Fwd_Self_Comp_Any", "Fwd_Self_Comp_3p", "Fwd_Penalty",
		                   "Rev_Seq", "Rev_Start", "Rev_End", "Rev_Size", "Rev_TM", "Rev_GC", "Rev_Hairpin", "Rev_End_Stability", "Rev_Self_Comp_Any", "Rev_Self_Comp_3p", "Rev_Penalty",
		                   "Intern_Seq", "Intern_Start", "Intern_End", "Intern_Size", "Intern_TM", "Intern_GC", "Intern_Hairpin", "Intern_End_Stability", "Intern_Self_Comp_Any", "Intern_Self_Comp_3p", "Intern_Penalty",
		                   "Overall_Penalty", "Product_Size", "Product_Seq"]
		rows = []
		for s, l, ir, er in zip(seqs, labels, included_regions, excluded_regions):
			res = primer3.bindings.designPrimers(seq_args={"SEQUENCE_ID":l, "SEQUENCE_TEMPLATE":s, "SEQUENCE_INCLUDED":ir, "SEQUENCE_EXCLUDED":er})
			
			processed_keys = set()
			
			total_hits = len(set([int(key.split("_")[2]) for key in res.keys() if len(key.split("_")) >= 4 and not "NUM" in key]))
			
			for i in range(total_hits):
				def rename(key):
					begin, end = key.split("_{0}".format(i))
					begin = begin.replace("PRIMER_", "")
					begin = begin.replace("LEFT", "Fwd")
					begin = begin.replace("RIGHT", "Rev")
					begin = begin.replace("INTERNAL", "Intern")
					
					end_sub_dict = {"_SEQUENCE":"_Seq",
					                "":"_Start",
					                "_TM":"_TM",
					                "_GC_PERCENT":"_GC",
					                "_HAIRPIN_TH":"_Hairpin",
					                "_END_STABILITY":"_End_Stability",
					                "_SELF_ANY_TH":"_Self_Comp_Any",
					                "_SELF_END_TH":"_Self_Comp_3p",
					                "_PENALTY":"_Penalty",
					                "_PRODUCT_SIZE":"_Product_Size"}
					
					if(end in end_sub_dict):
						end = end_sub_dict[end]
					
					new_key = begin + end
					
					if(new_key == "PAIR_Penalty"):
						new_key = "Overall_Penalty"
					if(new_key == "PAIR_Product_Size"):
						new_key = "Product_Size"
					
					return new_key
				# FUNCITON END
				sub_res = defaultdict(lambda: np.nan)
				sub_res.update({rename(k):v for k, v in res.iteritems() if str(i) in k})
				
				try:
					sub_res["Fwd_End"] = sub_res["Fwd_Start"][0] + sub_res["Fwd_Start"][1]
					sub_res["Fwd_Size"] = sub_res["Fwd_Start"][1]
					sub_res["Fwd_Start"] = sub_res["Fwd_Start"][0]
				except:
					pass
				
				try:
					sub_res["Rev_End"] = sub_res["Rev_Start"][0] - sub_res["Rev_Start"][1] + 1
					sub_res["Rev_Size"] = sub_res["Rev_Start"][1]
					sub_res["Rev_Start"] = sub_res["Rev_Start"][0] + 1
				except:
					pass
				
				try:
					sub_res["Intern_End"] = sub_res["Intern_Start"][0] - sub_res["Intern_Start"][1]
					sub_res["Intern_Size"] = sub_res["Intern_Start"][1]
					sub_res["Intern_Start"] = sub_res["Intern_Start"][0]
				except:
					pass
				
				sub_res["Rank"] = i+1
				sub_res["Sequence_ID"] = l
				sub_res["Product_Seq"] = s[sub_res["Fwd_Start"]:sub_res["Rev_Start"]]
				
				s1 = set(list(sub_res))
				s2 = set(headers_verbose)
				
				#print s1.intersection(s2)
				#print
				#print sorted(s1.difference(s2))
				#print
				#print sorted(s2.difference(s1))
				#
				#1/0
				rows.append([sub_res[x] for x in headers_verbose])
		
		df = pd.DataFrame(rows, columns=headers_verbose)
		
		return df[header]
	# FUNCTION END
# CLASS END

def dated_file(basename):
	now = datetime.datetime.now()
	return "{0:04}_{1:02}_{2:02}_{3}".format(now.year, now.month, now.day, basename)
# FUNCTION END

# Calculates distance RMSD between two pdb files or sets of np arrays of coordinates
# NOTE: dRMSD is based on a comparison of the pairwise distances within pdb1 compared
#       to those within pdb2, and DOES NOT actually perform an alignment between the
#       structures.
# The onus is on you to make sure that the provided PDB files or arrays of coordinates
# have a 1-1 coorespondence.
def calculate_dRMSD(pdb1, pdb2, chain1=None, chain2=None, atom_whitelist=["N", "CA", "C"], resi_blacklist=["*"], check_match=True, lower_cuttoff=None, upper_cuttoff=None):
	# For case where inputs are provided as two pdb files
	if(type(pdb1) == str and type(pdb2) == str and os.path.exists(pdb1) and os.path.exists(pdb2)):
		# Read in PDB, filter by chain, atom whitelist and resi blacklist
		# Default should select only protein and backbone (but drops any non-conventional AA)
		pdb1 = my.pdb2df(pdb1)
		if(not chain1 == None):
			pdb1 = pdb1[pdb1["Chain"].map(lambda x: x == chain1)]
		pdb1 = pdb1[pdb1["Atom Name"].map(lambda x: x in atom_whitelist)]
		pdb1 = pdb1[pdb1["Residue Sym"].map(lambda x: not x in resi_blacklist)]
		pdb1 = pdb1.drop_duplicates(["Atom Name", "Residue Name", "Chain", "Residue ID"])
		
		# Read in PDB, filter by chain, atom whitelist and resi blacklist
		# Default should select only protein and backbone (but drops any non-conventional AA)
		pdb2 = my.pdb2df(pdb2)
		if(not chain2 == None):
			pdb2 = pdb2[pdb2["Chain"].map(lambda x: x == chain2)]
		pdb2 = pdb2[pdb2["Atom Name"].map(lambda x: x in atom_whitelist)]
		pdb2 = pdb2[pdb2["Residue Sym"].map(lambda x: not x in resi_blacklist)]
		pdb2 = pdb2.drop_duplicates(["Atom Name", "Residue Name", "Chain", "Residue ID"])
		
		# Make sure that the atoms / order match up between structures
		if(check_match and not (pdb1[["Atom Name", "Residue Name", "Residue ID"]].values == pdb2[["Atom Name", "Residue Name", "Residue ID"]].values).all()):
			print "ERROR: Atoms Name, Residue Name, Residue ID do not match"
			raise valueError
		
		# Generate coordinate arrays from each structure
		coor1 = pdb1[["X", "Y", "Z"]].values
		coor2 = pdb2[["X", "Y", "Z"]].values
	
	# For case where inputs are provided as two np arrays
	else:
		# Make sure the shapes mathc
		if(pdb1.shape != pdb2.shape):
			print "ERROR: Provided coordinates are not the same shape"
			raise valueError
		# Make sure the shapes look like coordinate arrays
		if(len(pdb1.shape) != 2 or pdb1.shape[1] != 3):
			print "ERROR: Provided coordinates should be shape (X, 3)"
			raise valueError
		
		# Rename the variables
		coor1 = pdb1
		coor2 = pdb2
	
	# Calculate Pairwise distance within each set
	dm1 = np.linalg.norm(coor1 - coor1[:,None], axis=-1)
	dm2 = np.linalg.norm(coor2 - coor2[:,None], axis=-1)
	
	# Ignore the Diagonal
	diagonal_mask = ~np.eye(dm1.shape[0], dtype=bool)
	
	# Define lower cuttoff (remove pairs of atoms that are too close)
	lower_cuttoff_mask = np.ones(dm1.shape) == 1
	if(lower_cuttoff != None):
		lower_cuttoff_mask = (dm1 >= lower_cuttoff) | (dm2 >= lower_cuttoff)
	
	# Define upper cuttoff (remove pairs of atoms that are too far away)
	upper_cuttoff_mask = np.ones(dm1.shape) == 1
	if(upper_cuttoff != None):
		upper_cuttoff_mask = (dm1 <= upper_cuttoff) | (dm2 <= upper_cuttoff)
	
	# Apply cuttoffs
	dm1 = dm1[diagonal_mask&upper_cuttoff_mask&lower_cuttoff_mask]
	dm2 = dm2[diagonal_mask&upper_cuttoff_mask&lower_cuttoff_mask]
	# Calcular dRMSD
	
	#if((upper_cuttoff_mask*lower_cuttoff_mask).all()):
	#	 return np.sqrt((1.0/(len(pdb1)*(len(pdb1) - 1))) * ((dm1 - dm2)**2).sum())
	#else:
	if((dm1.shape[0] - len(pdb1)) <= 0):
		return np.nan
	else:
		return np.sqrt((1.0/(dm1.shape[0] - len(pdb1))) * ((dm1 - dm2)**2).sum())
# FUNCTION END

if __name__ == '__main__':
	'''
	def test(x, y, *args, **kwargs):
		for a in args:
			print("Additional Argument: {0}".format(a))
		for key, value in kwargs.iteritems():
			print("Additional Named Argument: {0} = {1}".format(key, value))
		return x + y, args, kwargs
	
	def tmpTest(x):
		vprint(x)
		return x, reserveTemp()
	
	def meanTest(x):
		print x
		return x, np.mean(x)
	
	#multiThreadFunction(test, 3, 1, 8, [1, 2, 3], ["A", "B", "C"], x=[1, 5, 10], y=1, j=[1, 2, 3, 4], apple=5)
	#out = multiThreadFunction(tmpTest, 100, 40, False, x=range(100))
	#out = multiThreadFunction(meanTest, 1000, 40, x=[random.sample(xrange(100), 10) for i in range(1000)])
	
	def test1(x):
		return x
	
	def test2(*args):
		try:
			x, y = args
		except:
			x, y = args[0]
		return x + y
	
	multiThreadFunction(test1, 5, 1, True, x=range(5))
	
	multiThreadFunction(test2, 1, 1, True, (1, 5))
	'''
	
	print("Finished")

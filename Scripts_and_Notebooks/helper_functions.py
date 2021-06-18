from config import *

import re
from re import split as regexsplit

import pandas as pd
import numpy as np

import scipy.stats as st

from collections import defaultdict

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo

import os
import warnings

import requests
import urllib
import urllib2

from multiprocessing.pool import ThreadPool
from subprocess import Popen, PIPE

from gzip import open as gzopen
from bitarray import bitarray


# Obtain the protein sequence from the provided UniProt ID
# - Accepts as parameters...
#   @ID - UniProt ID
# - Returns the sequence for the given UniProt ID
# - Written for general utility purposes 03/23/17 by Shayne Wierbowski
def get_Fasta(ID):
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
# FUNCTION END

# Reads in a fasta format file as an (ID --> Seq) Dictionary
# - Accepts as parameters...
#   @filename - The fasta file to read
#   @only_len - Option to only store the length of each entry
#   @key_transform - Option to parse / edit the fasta header when entering key to dictionary
#   @keep - Option to only keep a specified sub-set of entries from the fasta
#   Returns the dictionary representation of the file
# - Written for general utility purposes 1/17/18 by Shayne Wierbowski
def fasta2dict(filename, only_len=False, key_transform=str, keep=None):
    # Parse keep list as a set
    if(keep != None):
        keep = set(keep)
    
    # Iterate through specified fasta file
    with open(filename, "r") as f:
        # Output dictionary
        res = dict()
        
        # Current key
        cur = None
        
        # Store either the total entry length or the sequence
        if(only_len):
            tmp = 0
        else:
            tmp = ""
        
        # Iterate through specified fasta file
        for l in f:
            # If this is the start of a new fasta entry
            if(">" in l):
                # If this is the first entry, just update the current key
                if(cur is None):
                    cur = key_transform(l.split(">")[1].strip())
                # Otherwise, save current entry to dictionary, and update current key
                else:
                    # Save current entry (unless this entry is filtered by the keep set)
                    if(keep is None or cur in keep):
                        res[cur] = tmp
                    
                    # Update the current key
                    cur = key_transform(l.split(">")[1].strip())
                
                # Update the current value (length or sequence)
                if(only_len):
                    tmp = 0
                else:
                    tmp = ""
            # Otherwise this line continues the current entry
            # just update the value
            else:
                if(only_len):
                    tmp += len(l.strip())
                else:
                    tmp += l.strip()
        
        # Save the final entry to the dictionary
        if(keep is None or cur in keep):
            res[cur] = tmp
    
    # Return
    return res
# FUNCTION END



# Map between 3 letter AA and 1 letter AA and vice versa
AA_3to1d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
              'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
              'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
              'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
AA_1to3d = {v: k for k, v in AA_3to1d.iteritems()}

# Converts 3 character AA symbol to 1 character AA symbol
def AA_3to1(seq, raise_errors=True):
     if len(seq) % 3 != 0: 
        if(raise_errors):
            raise ValueError('Input length should be a multiple of three')
        else:
            return "*"
     
     return "".join([AA_3to1d.get(seq[3*i:3*i+3], "*") for i in range(len(seq)/3)])
# FUNCTION END

# Converts 1 character AA symbol to 3 character AA symbol
def AA_1to3(seq):
    return "".join([AA_1to3d[x] for x in list(seq)])
# FUNCTION END

# Reads in a PDB File as a pandas DataFrame
# - Accepts as parameters...
#   @filename - PDB file to read (can also be plain pdb ID or raw PDB text)
#   @whitelist - Set of row types to parse retain
#   Returns the atom entries in the PDB as a pandas DataFrame
# - Written for general utility purposes 7/10/18 by Shayne Wierbowski
def pdb2df(filename, whitelist=["ATOM", "HETATM"]):
    # Allow direct pdb lines to be used instead of a filename
    if(type(filename) is list and type(filename[1]) is str):
        lines = filename
    
    # Try accessing file using mjm_tools open_pdb
    elif(not os.path.exists(filename)):
        # Obtain file handler
        fh = open_pdb(filename)
        
        # Read lines
        lines = [x.strip() for x in fh.readlines()]
    else:
        # Handle unzipping gzipped files
        if(filename.split(".")[-1] == "gz"):
            tmp = reserveTemp()
            call("gunzip -c {0} > {1}".format(filename, tmp))
            filename = tmp
        
        # Read PDB File
        lines = [x.strip() for x in open(filename, "r").readlines()]
     
    # Filter lines to only interesting records saving header and tailer
    tmp = [i for i in range(len(lines)) if lines[i][:6].strip() in whitelist]
    header = lines[:min(tmp)]
    tailer = lines[max(tmp) + 1:]
    
    # Identify breaking point between different models (for our purposes we
    # ONLY look at the first model if there are every multile models)
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
    def apply_type(x, t):
        if(x == ""):
            return np.nan
        else:
            return t(x)
    lines = [[apply_type(x[points[i]:points[i+1]].strip(), types[i]) for i in range(len(points) - 1)] for x in lines]
    
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

# Writes a DF representation of a PDB to a PDB file
# - Accepts as parameters...
#   @filename - PDB file to write to
#   @df - DataFrame containing PDB atom rows
#   @columns - names of columns to use
#   Returns the atom entries in the PDB as a pandas DataFrame
# - Written for general utility purposes 1/15/19 by Shayne Wierbowski
def df2pdb(filename, df, columns=["Data Type", "Atom ID", "Atom Name", "Alternate", "Residue Name", "Chain", "Residue ID", "Insertions?", "X", "Y", "Z", "Occupancy", "Temperature Factor", "Segment ID", "Element"]):
    out = open(filename, "w+")
    indices = [list(df).index(c) for c in columns]
    s = ""
    for i in range(len(df)):
        atom, atmi, atmn, alt, resn, chain, resi, insertions, x, y, z, occupancy, temp, segid, esym, = df.iloc[i, indices].replace(np.nan, "")
        s += "{atom:<6}{atmi:>5} {atmn:<4}{alt:<1}{resn:>3} {chain}{resi:>4.4}{insertions:>1}    {x:>7.7} {y:>7.7} {z:>7.7}{occupancy:>6}{temp:>6}      {segid:<4}{esym:>2}\n".format(atom=atom, atmi=atmi, atmn=atmn, alt=alt, resn=resn, chain=chain, resi=str(resi), insertions=insertions, x=str(x)[:7], y=str(y)[:7], z=str(z)[:7], occupancy=occupancy, temp=temp, segid=segid, esym=esym)
    out.write(s)
    out.close()
# FUNCTION END



# Calculates an odds ratio for a provided set of case / exposures (formatted as boolean array)
# - Accepts as parameters...
#   @exposure_mask - array (usually numpy boolean array) of true / false conditions for the "exposed" status in the odds ratio
#   @case_mask - array (usually numpy boolean array) of true / false conditions for the "case" status in the odds ratio
#   @CI_value - p-value for significance defaults to 0.05 (i.e. return 95% confidence intervals)
#   @two-sided - option to have the significance test calculated as two-sided instead of one-sided (defaults to two-sided)
#   @log_odds - option to have the output automatically converted to log2 space (defaults to leave in non-log space)
#   @verbose - option to display the full case-exposure contigency matrix when run
#   @long_output - option to add the raw a, b, c, d counts from the contigency matrix to the output
#   @expose_label - label to use for exposed vs. non-exposed status in verbose output
#   @case_label - label to use for case vs. non-case status in verbose output
#   @error - option to return an upper / lower confidence interval vs. +/- standard error (defaults to CI)
#   Returns the odds ratio
# - Written for general utility purposes by Shayne Wierbowski
def odds_ratio(exposure_mask, case_mask, CI_value=0.05, two_sided=True, log_odds=False, verbose=False, long_output=False, expose_label="Exposed", case_label="Case", error="CI"):
    # Catch all for various types formats I've provided the masks in before (first case is for direct binary arrays which are
    # probably the most efficient, but I typically use numpy arrays instead)
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
    
    # Deside what the reported error bars should be (either confidence interval of standard error)
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

# - Flattens a nested list (e.g. [("A", "B"), ("C", "D")] --> ["A", "B", "C", "D"])
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




################################
#                              #
# Alignemnt Function + Helpers #
#                              #
################################

# Alignment Helper Functions
def alignPident(align1, align2, useShorter=True):
    align1 = align1.upper()
    align2 = align2.upper()
    
    # Convert alignments to array, sum up identical AA
    identical = sum(np.array((list(align1)) == np.array(list(align2)))&(np.array(list(align1)) != "-"))
    
    # Here the percent identity is calculated relative to the total length of the shorter
    # sequence by default. However, there are multiple approaches that could also be considered for the denominator...
    # - Length of the longer sequence
    # - Length of the full alignemnt
    # - Length of the aligned portion (exclude gaps)
    # - etc.
    if useShorter:
        totalLength = min([sum(np.array(list(align1)) != "-"), sum(np.array(list(align2)) != "-")])
    else:
        totalLength = max([sum(np.array(list(align1)) != "-"), sum(np.array(list(align2)) != "-")])
    
    # Return Percentage
    return float(identical)/totalLength
# FUNCTION END

# Alignment Helper Functions
def getBlosumScore(x, y, matrix=MatrixInfo.blosum62):
    if((x, y) in matrix):
        return matrix[(x, y)]
    else:
        return matrix[(y, x)]
# FUNCTION END

# Alignment Helper Functions
def alignPositives(align1, align2, matrix=MatrixInfo.blosum62, useShorter=True):
    # Calculate Blosum score over alilgnment
    blosumScores = [getBlosumScore(align1[i], align2[i], matrix=matrix) for i in range(len(align1)) if not (align1[i] == "-" or align2[i] == "-")]
    
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

# Alignment Helper Functions
def alignCoverage(align1, align2, matrix=MatrixInfo.blosum62, useShorter=True):
    # Sum up non-gapped regions of alignment
    non_gaps = len(align1) - sum([align1[i] == "-" or align2[i] == "-" for i in range(len(align1))])
    
    # Divide by sequence length of align1
    return float(non_gaps)/len(align1.replace("-", ""))
# FUNCTION END

def getDNASubMatrix():
    df = pd.read_csv("{0}/EDNAFULL.sub".format(resource_dir), sep="\t", index_col=0)
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

# Generate Positon Map from alignment info
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
    
    pos_map = pd.DataFrame(pos_map, columns=[name1 + "_Pos", name1 + "_AA", name2 + "_Pos", name2 + "_AA"])
    
    return pos_map
# FUNCTION END

# Pretty Prints an alignemnt
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

# Alignment Function
# - Accepts as parameters...
#   @s1 - First sequence to align
#   @s2 - Second seqence to align
#   @matrix - Substitution matrix to use for scoring (will default to standard BLOSUM 62 or DNA substitution matrix) if not provided
#   @useShorter - Calculate coverage with respect to the shorter sequence
#   @best - Return the first alignment option (set to false and instead re-ranks by overall Pident because from experience the "top"
#           alignment is not always the best
#   @show_align - Option to print the alignment
#   @genPosMap - Option to return a position map from s1_position --> s2_position
#   Returns the alignment between s1 and s2
# - Written for general utility purposes 7/17/18 by Shayne Wierbowski
def NWSeqAlignment(s1, s2, matrix=None, useShorter=True, best=False, show_align=False, genPosMap=False):
    keep = None
    
    # Fill in the appropriate default substitution matrix if none provided
    kind = "prot"
    if(matrix == None):
        # First assume this is a protein alignment and use BLOSUM62
        matrix = MatrixInfo.blosum62
        
        # Check if the sequences look like DNA
        if(len(set(list(s1.upper()) + list(s2.upper())).difference(set(["A", "C", "G", "T", "N"]))) == 0):
            kind = "nucl"
            # Pull out default DNA substitution matrix
            matrix = getDNASubMatrix()
    
    # Alignment parameters for protein alignment
    if(kind == "prot"):
        align = pairwise2.align.globalds(s1, s2, matrix, -10, -0.5)
    # Alignment parameters for DNA alignment
    else:
        align = pairwise2.align.globalds(s1, s2, matrix, -400, 0)
    
    # Reorder alignment outputs to break ties and fix poor gap placement.
    # e.g. Default orders might come up with something like this...
    #
    #		AAAAAG--------G
    #		AAAAAGGCCCCCCCG
    #
    # instead of this...
    #
    #		AAAAAGG--------
    #		AAAAAGGCCCCCCCG
    align = sorted(align, key=lambda x: (-x[2], len(re.findall("-[^-]", x[0])) + len(re.findall("-[^-]", x[1]))))
    
    # Iterate over all alignments and parse / select one to keep
    for a in align:
        a1, a2, score, begin, end = a
        
        # Format alignment as a description + add in other statistics (Pident, Positives, Coverage)
        r = dict([("Align1", a1), ("Align2", a2), ("Score", score), ("Begin", begin), ("End", end), ("Pident", alignPident(a1, a2, useShorter=useShorter)), ("Positives", alignPositives(a1, a2, matrix=matrix, useShorter=useShorter)), ("Coverage1", alignCoverage(a1, a2, matrix=matrix, useShorter=useShorter)), ("Coverage2", alignCoverage(a2, a1, matrix=matrix, useShorter=useShorter))])
        
        # Construct string describing alignment (e.g. identical "|", positive "+", negative "-", or gap " ")
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
        
        # If this is the first choice we've parsed or has the best pident so far keep this alignment
        if(keep == None or r["Pident"] >= keep["Pident"]):
            keep = r
        
        # Or if we just want the "best" alignment break out of the loop right away
        if(not best):
            break
    
    # If we want to display the alignment, display it now
    if(show_align):
        alignPrint(keep)
    
    # If we want to generate a position map, generate it and return alongside the alignment
    if(genPosMap):
        return keep, alignPosMap(keep)
    # Otherwise just return the alignment
    else:
        return keep
# FUNCTION END



############################
#                          #
# Extracted from mjm_tools #
#                          #
############################

# NOTE: The following functions are extracted from the mjm_tools.py library on the Yu Lab server
#       all code here is originally written / provided by Michael Meyer


def unzip_res_range(res_range):
    '''Converts ranges in the form: [2-210] or [3-45,47A,47B,51-67] into lists of strings including all numbers in these ranges in order'''
    res_ranges = res_range.strip()[1:-1].split(',')
    index_list = []
    for r in res_ranges:
        if re.match('.+-.+', r):
            a, b = r.split('-')
            index_list += [str(n) for n in range(int(a), int(b)+1)]
        else:
            index_list.append(r)
    
    if index_list == ['']:
        return []
    else:
        return index_list
# FUNCTION END

def zip_res_range(seq):
    '''Converts lists of residues as string in the form "1,2,2E,3,4,5,6,46,67,68,68A,68C,69,70" to
    zipped ranges in the form "[1-2,2E,3-6,46,67-68,68A,68C,69-70]" @author: haoran'''
    
    if type(seq) == list:
        seq = ','.join(seq)
    
    seqout = []
    tempst = ''
    for h in range(seq.count(',')+1):
        i = seq.split(',')[h]
        if i.isdigit() and h < seq.count(',') and str(int(i)+1) == seq.split(',')[h+1]:
            if tempst == '':
                tempst = i
        else:
            if tempst == '':
                seqout.append(i)
            else:
                seqout.append(tempst+'-'+i)
                tempst = ''
    return '['+','.join(seqout)+']'
# FUNCTION END


def natural_keys(text):
    '''natural sorting (adapted from http://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside)
          Useful for sorting PDB residues that have letters in them
       
       example:
       
       >> from mjm_tools import natural_keys
    
       >> my_list = ['1', '3', '2', '7', '2B']   
       >> my_list.sort(key=natural_keys)
       ['1', '2', '2B', '3', '7']
          
          '''
    def atoi(text): return int(text) if text.isdigit() else text
    return [ atoi(c) for c in regexsplit('(\d+)', text) ]
# FUNCTION END


def is_binary_file(filename):
    ''' check if file is binary (adapted from http://stackoverflow.com/questions/898669/how-can-i-detect-if-a-file-is-binary-non-text-in-python) '''
    textchars = bytearray([7,8,9,10,12,13,27]) + bytearray(range(0x20, 0x100))
    return bool(open(filename, 'rb').read(1024).translate(None, textchars))
# FUNCTION END


# NOTE: The way this function is used assumes that there is a local copy of the PDB
#       on the user's machine and path to this copy provided in the configure.py. If this
#       is not the case, requests will be made through the web and may be subject to considerable
#       slowdown and / or errors from too many requests.
def open_pdb(structure, verbose=True, try_web=True):
    '''Return an opened PDB file handle from STDIN, file, local PDB cache, or web'''

    # STDIN
    if "<open file '<stdin>', mode 'r' at" in str(structure):
        pdb_filehandle = structure

    # AS UNCOMPRESSED PDB FILE
    elif os.path.exists(structure) and is_binary_file(structure)==False:   #file exists and is a text-based file
        pdb_filehandle = open(structure, 'r')
        
    # AS GZIPPED PDB FILE
    elif os.path.exists(structure) and is_binary_file(structure)==True:    #file exists and is likely a gzipped file
        try:
            testopen = gzopen(structure, 'r')
            testopen.readline()
            testopen.close()
            pdb_filehandle = gzopen(structure, 'r')
        except IOError:
            if(verbose):
                print 'Invalid structure file-type. Structure file must be a plain-text PDB file or a gzipped PDB file.'
            return

    # AS PDB FILE FROM LOCAL COPY OF THE PDB -OR- FROM THE WEB
    elif len(structure)==4:
        
        pdb_storage_path = os.path.join(PDB_DATA_DIR, '%s/pdb%s.ent.gz' %(structure[1:3].lower(), structure.lower()))
        
        #local file
        if os.path.exists(pdb_storage_path):
            pdb_filehandle = gzopen(pdb_storage_path, 'r')
        #try the web
        elif(try_web):
            try:
                pdb_filehandle = urlopen('http://www.rcsb.org/pdb/files/%s.pdb' %(structure.upper()))
            except HTTPError:
                if(verbose):
                    print 'Invalid structure input: %s. Not found as local file, as PDB structure in %s, or on the web.' %(structure, PDB_DATA_DIR)
                return
        else:
            return
    else:
        if(verbose):
            print 'Invalid structure input: %s. Not found as local file, and wrong number of characters for direct PDB reference.' %(structure)
        return
    
    return pdb_filehandle
# FUNCTION END



def parse_naccess(naccess_output):
    '''	Parse naccess output into dictionary.
    
        Example head of naccess output file (.rsa):
    
        REM  Relative accessibilites read from external file "/usr/local/naccess2.1.1/standard.data"
        REM  File of summed (Sum) and % (per.) accessibilities for 
        REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
        REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
        RES SER A  14    28.88  24.8  22.22  28.5   6.66  17.3  18.69  38.5  10.19  15.0
        RES ALA A  15    10.35   9.6   0.00   0.0  10.35  26.8   2.33   3.3   8.01  21.9
        RES ASN A  16    37.02  25.7  36.33  34.2   0.69   1.8  25.49  55.1  11.53  11.8
        RES LEU A  17   130.51  73.1 105.12  74.5  25.39  67.7 105.51  74.1  25.00  68.8
        RES ASP A  18   103.34  73.6 100.05  97.4   3.29   8.7  33.55  68.1  69.79  76.6
        RES HIS A  19   143.10  78.2 123.76  84.1  19.34  54.0  75.48  77.7  67.62  78.9
        
    '''
    
    parsed_output = []
    
    for l in naccess_output:
        if l[:3] != 'RES':
            continue
        
        res, chain, res_num = l[3:7].strip(), l[7:9].strip(), l[9:14].strip()
        combined_key = l[4:14]
        
        all_atoms_abs, all_atoms_rel = float(l[14:22].strip()), float(l[22:28].strip())
        total_side_abs, total_side_rel = float(l[28:35].strip()), float(l[35:41].strip())
        main_chain_abs, main_chain_rel = float(l[41:48].strip()), float(l[48:54].strip())
        non_polar_abs, non_polar_rel = float(l[54:61].strip()), float(l[61:67].strip())
        all_polar_abs, all_polar_rel = float(l[67:74].strip()), float(l[74:80].strip())
        
        parsed_output.append({	'res': res, 'chain': chain, 'res_num': res_num, 'combined_key': combined_key,
                            'all_atoms_abs': all_atoms_abs, 'all_atoms_rel': all_atoms_rel,
                            'total_side_abs': total_side_abs, 'total_side_rel': total_side_rel,
                            'main_chain_abs': main_chain_abs, 'main_chain_rel': main_chain_rel,
                            'non_polar_abs': non_polar_abs, 'non_polar_rel': non_polar_rel,
                            'all_polar_abs': all_polar_abs, 'all_polar_rel': all_polar_rel })
        
    return parsed_output
# FUNCTION END


def naccess(pdb_file, parsed=False):
    '''Run NACCESS and return the results.'''
    cwd = os.getcwd()     #naccess writes to the current working directory, so save current directory, and move to scatchDir to write output files
    os.chdir(os.path.dirname(pdb_file))   #write naccess output files to directory of PDB file

    raw_naccess_output = []
    
    # Delete output files if they already exist (EDIT made by Shayne July 12, 2018)
    # Leaving these files can be a problem when running naccess multiple times (i.e. irescalc.py)
    # Since if the naccess command fails, this function will see the previous outputs and not detect the failure
    if(os.path.exists(os.path.splitext(pdb_file)[0]+'.rsa')):
        os.remove(os.path.splitext(pdb_file)[0]+'.rsa')
    if(os.path.exists(os.path.splitext(pdb_file)[0]+'.asa')):
        os.remove(os.path.splitext(pdb_file)[0]+'.asa')
    
    _, _ = Popen([NACCESS_PATH, pdb_file], stdout=PIPE, stderr=PIPE).communicate()
    try:
        raw_naccess_output += open(os.path.splitext(pdb_file)[0]+'.rsa').readlines()
        
        # Make sure output files are not empty (EDIT made by Shayne July 12, 2018)
        # See above
        if(os.stat(os.path.splitext(pdb_file)[0]+'.rsa').st_size == 0):
            print 'ERROR: Naccess .rsa file is empty. We suspect this is an edge case where Naccess cannot calculate ASA for extremely large chains. The following command was attempted: %s %s' %(NACCESS_PATH, pdb_file)
            print 'WARNING: This error message is only a temporary measure. The longterm solution will involve identifying an alternative SASA calculator. If you see this message, please contact Shayne.'
            exit()
    except IOError:
        print 'ERROR: Naccess .rsa file was not written. The following command was attempted: %s %s' %(NACCESS_PATH, pdb_file)
        exit()

    os.chdir(cwd)  #move back to old cwd
    
    if parsed:
        return parse_naccess(raw_naccess_output)
    else:
        return raw_naccess_output
# FUNCTION END



def fetch_uniprot(uniprots, attributes):
    '''Fetch information about a UniProt ID from uniprot.org
    UniProt attributes can be a single string, which will return a single string value,
    or a list of attributes, which will return a dictionary of results attribute => value
    
    accepted attributes:
    
    citation | clusters | comments | domains | domain | ec | id
    entry name | existence | families | features | genes | go
    go-id | interactor | keywords | last-modified | length
    organism | organism-id | pathway | protein names | reviewed
    sequence | 3d | version | virus hosts
    
    For more info, see: http://www.uniprot.org/help/programmatic_access
    '''
    
    if type(uniprots) == str:
        uniprots = [uniprots,]
    
    if type(attributes) == str:
        attributes = [attributes,]
    
    if 'id' not in attributes:
        attributes.append('id')
    
    # EDIT: Shayne Wierbowski July 24, 2018
    # UniProt API seems to have changed, joining uniprot IDs by commas no longer seems to work
    # I have modified this code to join uniprot IDs by "+OR+" instead as if you were building a query
    # directly as described at https://www.uniprot.org/help/api_queries
    # I am not sure why a similar change was not required to the columns field
    # Changed from...
    #uniprots = ','.join(uniprots)
    # To...
    uniprots = '+OR+'.join(uniprots)
    columns = ','.join(attributes)
    
    #build and send API request
    # EDIT: Shayne Wierbowski July 24, 2018
    # UniProt API now seems to require https rather than http url
    #Changed from...
    #url = 'http://www.uniprot.org/uniprot/' by Shayne Wierbowski July 24, 2018 because
    # To...
    url = 'https://www.uniprot.org/uniprot/'
    params = {'query': uniprots,
              'format': 'tab',
              'columns': columns}
    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    request.add_header('User-Agent', 'Python mjm659@cornell.edu')
    response = urllib2.urlopen(request)
    data = response.readlines()
    ###

    if data==[] or len(data)<2:
        return {}
    
    # EDITS MADE BY SHAYNE WIERBOWSKI ON OCTOBER 30m 2019
    # Added extra_dict / warnings for extra IDs that get fetched but are not
    # requested. Also added handling for "Merged Into" rows (e.g.)
    #
    # P18463	1B37_HUMAN							Merged into P01889.	
    #
    return_dict = {}
    extra_dict = {}
    redo = []
    
    for line in data[1:]:
        
        line_data = line.strip().split('\t')
        if len(line_data) != len(attributes):
            if("Merged into" in line):
                redo.append((line.strip().split("\t")[0], line.strip().split("\t")[-1].split()[-1].replace(".", "")))
            continue
        
        entry_dict = {}
        for i, a in enumerate(attributes):
            entry_dict[a] = line_data[i]
        if entry_dict['id'] in uniprots:
            return_dict[entry_dict['id']] = entry_dict
        else:
            extra_dict[entry_dict['id']] = entry_dict
    
    if(len(extra_dict) != 0):
        # Apparrently this just happens all the time, but I can't figure out why?
        #warnings.warn("Unrequested UniProts ({0}) were returned. This should not happen.".format(";".join(extra_dict.keys())))
        pass
    
    for p1, p2 in redo:
        if(p2 in return_dict):
            return_dict[p1] = copy.deepcopy(return_dict[p2])
            return_dict[p1]["id"] = p1
        elif(p2 in extra_dict):
            return_dict[p1] = copy.deepcopy(extra_dict[p2])
            return_dict[p1]["id"] = p1
        warnings.warn("Requested UniProt ID, {0}, appears to have been replaced by {1}".format(p1, p2))
    
    # END EDITS
    
    return return_dict
# FUNCTION END


############################
#                          #
# Extracted from jfb_tools #
#                          #
############################

# NOTE: The following functions are extracted from the jfb_tools.py library on the Yu Lab server
#       all code here is originally written / provided by Juan Felipe Beltran

# Also see the UniProt API documentation
# https://www.uniprot.org/help/api_idmapping
def callUniprotAPI(package):
    
    source_id, target_id, contents = package
    url = 'https://www.uniprot.org/mapping/'
    
    params = {
        'from': source_id,
        'to': target_id,
        'format': 'tab',
        'query': contents
    }
    
    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    contact = uniprot_contact_email
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    page = response.read(200000)
    return [line.split('\t') for line in page.strip().split('\n')[1:]]
# FUNCTION END

def batchUniProtAPI(to_translate, source_id, target_id, n_threads=50, batch_size=20):
    
    pool = ThreadPool(n_threads)
    mapper = defaultdict(str)
    
    batched_contents = [(source_id, target_id, ' '.join(to_translate[i:i + batch_size])) for i in range(0, len(to_translate), batch_size)]
    responses = pool.map(callUniprotAPI, batched_contents)
    
    for key, val in [(x, y) for group in responses for (x, y) in group]:
        mapper[key] = val
    return [mapper[val] for val in to_translate]
# FUNCTION END















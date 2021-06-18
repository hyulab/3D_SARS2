from config import *

import sys
import os
import time

import subprocess as sp
import math

start = time.time()

# NOTE: I had to mess around with the PATH variables for some reason. This script worked in the command line, but didn't work in when
#       called through jupyter unless I manually tweaked the paths. Likely an issue related to installation (e.g. jupyter server was
#       already running when HADDOCK installed so paths from subprocesses spawned by my jupyter shell were not updated correctly?)
#
#       Hopefully this will not be necessary for anyone trying to replicate this...
#os.environ["PYTHONPATH"] = os.environ["PYTHONPATH"] + ":" + haddock_path

# Logger
def pprint(s):
	print s
	log.write(str(s) + "\n")
# FUNCTION END

# Constants
starting_dir = raw_docking_dir
os.chdir(starting_dir)

chain1 = "A"
chain2 = "B"

pymol_sele_path = "{0}/execute_pymol_sele.py".format(script_dir)
sres_path = "{0}/srescalc.py".format(script_dir)

sres_low_thresh = 15
sres_high_thresh = 40

passive_thresh = 6

haddock_cmd = "{0} {1}/Haddock/RunHaddock.py".format(python_path, haddock_path)



# Parse Inputs
name = sys.argv[1]
os.system("mkdir -p {0}".format(name))
os.chdir(name)
log = open("wrapper.log", "w+")

pdb1 = sys.argv[2]
pdb2 = sys.argv[3]

ires1 = sys.argv[4]
ires1 = [int(x) for x in ires1.split(",")]
ires2 = sys.argv[5]
ires2 = [int(x) for x in ires2.split(",")]
pprint("Ires1 (Init): " + ",".join([str(x) for x in ires1]))
pprint("Ires2 (Init): " + ",".join([str(x) for x in ires2]))
pprint("")

passive1 = sys.argv[6]
if(passive1 == "None"):
	passive1 = []
else:
	passive1 = [int(x) for x in passive1.split(",")]
passive2 = sys.argv[7]
if(passive2 == "None"):
	passive2 = []
else:
	passive2 = [int(x) for x in passive2.split(",")]
pprint("Passive1 (ECLAIR Medium): " + ",".join([str(x) for x in passive1]))
pprint("Passive2 (ECLAIR Medium): " + ",".join([str(x) for x in passive2]))
pprint("")



# Calculate surface residues for both chains
pprint("Calculating Surface Residues")
sres1_low = set([int(x) for x in sp.check_output("{0} -c {1} {2} -uSASA {3}".format(sres_path, chain1, pdb1, sres_low_thresh), shell=True).split("\n")[0].split(",")])
sres1_high = set([int(x) for x in sp.check_output("{0} -c {1} {2} -uSASA {3}".format(sres_path, chain1, pdb1, sres_high_thresh), shell=True).split("\n")[0].split(",")])
sres2_low = set([int(x) for x in sp.check_output("{0} -c {1} {2} -uSASA {3}".format(sres_path, chain2, pdb2, sres_low_thresh), shell=True).split("\n")[0].split(",")])
sres2_high = set([int(x) for x in sp.check_output("{0} -c {1} {2} -uSASA {3}".format(sres_path, chain2, pdb2, sres_high_thresh), shell=True).split("\n")[0].split(",")])
pprint("Sres1_low: " + ",".join([str(x) for x in sres1_low]))
pprint("Sres1_high: " + ",".join([str(x) for x in sres2_high]))
pprint("Sres2_low: " + ",".join([str(x) for x in sres1_low]))
pprint("Sres2_high: " + ",".join([str(x) for x in sres2_high]))
pprint("")

# Filter active residues to only the surface ones
ires1 = [x for x in ires1 if x in sres1_low]
ires2 = [x for x in ires2 if x in sres2_low]
pprint("Ires1 (Filtered): " + ",".join([str(x) for x in ires1]))
pprint("Ires2 (Filtered): " + ",".join([str(x) for x in ires2]))
pprint("")

passive1 = [x for x in passive1 if x in sres1_low]
passive2 = [x for x in passive2 if x in sres2_low]
pprint("Passive1 (Filtered): " + ",".join([str(x) for x in passive1]))
pprint("Passive2 (Filtered): " + ",".join([str(x) for x in passive2]))
pprint("")

# If no high or very high confidence ires, use the medium (passive) as actives
if(len(ires1) == 0):
	ires1 = passive1
	passive1 = []
if(len(ires2) == 0):
	ires2 = passive2
	passive2 = []

# Identify passive residues
pprint("Identifying Passive Residues")
if(len(ires1) != 0):
    sele_string = "byresi obj within {0} of (resi {1})".format(passive_thresh, " or resi ".join([str(x) for x in ires1]))
    passive1 = set(passive1 + [int(x) for x in sp.check_output("python {0} {1} \"{2}\"".format(pymol_sele_path, pdb1, sele_string), shell=True).split("\n")[0].split(",")])
    pprint("Pres1 (Init): " + ",".join([str(x) for x in passive1]))
    passive1 = [x for x in passive1 if x in sres1_low and not x in ires1]
    pprint("Pres1 (Filtered): " + ",".join([str(x) for x in passive1]))
else:
    passive1 = sres1_high

if(len(ires2) != 0):
    sele_string = "byresi obj within {0} of (resi {1})".format(passive_thresh, " or resi ".join([str(x) for x in ires2]))
    passive2 = set(passive2 + [int(x) for x in sp.check_output("python {0} {1} \"{2}\"".format(pymol_sele_path, pdb2, sele_string), shell=True).split("\n")[0].split(",")])
    pprint("Pres2 (Filtered): " + ",".join([str(x) for x in passive2]))
    passive2 = [x for x in passive2 if x in sres2_low and not x in ires2]
    pprint("Pres2 (Filtered): " + ",".join([str(x) for x in passive2]))
else:
    passive2 = sres2_high
pprint("")

# Write active / passive lists and create AIR restraints file
pprint("Generating restraint files")
pprint("")
act_pass1 = os.path.basename(pdb1).replace(".pdb", "act_pass.txt")
out = open(act_pass1, "w+")
out.write(" ".join([str(x) for x in ires1]) + "\n")
out.write(" ".join([str(x) for x in passive1]))
out.close()

act_pass2 = os.path.basename(pdb2).replace(".pdb", "act_pass.txt")
out = open(act_pass2, "w+")
out.write(" ".join([str(x) for x in ires2]) + "\n")
out.write(" ".join([str(x) for x in passive2]))
out.close()

#os.system("active-passive-to-ambig.py {0} {1} > {2}_AIR.txt".format(act_pass1, act_pass2, name))
# Only use CA for the restraint
os.system("active-passive-to-ambig.py {0} {1} | sed s/segid/name\ CA\ and\ segid/g | sed s/2.0/3.0/g > {2}_AIR.txt".format(act_pass1, act_pass2, name))


# Write run.params file
s = '''AMBIG_TBL=./{0}_AIR.txt
HADDOCK_DIR={1}
N_COMP=2
PDB_FILE1={2}
PDB_FILE2={3}
PROJECT_DIR=.
PROT_SEGID_1=A
PROT_SEGID_2=B
RUN_NUMBER=1
'''.format(name, haddock_path, pdb1, pdb2)
out = open("run.param", "w+")
out.write(s)
out.close()


# Run Haddock (Sets up the run)
pprint("Initializing Haddock")
pprint("")
pprint(os.system("{0} >&init.log".format(haddock_cmd)))

# Move into run directory
os.chdir("run1")

# Run Haddock again (does the docking)
pprint("Running Haddock Docking")
pprint("")
pprint(os.system("{0} >&haddock.out".format(haddock_cmd)))

end = time.time()
time = end - start
pprint("Finished" + "{0:02.0f}:{1:02.0f}:{2:02g}".format(int(math.floor(float(time) / 60 / 60)), int(math.floor(float(time) / 60 % 60)), float(time)%60))

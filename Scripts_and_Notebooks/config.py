import os

#############################
# User Specified Parameters #
#############################

# Base directory for the whole project (i.e. the path to where you've
# cloned this repository)
base_dir = ""
# If not specified try to guess based on where this script is being
# run from
if(base_dir == ""):
    base_dir = "/".join(os.getcwd().split("/")[:-1])


# Contact email for UniProt Batch Requests
uniprot_contact_email = ""

# Local Copy of the PDB
# I'm unsure if the pipeline runs without error if this is left
# blank. Should Theoretically be uneccessary.
PDB_DATA_DIR = "/"

# HADDOCK Path
# Should be something like...
# /home/sdw95/bin/haddock2.4-2020-09
haddock_path = ""

# Naccess Path
# after install naccess, find this by adding to your PATH and running "which naccess"
NACCESS_PATH = ""

# Python Path
# path to whicever version of python you want to run with
python_path = ""



##################
# File Locations #
##################

# These parameters should all be constants that can be derived as long as a base directory
# is specified. But if you wanted to change the names of any directories you could tweak
# these values. Not recommended.


# Location of Scripts / Resources
script_dir = "{0}/Scripts_and_Notebooks".format(base_dir)
resource_dir = "{0}/static_resources".format(base_dir)

# Location of all Input files (for the purposes of this Notebook,
# this is fixed to a small set of representative examples)
input_dir = "{0}/Inputs".format(base_dir)

# Location of all output / data produced from this pipeline
output_dir = "{0}/Output".format(base_dir)
if(not os.path.exists(output_dir)):
	os.mkdir(output_dir)


# Location to save raw / parsed ECLIAR predictions
# NOTE: Code to run ECLIAR pipeline not provided, example outputs
#       provided as part of the static resources
raw_eclair_dir = "{0}/Precomputed_ECLAIR_Predictions".format(resource_dir)
eclair_dir = "{0}/Eclair_Predictions".format(output_dir)
if(not os.path.exists(eclair_dir)):
	os.mkdir(eclair_dir)


# Location of raw / parsed HADDOCK PPI Docking Runs
raw_docking_dir = "{0}/Docking_Runs".format(base_dir)
docking_dir = "{0}/Docked_Structures".format(output_dir)
if(not os.path.exists(docking_dir)):
	os.mkdir(docking_dir)


# Location of all ddG SARS1 --> SARS2 Mutated Sturctures
ddg_dir = "{0}/ddG_Mutated_Structures/".format(output_dir)
if(not os.path.exists("{0}".format(ddg_dir))):
    os.mkdir("{0}".format(ddg_dir))
if(not os.path.exists("{0}Structures/".format(ddg_dir))):
    os.mkdir("{0}Structures/".format(ddg_dir))
if(not os.path.exists("{0}Summary_Logs/".format(ddg_dir))):
    os.mkdir("{0}Summary_Logs/".format(ddg_dir))

# Location of all ddG interface scanning single mutant structures
ddg_singles_dir = "{0}/ddG_Single_Mutants/".format(output_dir)
if(not os.path.exists("{0}".format(ddg_singles_dir))):
    os.mkdir("{0}".format(ddg_singles_dir))
if(not os.path.exists("{0}Raw_Outputs".format(ddg_singles_dir))):
    os.mkdir("{0}Raw_Outputs".format(ddg_singles_dir))
if(not os.path.exists("{0}Summaries".format(ddg_singles_dir))):
    os.mkdir("{0}Summaries".format(ddg_singles_dir))
if(not os.path.exists("{0}Structures".format(ddg_singles_dir))):
    os.mkdir("{0}Structures".format(ddg_singles_dir))
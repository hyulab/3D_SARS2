from collections import defaultdict
import pandas as pd
import os

ppi_docking = pd.read_csv("Data/Docking_Summary.txt", sep="\t")
ppi_docking = ppi_docking.groupby(["P1", "P2"]).head(10)

viral2muts = defaultdict(set)

viral_muts = pd.read_csv("Data/Viral_Muts.txt", sep="\t")
for prot, pos, aa_ref, aa_alt in viral_muts[["COVID_ID", "COVID_Pos", "COVID_AA", "SARS_AA"]].values:
    if(not (aa_alt == "" or pd.isnull(aa_alt))):
        viral2muts[prot].add((pos, aa_alt))

from tqdm import tqdm, tqdm_notebook
from tqdm._tqdm_notebook import tqdm_notebook

tqdm.pandas(tqdm_notebook)
tqdm_notebook.pandas()

import helper as my
import time

def run_ddg_calc(p1, p2, pdb, rank, muts, trials, interface_cutoff):
    cmd = "nice python Mutant_ddG_Calc.py {0} {1} {2} {3} {4} {5} {6}".format(p1, p2, pdb, rank, muts, trials, interface_cutoff)
    p = sp.Popen(cmd, shell=True, stdout=FNULL, stderr=sp.STDOUT)
    return {"p":p, "cmd":cmd, "start_time":time.time(), "end_time":None, "in":(p1, p2, rank)}
# FUNCTION END

import subprocess as sp
# Run ddG in Loop
i_num = 1
finished_processes = []
processes = []

max_processes = 50

trials = 50
interface_cutoff = 8.0
FNULL = open(os.devnull, 'w')
# Iterate over all Docking Inputs
for p1, p2, pdb, rank in tqdm(ppi_docking.sort_values(["Rank", "P1", "P2"])[["P1", "P2", "File", "Rank"]].values):
    if(os.path.exists("Data/ddG_Mutated_Structures/Summary_Logs/{0}_{1}_{2}".format(p1, p2, rank))):
        continue
    
    muts = viral2muts[p1]
    
    pdb_df = my.pdb2df(pdb)
    resis = set(pdb_df[pdb_df["Chain"] == "A"]["Residue ID"].unique())
    
    muts = [x for x in muts if x[0] in resis]
    
    if(len(muts) == 0):
        continue
    
    muts = ",".join("_".join([str(y) for y in x]) for x in muts)
    
    # Make Block new jobs if too many running already
    while(True):
        if(len(processes) <= max_processes):
            p = run_ddg_calc(p1, p2, pdb, rank, muts, trials, interface_cutoff)
            processes.append(p)
            break
        else:
            new_processes = []
            for p in processes:
                if(p["p"].poll() is None):
                    new_processes.append(p)
                elif(p["p"].poll() != 0):
                    p["end_time"] = time.time()
                    print "Error", p["p"].poll()
                    print "cmd:", p["cmd"]
                    print "RunTime:", (p["end_time"] - p["start_time"])
                    print
                    finished_processes.append(p)
                else:
                    p["end_time"] = time.time()
                    print "Finished ddG", p["in"], "in", (p["end_time"] - p["start_time"])
                    finished_processes.append(p)
            processes = new_processes
            time.sleep(5)
# Wait to completion
while(True):
    if(len(processes) == 0):
        break
    else:
        new_processes = []
        for p in processes:
            if(p["p"].poll() is None):
                new_processes.append(p)
            elif(p["p"].poll() != 0):
                p["end_time"] = time.time()
                print "Error", p["p"].poll()
                print "cmd:", p["cmd"]
                print "RunTime:", (p["end_time"] - p["start_time"])
                print
                finished_processes.append(p)
            else:
                p["end_time"] = time.time()
                print "Finished ddG", p["in"], "in", (p["end_time"] - p["start_time"])
                finished_processes.append(p)
        processes = new_processes
        time.sleep(5)

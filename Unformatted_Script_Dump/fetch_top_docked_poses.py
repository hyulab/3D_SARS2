import glob
import os

from collections import defaultdict
import helper as my
import pandas as pd
import numpy as np
import tqdm


# Params
base_dir = "[REDACTED_PATH]/Collaborators/Lab_Member_Requests/Charles/2020_05_04_COVID19_Website"
orig_dir = "[REDACTED_PATH]/Collaborators/Resource_Maintenance/2020_04_27_Guided_Docking_Test/docking_outputs"
dest_dir = "{0}/Data/Docked_Structures".format(base_dir)

# Collect all docking attempts and associated energy scores
interaction2docks = defaultdict(list)
summary = []
for f in glob.glob(orig_dir + "/*.pdb"):
    p1, p2, attempt_num =  os.path.basename(f).split("_")[:3]
    
    data_dict = eval(open(f.replace("_0.pdb", ".fasc"), "r").read())
    score = data_dict["total_score"]
    
    interaction2docks[(p1, p2)].append((f, score))
    
    summary.append([p1, p2, attempt_num, score, f])



# Sort docking attempts and save the best (lowest socre) attempt
file2rank = dict()
for interaction in interaction2docks.keys():
    interaction2docks[interaction] = sorted(interaction2docks[interaction], key=lambda x: x[1])
    for i, f in enumerate(interaction2docks[interaction]):
        file2rank[f[0]] = i+1
    
    best_interaction, best_score = interaction2docks[interaction][0]
    
    os.system("cp {0}/{1} {2}/{3}_{4}_top_dock.pdb".format(orig_dir, os.path.basename(best_interaction), dest_dir, interaction[0], interaction[1]))



# Generate Summary file describing all attempts
summary = pd.DataFrame(summary, columns=["P1", "P2", "Attempt", "Score", "File"])
summary["Rank"] = summary["File"].map(lambda x: file2rank[x])
summary = summary.sort_values(["P1", "P2", "Score"])



# Calculate Ires for each attampt
pbar = tqdm.tqdm(total=len(summary))
def calc_ires(f, c1="A", c2="B"):
    try:
        pbar.update()
        ires1, ires2 = my.call("python [REDACTED_PATH]/mjm_tools/irescalc.py {0} -c1 {1} -c2 {2}".format(f, "A", "B"))
        
        return ires1, ires2
    except KeyboardInterrupt:
        raise
    except:
        return np.nan, np.nan
# FUNCTION END
tmp = summary["File"].map(calc_ires)
summary["P1_Ires"] = [x[0] for x in tmp]
summary["P2_Ires"] = [x[1] for x in tmp]



# Compare Against ECLIAR Ires
ires_summary = pd.read_csv("{0}/Data/Interface_Summary.txt".format(base_dir), sep="\t")
ires_summary = ires_summary[ires_summary["Source"] == "ECLAIR"]
interaction2ires = ires_summary.set_index(["P1", "P2"])[["P1_Ires", "P2_Ires"]].to_dict(orient="index")

def calc_stats(args):
    #print len(args)
    #print args
    p1, p2, ires1, ires2 = args
    
    # Format Sets
    if(pd.isnull(ires1)):
        ires1 = set()
    else:
        ires1 = set(ires1.split(","))
    if(pd.isnull(ires2)):
        ires2 = set()
    else:
        ires2 = set(ires2.split(","))
    
    # Fetch Eclair Ires / Format Sets
    real_ires1 = interaction2ires[(p1, p2)]["P1_Ires"]
    real_ires2 = interaction2ires[(p1, p2)]["P2_Ires"]
    
    if(pd.isnull(real_ires1)):
        real_ires1 = set()
    else:
        real_ires1 = set(real_ires1.split(","))
    if(pd.isnull(real_ires2)):
        real_ires2 = set()
    else:
        real_ires2 = set(real_ires2.split(","))
    
    # Calculate Jaccard Similarity
    j1 = np.nan
    try:
        j1 = len(ires1.intersection(real_ires1)) / float(len(ires1.union(real_ires1)))
    except ZeroDivisionError:
        pass
    
    j2 = np.nan
    try:
        j2 = len(ires2.intersection(real_ires2)) / float(len(ires2.union(real_ires2)))
    except ZeroDivisionError:
        pass
    
    # Calculate Recall
    r1 = np.nan
    try:
        r1 = len(ires1.intersection(real_ires1)) / float(len(real_ires1))
    except ZeroDivisionError:
        pass
    
    r2 = np.nan
    try:
        r2 = len(ires2.intersection(real_ires2)) / float(len(real_ires2))
    except ZeroDivisionError:
        pass
    
    
    return j1, j2, r1, r2
# FUNCTION END
tmp = summary[["P1", "P2", "P1_Ires", "P2_Ires"]].apply(calc_stats, axis=1)
summary["P1_Jaccard"] = [x[0] for x in tmp]
summary["P2_Jaccard"] = [x[1] for x in tmp]
summary["P1_Recall"] = [x[2] for x in tmp]
summary["P2_Recall"] = [x[3] for x in tmp]



# Save Summary
summary[["P1", "P2", "Attempt", "File", "Rank", "Score", "P1_Jaccard", "P1_Recall", "P1_Ires", "P2_Jaccard", "P2_Recall", "P2_Ires"]].to_csv("{0}/Data/Docking_Summary.txt".format(base_dir), sep="\t", index=None)

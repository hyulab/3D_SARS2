import pandas as pd
import glob

# Params
base_dir = "[REDACTED_PATH]/Collaborators/Lab_Member_Requests/Charles/2020_05_04_COVID19_Website"
orig_dir = "[REDACTED_PATH]/Collaborators/Eclair_Runs/2020_04_22_COVID19_Human_Interactome/Predictions"
dest_dir = "{0}/Data/Eclair_Predictions".format(base_dir)

# Collect all docking attempts and associated energy scores
summary = []
for f in glob.glob(orig_dir + "/*.pkl"):
    # Read ECLAIR Predictions
    preds = pd.read_pickle(f)
    
    # Label by confidence tier
    def label_tier(p):
        if(p < 0.12):
            return "Very Low"
        elif(p < 0.24):
            return "Low"
        elif(p < 0.36):
            return "Medium"
        elif(p < 0.48):
            return "High"
        else:
            return "Very High"
    # FUNCTION END
    preds["Tier"] = preds["Pred"].map(label_tier)
    
    # Reorder so COVID protein is P1 (for consistency with other naming)
    p1, p2 = preds[["P1", "P2"]].values[0]
    if(not "COVID" in p1):
        preds["P1"] = p2
        preds["P2"] = p1
        p1, p2 = p2, p1
        preds["Prot"] = 1 - preds["Prot"]
        preds = preds.sort_values(["Prot", "Pos"])
    
    # Save Preds to dest folder
    preds.to_csv("{0}/{1}_{2}.txt".format(dest_dir, p1, p2), sep="\t", index=None)
    
    # Store High / Very High Ires for Compiled Summary
    p1_ires = sorted(preds[(preds["Pred"] >= 0.36)&(preds["Prot"] == 0)]["Pos"].to_list())
    p2_ires = sorted(preds[(preds["Pred"] >= 0.36)&(preds["Prot"] == 1)]["Pos"].to_list())
    
    summary.append([p1, p2, "ECLAIR", sum(preds["Prot"] == 0), len(p1_ires), ",".join([str(x) for x in p1_ires]), sum(preds["Prot"] == 1), len(p2_ires), ",".join([str(x) for x in p2_ires])])
summary = pd.DataFrame(summary, columns=["P1", "P2", "Source", "P1_Len", "P1_N_Ires", "P1_Ires", "P2_Len", "P2_N_Ires", "P2_Ires"])

# Save summary file
summary.to_csv("{0}/Data/Interface_Summary.txt".format(base_dir), sep="\t", index=None)

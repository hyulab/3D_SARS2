import pandas as pd

# Uniprot Info Generated in ECLIAR
uniprot_info = pd.read_csv("/home/adr66/eclair/data/uniprot_info.txt", sep="\t")

# Original Interaciton List from Krogan Paper
interactions = pd.read_csv("/home/sdw95/Collaborators/Lab_Member_Requests/Haiyuan/2020_03_27_COVID19_3DInteractome/COVID19_Interactome.txt", sep="\t")
id2gene = interactions.set_index("Preys")["PreyGene"].to_dict() # Map Human UniProt to Prefered Gene Name

# Generate Full Set of Identifiers submitted to ECLIAR
interactions2 = pd.read_csv("/home/sdw95/Collaborators/Eclair_Runs/2020_04_22_COVID19_Human_Interactome/Interactions.txt", sep="\t", names=["P1", "P2"])
all_ids = set(interactions2["P1"].to_list() + interactions2["P2"].to_list())

# Pull out / Reformat Lines in UniProt Info we care about
protein_summary = uniprot_info[uniprot_info["id"].map(lambda x: x in all_ids)][["id", "reviewed", "genes", "protein names", "length", "sequence"]]
protein_summary["Is_Viral"] = protein_summary["id"].map(lambda x: "COVID" in x)
protein_summary["reviewed"] = protein_summary[["id", "reviewed"]].apply(lambda x: True if x[1] == "reviewed" and not "COVID" in x[0] else False, axis=1)
protein_summary["genes"] = protein_summary[["id", "genes"]].apply(lambda x: x[1] if not x[0] in id2gene else id2gene[x[0]], axis=1)

# Save
#protein_summary.sort_values(["Is_Viral", "genes"])[["id", "Is_Viral", "reviewed", "genes", "protein names", "length", "sequence"]]
protein_summary.sort_values(["Is_Viral", "genes"])[["id", "Is_Viral", "genes", "length", "sequence"]].to_csv("Data/Proteins.txt", sep="\t", header=["ID", "Is_Viral", "Gene_Name", "Length", "Sequence"], index=None)


# Read ECLAIR Domain Info
pfam_doms = pd.read_csv("/home/adr66/eclair/features/per_feature/pfam_domains.txt", names=["ID", "Is_Domain"], sep="\t")
pfam_doms = pfam_doms[pfam_doms["ID"].map(lambda x: x in all_ids)]

pfam_doms["Is_Viral"] = pfam_doms["ID"].map(lambda x: "COVID" in x)

pfam_doms.sort_values(["Is_Viral", "ID"])[["ID", "Is_Viral", "Is_Domain"]].to_csv("Data/Protein_Domains.txt", sep="\t", index=None)


# Add in COVID UniProt IDs
proteins = pd.read_csv("Data/Proteins.txt", sep="\t")

covid_fasta = my.fasta2dict("uniprot_covid_19.fasta")

covid_fasta = {k.split("|")[1]:v for k, v in covid_fasta.iteritems() if "OX=2697049" in k}

#covid2best = defaultdict(lambda: ["None", {"Pident":0}])
for uniA, seqA in tqdm_notebook(proteins[["ID", "Sequence"]].values):
    if(not "COVID19" in uniA):
        continue
    if(not "orf8" in uniA):
        continue
    for uniB, seqB in covid_fasta.iteritems():
        if(not uniB == "P0DTC8"):
            continue
        print uniA, uniB
        align = my.NWSeqAlignment(seqB, seqA)
        print align
        if(align["Pident"] > covid2best[uniB][1]["Pident"]):
            covid2best[uniB] = [uniA, align]



proteins[proteins["ID"].map(lambda x: "orf7" in x.lower())]

blacklist = ["P0DTD8", "P0DTC1", "P0DTD1"]
krogan2uni = defaultdict(str)
for uni, v in covid2best.iteritems():
    uniB, align = v
    print uni, uniB, len(align["Alignment"].replace("|", ""))
    if(not uni in blacklist):
        krogan2uni[uniB] = uni
    
    #if(align["Pident"] < 1):
    #    pass#print uni
    #    #my.alignPrint(align)

proteins["UniProt"] = proteins["ID"].map(lambda x: krogan2uni[x] if "COVID" in x else x)

from jfb_tools import batchUniProtAPI
s = proteins["UniProt"].to_list()
tmp = dict(zip(s, batchUniProtAPI(s, source_id="ACC", target_id="GENENAME")))
proteins["Gene Name"] = proteins[["ID", "UniProt"]].apply(lambda x: tmp[x[1]] if not "COVID" in x[0] else x[0].replace("COVID19", "").replace("orf9c", "orf14").replace("Spike", "S").replace("C145A", "").upper(), axis=1)

proteins.to_csv("Data/Proteins.txt", sep="\t")

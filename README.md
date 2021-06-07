# 3D_SARS2
Compilation of relevant scripts and analyses for our paper "A 3D Structural Interactome to Explore the Impact of Evolutionary Divergence, Population Variation, and Small-molecule Drugs on SARS-CoV-2-Human Protein-Protein Interactions"

For questions or clarifications please address Shayne Wierbowski at sdw95@cornell.edu

# Recommended Run Order

In the coming days we should update this repository to better organize and label which scripts should be run in which order to produce all the necessary output files in order. A work-in-progress summary is provided...

# Assorted Data Prep
1. Generate Proteins Summary
   + NEEDS TO BE RUN AFTER FETCH_ECLAIR_PREDS
   - Dependencies
     + Must have run ECLIAR pipelien
     + Currently only handles Krogan interactions (fetched some gene names from their data)
     + Requires UniProt Download for Viral Name - UniProt Mapping
   - Outputs
     + Proteins.txt

2. Select_Models - Set up interactions / select models
   + NEEDS TO BE RUN AFTER FETCH_ECLAIR_PREDS
   - Dependencies
     + Fetch_Eclair_Preds
   - Inputs
     + Interactions.txt
     + Manually created Viral Homology Models (in ECLAIR feature pipeline)
     + Manually selected / reindexed Viral PDB Structures (generated separately)
     + SIFTS / Modbase resources
   - Outputs
     + All selected Undocked Structures copied over locally
     + All structures reoriented / paired for docking
     + Models.txt

3. Fetch_Population_Variants
   - Inputs
     + Proteins.txt
     + Human_Population_Varaints_VEP_Mapped.txt (REQUIRE MANUAL VEP SUBMISSION)
   - Outputs
     + VEP_Input
     + Pop_Vars.txt

4. Compile Viral Muts
   - Inputs

   - Outputs


# Run Predictions
5.  XXX - run_eclair

6.  Fetch_Eclair_Preds
    - Dependencies
      + run_eclair
    - Inputs
       + ECLAIR .pkl files
    - Outputs
       + ECLAIR .pkl files relocated (optional)
       + Interface_Summary.txt (created with ECLAIR results to start)

7.  XXX - Run guided docking

8.  Fetch_Top_Docked_Poses
    - Dependencies
      + XXX
      + Fetch_Eclair_Preds
      + irescalc.py (sub-process)
    - Inputs
       + HADDOCK docking outputs
       + Interface_Summary.txt
       + Proteins.txt
       + Models.txt
    - Outputs
       + Top-ranked docking outputs relocated
       + Interface_Summary.txt (updated to add Docked interface entries)
       + Docking_Summary.txt
       + Docking_Summary_Charles.txt (for web purposes?)




# Run Analyses
9.  Calculate Interface Variant Enrichment
    -
    - Inputs
      + Proteins
      + Interface Summary
      + Undocked Structures
      + Pop Vars
      + Viral Muts
      + Some static resources (INSIDER, SIFTS, IRES)
    - Outputs
      + Pop_Var_Enrichments.txxt
      + Viral_Mut_Enrichments.txt

10.  Calculate Interface Similarity
    -
    - Inputs
      + Interface_Summary
      + INSIDER H_Sapiens
    - Outputs
      + Interface Similarity Viral
      + Interface Similarity Human
      + Ordered Ires Comparisons

11.  Drug Docking
    -
    - Inputs
    - Outputs
      + Krogan Drug Candidates
      + Removes Hydrogens from Undocked Structures (should be unnecessary? But also should be done when these are first created since now they could have changed between docking and now)
      + Creates docking sub-batches (N trials per drug-target pair)
      + Creates docking ranked_poses (100 top ranked per drug-target pair)
      + Createes Docked_Ligands/*.pdb (combined top 100 ranked poses)

12.  Calculate Drug Interfaces
    - Dependencies
      + Calls irescalc_ligand.py
    - Inputs
      +
    - Outputs
      + Drug Docking Ires Summary
      + Drug Interface Enrichment

13. Run SARS1 SARS2 ddG Predicitons
    - Dependencies
      + Calls Mutant_ddG_Calc.py
    - Inputs
      + Docking Summary
      + Viral Muts
    - Outputs
      + Creates ddG_Mutants/Structures and ddG_Mutants/Summary_Logs for each docked inter
      + ddG_All.txt
      + ddG_Summary.txt

14. Run Interface ddG Scanning
    - Dependencies
      + Calls Ires_ddG_Scanning.py
    - Inputs
      + Docking Summary
      + Pop vars
    - Outputs
      + Creates ddG_Single_Mutants/Structures ddG_Single_Mutants/Raw_Outputs and ddG_Mutants/Summaries for each docked inter
      + ddG_Single_Mutants/Hotspot Score Mutants


15. Calculate Disease Enrichment
    - Inputs
      + Interface Summary
      + HGMD Data
      + ClinVar Data
      + MEDGEN Relationships
      + MEDGEN Term Names
    - Outputs
      + Medgen Term Enrichments.txt

16. Prep Tables

17. Prep Downloads

# NOTE
Our paper is currently under in revision and the best way to present and format this code for ease of use may be subject to change. The current repository contains a set of Jupyter Notebooks and python scripts each intended to perform a single part of our predictions or analyses. At this stage the repository is primarily provided as a means of transparently seeing the process and code used. The exact file-structure to be able to fully re-run this code on a separate machine is not fully provided at this time. Several aspects of the resources and dependencies associated with this project are so engrained in our local system that they cannot easily be extracted.

If there are any elements of the project not adequatley covered or reproducible by this repository and that downstream users and interested in, please contact us and we will work with you make these elements available as soon as possible.

# ECLAIR Dependencies

Code related to the initial generation of ECLAIR interface features and predictions can be found in our ECLAIR repository (https://github.com/hyulab/ECLAIR). Our previous code release only included code related to training and setting up the classifiers used. We are in the process of adding the code for the whole prediction pipeline to this repository as well, but the code is similarly difficult to fully compartmentalize and release.

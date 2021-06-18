# 3D_SARS2
Compilation of relevant scripts and analyses for our paper "A 3D Structural Interactome to Explore the Impact of Evolutionary Divergence, Population Variation, and Small-molecule Drugs on SARS-CoV-2-Human Protein-Protein Interactions"

For questions or clarifications please address Shayne Wierbowski at sdw95@cornell.edu

# Installing and Running

### Disclaimer
The provided scripts are not intended to be used as a definitive pipeline or piece of generalized software, but are rather included to provide transparency as to what was done in our analyses and a potential starting point for similar analyses.

Several components of this project depend on pre-existing tools and resources built into the Yu Lab's servers and are not currently compartmentalized in a publicly accessible format. Specifically, the full ECLAIR pipeline necessary to construct features and run new interface predictions is not fully released. Our previous ECLAIR code release (https://github.com/hyulab/ECLAIR) only included code related to training and setting up the classifiers used.

### Limited Demo
As such a version of our pipeline that is fully operational from start to finish for any input is not currently possible. To provide the best possible representation of our workflow, we have included with this repository a small demo  set consisting of 10 interactions. We have provided static copies of all necessary intermediate dependencies / resources (whose fresh generation could not be generalized in this repository) so that users can work through our analysis on this small demo set.

All outputs from this demo have been included in the repository where possible, but any exceedingly large results / intermediate files have been removed. The full reproduced output from this demo is expected to require ~30 GB of space.

If as an end user you feel there are any elements of the project you are interested in but that are not adequately covered or reproducible by this repository, please contact us and we would be happy to work with you make these elements available as soon as possible.

### Dependencies and Conda Environment
This project is run in python 2.7.15. I am fairly certain that the full list of non-standard packages used is as follows (but have not explicitly tested a minimal environment to be certain)...

- biopython
- bitarray
- matplotlib
- networkx
- numpy
- pandas
- pymol
- requests
- scikit-learn
- scipy
- seaborn
- setproctitle
- tqdm

We have also included a yaml file describing the local conda environment that this project was carried out in. This environment is excessive (more packages than needed here), but should be sufficient to successfully run all of the scripts and notebooks contained here.

Additionally, several additional pieces of software must be installed...

 - HADDOCK - see installation guide [here](https://www.bonvinlab.org/software/haddock2.2/installation/)
 - PyRosetta - follow download and instructions [here](http://www.pyrosetta.org/dow)
 - NACCESS - follow instructions [here](http://www.bioinf.manchester.ac.uk/naccess/)
    - If an access key cannot be acquired [FreeSASA](https://freesasa.github.io/) may be used an as alternative. But this will require re-working the srescalc.py and irescalc.py scripts as well as corresponding functions in the helper_functions.py.
       -  (we have developed but not fully tested a draft version using FreeSASA, contact us if necessary)
 - Smina
	 - Can be installed through conda with the command...
        - > conda install -c conda-forge smina
 - OpenBabel  
	 - Can be installed through conda with the command...
        - > conda install -c conda-forge smina
  - PyMol
	 - Can be installed through conda with the command...
        - > conda install -c schrodinger pymol

# Recommended Run Order
There is no exactly correct way to run through our provided demo, but the notebooks have been numbered in a recommended run order. Each separate notebook has been divided to perform a single task (or series of related tasks) so that the code for any particular step in our analysis can easily be identified.

The individual Notebooks alongside any inputs, outputs, dependent resources, and code pre-requisites are described below...

## 01_Run_ECLAIR
This notebook contains the command that WOULD run the ECLAIR pipeline on the provided set of interactions. The full ECLAIR pipeline is not currently provided as a standalone piece of software. For demonstrative purposes, the raw ECLAIR prediction outputs for a handful of interactions are provided (so that the rest of the pipeline is here is runnable in some form).

- Inputs:
  - Interactions.txt


- Outputs:
  - [P1]\_[P2]\_[prot].pkl (Raw Eclair Prediction outputs)
  - Updates several feature files some of which are used downstream and are provided as static copies. Specifically...
    - uniprot_info.txt
    - pfam_domains.txt


- Dependencies:
  - Should be run AFTER all interactions have been fed through ECLAIR pipeline for interface prediction
    - **NOTE:** The ECLAIR pipeline is not included in this repository, and treats any output from this pipeline as a static result that is already available

## 02_Fetch_ECLAIR_Preds
This notebook reads in all raw outputs from the ECLAIR pipeline, parses them, and resaves them as plain text files. Additionally generates the initial data for the Interface_Summary.txt.

Raw ECLAIR outputs are saved as .pkl objects and separated out as [P1]\_[P2]\_0.pkl and [P1]\_[P2]\_1.pkl (for the interface predictions on the P1-P2 interaction for P1 and P2 respectively). The parsing done here merges the two files into one and additionally adds a confidence tier interpretation to the raw prediction values (only High and Very High) predictions are
retained in the binary interface prediction output.

- Inputs:
  - [P1]\_[P2]\_[prot].pkl (Raw Eclair Prediction outputs)


- Outputs:
  - [P1]\_[P2].txt (Parsed Eclair Prediciton outputs)
  - Interface_Summary.txt


- Dependencies:
  - Should be run AFTER all interactions have been fed through ECLAIR pipeline for interface prediction (01_Run_ECLIAR.ipynb)
    - **NOTE:** The ECLAIR pipeline is not included in this repository, and treats any output from this pipeline as a static result that is already available

## 03_Generate_Proteins_Summary
This notebook generates a summary of all of the proteins included in the input interaction set and outputs it as Proteins.txt. The primary purpose of this is just to have convenient access to gene names / sequences / lengths for all human UniProt / SARS-CoV-2 protein entries for later analyses and result formatting.

- Inputs:
  - Interactions.txt
  - Covid19_Interactome.txt


- Static Resource Dependencies:
  - uniprot_info.txt
  - pfam_domains.txt
  - uniprot_covid_19.fasta


- Outputs:
  - Proteins.txt
  - Protein_Domains.txt


- Dependencies:
  - Should be run AFTER all interactions have been fed through ECLAIR pipeline for interface prediction (01_Run_ECLIAR.ipynb)
    - **NOTE:** The ECLAIR pipeline is not included in this repository, and treats any output from this pipeline as a static result that is already available

## 04_Select_Models
This notebook selects the best available homology model or PDB structure to use for protein-protein docking for each interaction and saves them in the Output/Undocked_Structures. It additionally creates Models.txt summarizing the selected models for all proteins. A set of Oriented_Structures with their most likely interface sides oriented towards each other (and separated 5 angstroms) is also produced, but is not currently used in any downstream analysis.

- Inputs:
  - Interactions.txt
  - [P1]\_[P2].txt (Parsed Eclair Prediciton outputs)
  - [P1].pdb (SARS-CoV-2 Homology Models)
  - [P1].pdb (SIFTS Mapped SARS-CoV-2 PDB Structures)


- Static Resource Dependencies:
  - pdbresiduemapping.txt


- Outputs:
  - [Prot]\_[Source].pdb (Undocked Structures)
  - [P1]\_[P2]\_[Chain].pdb (Oriented Structures)
  - Models.txt


- Dependencies:
  - Must be run after 02_Fetch_Eclair_Preds
  - Assumes that SARS-CoV-2 Homology Models and / or PDB Structures have already been prepared

## 05_Fetch_Population_Variants
This notebook fetches all known population variants for the human interactors included in Proteins.txt using the GnomAD API. Then (requiring some manual input) these variants are mapped through the Variant Effect Predictor (VEP) to get corresponding UniProt positions / other information.

**NOTE:** When returning to this step to reformat the code for this repository I encountered a lot of errors because some of the external API's and tools that I used seem to have changed since I originally ran the code. I've reproduced the original workflow as best as possible, but still encounter some troubles I did not originally encounter with genes failing the GnomAD API query. The example output from this notebook included in this repository is just a sub-set of my original population variants list for the human genes included in the demo interactions list. Results could vary if you try to reproduce this from scratch.

- Inputs:
  - Proteins.txt


- Outputs:
  - VEP_Input.txt
  - Human_Population_Variants_VEP_Mapped.txt
  - Pop_Vars.txt


- Dependencies:
  - Must be run after 03_Generate_Proteins_Summary

## 06_Compile_Viral_Mutations
This notebook compares the sequences of all SARS-CoV-2 proteins with their SARS-CoV homologs where available and compiles a list of all sequence deviations (Viral_Muts.txt)

- Inputs:
  - Proteins.txt


- Outputs:
  - Viral_Muts.txt


- Dependencies:
  - Must be run after 03_Generate_Proteins

## 07_Run_PPI_Docking
This notebook is a wrapper to call the HADDOCK protein-protein interaction docking protocol for all interactions.

- Inputs:
  - Interactions.txt
  - [P1]\_[P2]\_[Chain].pdb (Oriented Structures)
  - [P1]\_[P2].txt (Parsed ECLAIR Predictions)


- Outputs:
  - [P1]\_[P2] (Interaction Haddock Run Directory created under "Docking_Runs")


- Dependencies:
  - Must be run after 02_Fetch_Eclair_Preds and 04_Select_Models
  - Must have HADDOCK installed locally
  - Calls run_haddock.py
    - Which itself calls srescalc.py
      - **NOTE:** srescalc.py *may not* be currently properly extracted from the Yu Lab's server and may not run successfully in this repository. The raw code is provided, but it itself calls several separate dependencies, and I have not been able to thoroughly confirm there are no specifics to our machine still linked to it.
      - I *believe* it should be functional, but if any end user encounters errors running srescalc.py from this repository please contact the authors.
      - Requires NACCESS installed locally

## 08_Fetch_Top_Docks
This notebook selects the top-scored docking output from each of the HADDOCK runs to use in downstream analysis. It additionally calculates the interface residues from these top-scored docking outputs, compares them against the original ECLAIR predictions, and produces the Docking_Summary.txt. It also updates the Interface_Summary.txt to add interface annotations from the docking results.

**NOTE:** At this stage AFTER docking is already run, docked outputs are filtered to only retain for downstream analysis docked interactions whose structures had sufficient coverage, or included a high-confidence interface prediction that was used to guide docking. Practically speaking this filtering step should instead be applied before docking is run at all to save time.

- Inputs:
  - [P1]\_[P2] (Interaction Haddock Run Directory created under "Docking_Runs")
  - Interface_Summary.txt
  - Proteins.txt
  - Models.txt


- Outputs:
  - [P1]\_[P2]_top_dock.pdb
  - Docking_Summary
  - Interface_Summary.txt (updated)


- Dependencies:
  - Must be run after 07_Run_PPI_Docking and 03_Generate_Proteins
  - Calls irescalc.py
    - **NOTE:** irescalc.py *may not* be currently properly extracted from the Yu Lab's server and may not run successfully in this repository. The raw code is provided, but it itself calls several separate dependencies, and I have not been able to thoroughly confirm there are no specifics to our machine still linked to it.
    - I *believe* it should be functional, but if any end user encounters errors running irescalc.py from this repository please contact the authors.
    - Requires NACCESS installed locally

## 09_Calculate_Interface_Variation_Enrichment
This notebook calculates the Log Odds enrichment for occurrence of human population variants or SASR-CoV-1 to SARS-CoV-2 sequence divergences along the ECLAIR-predicted or docked interfaces. Output is summarized in Pop_Var_Enrichments.txt and Viral_Mut_Enrichments.txt.

- Inputs:
  - Interface_Summary.txt
  - Pop_Vars.txt
  - Viral_Muts.txt
  - Proteins.txt


- Static Resource Dependencies:
  - H_sapiens_interfacesAll.txt
  - pdbresiduemapping.txt
  - ires_perpdb_alltax.txt


- Outputs:
  - Pop_Var_Enrichments.txt
  - Viral_Mut_Enrichments.txt


- Dependencies:
  - Must be run after 05_Fetch_Population_Variants, 06_Compile_Viral_Mutations, and 08_Fetch_Top_Docks

## 10_Calculate_Interface_Similarity
This notebook calculates the similarity (Jaccard index) between interface annotations on different interactions involving the same protein. Interface similarities in viral proteins are calculated relative to the interfaces predicted for other human interactors of the same viral protein and are summarized in Viral_Interface_Similarities.txt (minimally applicable in the example set selected for this demo). Interface similarities in human proteins are calculated relative to known and predicted interfaces for other human-human interactions involving the same human protein and are summarized in Human_Interface_Similarities.txt.

The primary purpose of these similarities are for ordering and displaying similarity to other interfaces on our web site (ordering selected based on Ordered_Ires_Comparisons.txt). No further analysis was done based on similarity between interfaces.

- Inputs:
  - Interface_Summary.txt
  - Proteins.txt


- Static Resource Dependencies:
  - H_sapiens_interfacesHQ.txt


- Outputs:
  - Viral_Interface_Similarities.txt
  - Human_Interface_Similarities.txt
  - Ordered_Ires_Comparisons.txt


- Dependencies:
  - Must be run after 08_Fetch_Top_Docks

## 11_Calculate_Disease_Enrichment
This notebook calculates the similarity log odds ratios for enrichment of disease mutations on human interactors of SARS-CoV-2 and on specific interfaces (either with SARS-CoV-2 proteins or with other human proteins). A summary of the most enriched disease terms in the human interactor set is retained (MedGen_Term_Enrichments.txt).


- Inputs:
  - Interface_Summary.txt
  - Proteins.txt


- Static Resource Dependencies:
  - HGMD_201801_uniprotmapped.txt
    - **NOTE:** Database publicly available without, a blank header describing the layout is provided instead
  - ClinVar2017_missense_uniprotmapped.txt
  - H_sapiens_interfacesALL.txt
  - keep_uniprot.fasta
  - MGREL.RRF
    - Downloaded on first run
  - MGCONSO.RRF
    - Downloaded on first run


- Outputs:
  - MedGen_Term_Enrichments.txt


- Dependencies:
  - Must be run after 08_Fetch_Top_Docks

## 12_Run_SARS1_SARS2_ddG_Predictions
This notebook is a wrapper to that calls a PyRosetta protocol to predict SARS1 --> SARS2 ddG value for all docked interactions. The prediction is based on using the Rosetta energy function to model the change in energy between the bound and unbound forms of the SARS-CoV and SARS-CoV-2 interaction. SARS-CoV docked structures are generated by applying all of the SARS1 --> SARS2 sequence deviations (in reverse) to the SARS-CoV-2 docked structure.


- Inputs:
  - Docking_Summary.txt
  - Viral_Muts.txt
  - [P1]\_[P2] (Interaction Haddock Run Directory created under "Docking_Runs")


- Outputs:
  - [P1]\_[P2]  (directory and summary logs per interaction created under ddG_Mutated_Structures)
  - ddG_All.txt
  - ddG_Summary.txt


- Dependencies:
  - Must be run after 06_Compil_Viral_Mutations and 07_Run_PPI_Docking
  - Calls Mutant_ddG_Calc.py

## 13 Run_Interface_ddG_Scanning
This notebook is a wrapper to that calls a PyRosetta perform interface scanning mutagenesis across all docked interactions to predict the ddG impact of all possible mutations along the interface.


- Inputs:
  - Docking_Summary.txt
  - Pop_Vars.txt
  - [P1]\_[P2] (Interaction Haddock Run Directory created under "Docking_Runs")


- Outputs:
  - [P1]\_[P2]  (directory and summary logs per interaction created under ddG_Single_Mutants)
  - Hotspot_Scored_Mutants.txt


- Dependencies:
  - Must be run after 05_Fetch_Population_Variants and 07_Run_PPI_Docking
  - Calls Ires_ddG_Scanning.py

## 14_Run_Drug_Docking
This notebook first parses a provided list of candidate drugs that target one or more human interactors of SARS-CoV-2 (this is a static input, but this setup could theoretically be generalized to any list of protein-drug pairs). Then protein-ligand docking for each pair is run in smina and the top-ranked docked conformation from each pair is retained.


- Inputs:
  - Krogan_Drug_Candidates.txt
  - COVID_19_Interactome.txt
  - Models.txt
  - Proteins.txt
  - [Prot]\_[Source].pdb (Undocked Structures)


- Outputs:
  - [Drug].pdb (Undocked Drug Structures)
  - [Drug].svg (2D Image of Drug Structure)
  - [P1]\_[Drug].pdb (top ranked drug docking for each protein-drug pair)


- Dependencies:
  - Must be run after 04_Select_Models
  - Must have smina installed locally (available through conda)
  - Must have openbabel installed locally (available through conda)

## 15
This notebook calculates the interface residues from all drug-protein docked outputs (summarized in Drug_Docking_Ires_Summary.txt). It then calculates the log odds enrichment for co-occurence of ligand binding site residues for the drug-protein pair with interface residues for the human-viral protein-protein interaction (summarized in Drug_Interface_Enrichment.txt).


- Inputs:
  - Krogan_Drug_Candidates.txt
  - [P1]\_[Drug].pdb (top ranked drug docking for each protein-drug pair)
  - [Prot]\_[Source].pdb (Undocked Structures)
  - Proteins.txt


- Outputs:
  - Drug_Docking_Ires_Summary.txt
  - Drug_Interface_Enrichment.txt


- Dependencies:
  - Must be run after 14_Run_Drug_Docking
  - Calls irescalc_ligand.py
    - **NOTE:** irescalc_ligand.py *may not* be currently properly extraced from the Yu Lab's server and may not run successfully in this repository. The raw code is provided, but it itself calls several separate dependencies, and I have not been able to thoroughly confirm there are no specifics to our machine still linked to it.
    - I *believe* it should be functional, but if any end user encounters errors runngin irescalc_ligand.py from this repository please contact the authors.
    - Requires NACCESS installed locally
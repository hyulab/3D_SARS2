# 3D_SARS2
Compilation of relevant scripts and analyses for our paper "A 3D Structural Interactome to Explore the Impact of Evolutionary Divergence, Population Variation, and Small-molecule Drugs on SARS-CoV-2-Human Protein-Protein Interactions"

# NOTE
Our paper is currently under review. To provide immediate code transparency, we have temporarilly provided all raw scripts used in completing this paper. Most of the work is done in one large interactive jupyter notebook. We plan to reformat and reorganize this code into a series of executable scripts. We are currently in the process of editing, commenting, and reorganizing these scripts into a format appropraite for a public code repository. These formatted scripts will be made available as soon as possible. At the latest we will update this repository with formatted scripts during a revision process.

Most of the computational experiments and analyses were set up and run in the Interactive_Testing.ipynb jupyter notebook or were run using the supplemental python scripts (e.g. Docking_Script.py, Mutant_ddG_Calc.py, or Ires_ddG_Scanning.py). If extensive review of this code is desired or necessary before we are able to provide an organized set of scripts, these would be the most relevant files. Some of the setup for the initial structures used for guided docking was done in the Interactive_Testing.ipynb, but the actual docking experiments were called using the Guided_Docking_Test.ipynb notebook and through the Docking_Script.py.

The Prep_*.ipynb notebooks were used to reformat data / analyses relevant to particular items for the manuscript (Figures, Tables, additional Downloads), and are largely not important to understanding the methods and experiments performed here.

The helper.py script is largely unreleated to this project, but contains some custom functions that were imported and recyled here (mainly for parsing and working with pdb files). It is included in this repository so that methods from this library can be explored when encountered in the other notebooks if necessary.

The Docking_Script.py, Mutant_ddG_Calc.py, and Ires_ddG_Scanning.py are modified versions of the provided PyRosetta pipelines for Docking and Ala_Scanning provided in the documentation.

Code related to the initial generation of ECLAIR interface features and predictions can be found in our ECLAIR repository (https://github.com/hyulab/ECLAIR). Our previous code release only included code related to training and setting up the classifiers used. We are in the process of adding the code for the whole prediction pipeline to this repository as well.

# This brief script is just used to call a provided selection command in PyMol
# on a provided PDB file, and return the list of residues included in the selection.
#
# This is called as a sub-script because in my experience its easiest and least error-prone
# to implement PyMol functions as their own standalone script. But this could theoretically
# also be injected as a function in the code calling it.
#
# Used specifically by run_haddock.py to determine nearby passive residues for AIR defining.

import __main__
__main__.pymol_argv = ['pymol','-cqk'] # Pymol: quiet and no GUI

import pymol
pymol.finish_launching()
from pymol import cmd as pymolCmd

import sys

pdb = sys.argv[1]
sele = sys.argv[2]

pymolCmd.load(pdb, "obj")
pymolCmd.select("tmp", sele)

mylist = []
myspace = {'mylist': mylist}
pymolCmd.iterate("tmp", 'mylist.append(resi)', space=myspace)

print ",".join(set(mylist))

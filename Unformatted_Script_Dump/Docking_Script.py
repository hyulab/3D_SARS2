#!usr/bin/env python

from __future__ import print_function

################################################################################
# A GENERAL EXPLANATION

"""
docking.py

This script performs the general Rosetta DockingProtocol, for predicting the
conformations of two binding proteins. Rosetta docking consists of a low-
resolution (centroid) stage of many random rigid-body (fixed conformation)
perturbations to locate interfaces and a high-resolution (fullatom) stage
consists of smaller rigid-body perturbations, sidechain packing, and
minimization to accurately predict the protein conformations.

Instructions:

1) ensure that your PDB file is in the current directory
2) run the script:
    from commandline                        >python D100_Docking.py

    from within python/ipython              [1]: run D100_Docking.py

Author: Evan H. Baugh
    based on an original script by Sid Chaudhury
    revised and motivated by Robert Schleif

Updated by Boon Uranukul, 6/9/12
Simplified special constant seed initialization ~ Labonte

References:
    J. J. Gray, "High-resolution protein-protein docking," Curr. Opinions in
        Struct. Bio. 16 (2) 183-193 (2006).

"""

################################################################################
# THE BASIC PROTOCOL, sample_docking

"""
This sample script is setup for usage with
    commandline arguments,
    default running within a python interpreter,
    or for import within a python interpreter,
        (exposing the sample_docking method)

The method sample_docking:
1.  creates a pose from the desired PDB file
2.  sets up the pose FoldTree for docking
3.  creates Movers for switching between fullatom and centroid and for
        recovering the original sidechain conformations of the fullatom pose
4.  converts the pose to centroid
5.  creates a (centroid) copy of the pose to be modified
6.  creates relevant ScoreFunctions
        -for centroid docking
        -for fullatom docking
        -for fullatom docking minimization
7.  sets up perturbation Movers:
        -RigidBodyRandomizeMovers
        -RigidBodyPerturbMover
        -RigidBodySpinMover
        -DockingSlideIntoContactMover
8.  sets up a MinMover to minimize on the docking jump
9.  creates a SequenceMover for the perturbation step
10. sets up the DockingProtocol
11. creates a (Py)JobDistributor for managing multiple trajectories
12. create a PyMOL_Observer for viewing intermediate output
13. perform protein-protein docking:
        a. set necessary variables for the new trajectory
            -reset the pose structure to the input (centroid) conformation
            -change the pose's PDDInfo.name, for exporting to PyMOL
        b. perturb the starting structure
            >randomize the upstream partner
            >randomize the downstream partner
            >spin the pose
            >slide the docking partners into contact
            >convert to fullatom
            >(fullatom) minimize the score on the inter-chain distance (jump)
        c. perform Rosetta docking
        d. output the decoy structure
            -to PyMOL using the PyMOL_Observer.pymol
            -to a PDB file using the PyJobDistributor.output_decoy

"""

import optparse    # for sorting options

from rosetta import *
from pyrosetta import *
from rosetta.protocols.rigid import *
from pyrosetta.rosetta import protocols

init(extra_options = "-constant_seed")

import os; os.chdir("[REDACTED_PATH]Collaborators/Resource_Maintenance/2020_04_27_Guided_Docking_Test/docking_outputs")

#########
# Methods

def sample_docking(pdb_filename, partners,
        translation = 3.0, rotation = 8.0,
        jobs = 1, job_output = 'dock_output'):
    """
    Performs protein-protein docking using the Rosetta standard DockingProtocol
        on the proteins in  <pdb_filename>  using the relative chain
        <partners>  with an initial perturbation using  <translation>
        Angstroms and  <rotation>  degrees.  <jobs>  trajectories are performed
        with output structures named  <job_output>_(job#).pdb.
        structures are exported to a PyMOL instance.

    """
    # 1. creates a pose from the desired PDB file
    pose = Pose()
    pose_from_file(pose, pdb_filename)

    # 2. setup the docking FoldTree
    # using this method, the jump number 1 is automatically set to be the
    #    inter-body jump
    dock_jump = 1
    # the exposed method setup_foldtree takes an input pose and sets its
    #    FoldTree to have jump 1 represent the relation between the two docking
    #    partners, the jump points are the residues closest to the centers of
    #    geometry for each partner with a cutpoint at the end of the chain,
    # the second argument is a string specifying the relative chain partners
    #    such as "A_B" of "LH_A", ONLY TWO BODY DOCKING is supported and the
    #    partners MUST have different chain IDs and be in the same pose (the
    #    same PDB), additional chains can be grouped with one of the partners,
    #    the "_" character specifies which bodies are separated
    # the third argument...is currently unsupported but must be set (it is
    #    supposed to specify which jumps are movable, to support multibody
    #    docking...but Rosetta doesn't currently)
    # the FoldTrees setup by this method are for TWO BODY docking ONLY!
    protocols.docking.setup_foldtree(pose, partners, Vector1([dock_jump]))

    # 3. create centroid <--> fullatom conversion Movers
    to_centroid = SwitchResidueTypeSetMover('centroid')
    to_fullatom = SwitchResidueTypeSetMover('fa_standard')
    # and a Mover to recover sidechain conformations
    #    when a protocol samples backbone torsion space in centroid,
    #    the sidechain conformations are neglected, when it is transferred
    #    to fullatom, we typically set the sidechain conformations to their
    #    "original" values and perform sidechain packing,
    #    a ReturnSidechainMover saves a pose's sidechains (in this case
    #    staring_pose) and when applied, inserts these conformations
    #    into the input pose
    recover_sidechains = protocols.simple_moves.ReturnSidechainMover(pose)

    # 4. convert to centroid
    to_centroid.apply(pose)

    # 5. create a (centroid) test pose
    test_pose = Pose()
    test_pose.assign(pose)

    # 6. create ScoreFunctions for centroid and fullatom docking
    scorefxn_low = create_score_function('interchain_cen')
    scorefxn_high = create_score_function('docking')

    # PyRosetta3: scorefxn_high_min = create_score_function_ws_patch('docking', 'docking_min')
    scorefxn_high_min = create_score_function('docking', 'docking_min')

    # 7. create Movers for producing an initial perturbation of the structure
    # the DockingProtocol (see below) can do this but several Movers are
    #    used to demonstrate their syntax
    # this Mover uses the axis defined by the inter-body jump (jump 1) to move
    #    the docking partners close together
    slide_into_contact = protocols.docking.DockingSlideIntoContact(dock_jump)

    # 8. setup the MinMover
    # the MoveMap can set jumps (by jump number) as degrees of freedom
    movemap = MoveMap()
    movemap.set_jump(dock_jump, True)
    # the MinMover can minimize score based on a jump degree of freedom, this
    #    will find the distance between the docking partners which minimizes
    #    the score
    minmover = protocols.minimization_packing.MinMover()
    minmover.movemap(movemap)
    minmover.score_function(scorefxn_high_min)

    # 9. create a SequenceMover for the perturbation step
    perturb = protocols.moves.SequenceMover()
    perturb.add_mover(slide_into_contact)
    perturb.add_mover(to_fullatom)
    perturb.add_mover(recover_sidechains)
    perturb.add_mover(minmover)

    # 10. setup the DockingProtocol
    # ...as should be obvious by now, Rosetta applications have no central
    #    standardization, the DockingProtocol object can be created and
    #    applied to perform Rosetta docking, many of its options and settings
    #    can be set using the DockingProtocol setter methods
    # here, on instance is created with all default values and the movable jump
    #    is manually set to jump 1 (just to be certain), the centroid docking
    #    ScoreFunction is set and the fullatom docking ScoreFunction is set
    dock_prot = protocols.docking.DockingProtocol()    # contains many docking functions
    dock_prot.set_movable_jumps(Vector1([1]))    # set the jump to jump 1
    dock_prot.set_lowres_scorefxn(scorefxn_low)
    dock_prot.set_highres_scorefxn(scorefxn_high_min)
    #### you can alternatively access the low and high resolution sections of
    ####    the DockingProtocol, both are applied by the DockingProtocol but
    ####    a novel protocol may only require centroid (DockingLowRes) or
    ####    fullatom (DockingHighRes), uncomment the lines below and their
    ####    application below
    #docking_low = DockingLowRes()
    #docking_low.set_movable_jumps(Vector1([1]))
    #docking_low.set_scorefxn(scorefxn_low)
    #docking_high = DockingHighRes()
    #docking_high.set_movable_jumps(Vector1([1]))
    #docking_high.set_scorefxn(scorefxn_high)

    # 11. setup the PyJobDistributor
    jd = PyJobDistributor(job_output, jobs, scorefxn_high)
    temp_pose = Pose()    # a temporary pose to export to PyMOL
    temp_pose.assign(pose)
    to_fullatom.apply(temp_pose)    # the original pose was fullatom
    recover_sidechains.apply(temp_pose)    # with these sidechains
    jd.native_pose = temp_pose    # for RMSD comparison

    # 12. setup a PyMOL_Observer (optional)
    # the PyMOL_Observer object owns a PyMOLMover and monitors pose objects for
    #    structural changes, when changes are detected the new structure is
    #    sent to PyMOL
    # fortunately, this allows investigation of full protocols since
    #    intermediate changes are displayed, it also eliminates the need to
    #    manually apply the PyMOLMover during a custom protocol
    # unfortunately, this can make the output difficult to interpret (since you
    #    aren't explicitly telling it when to export) and can significantly slow
    #    down protocols since many structures are output (PyMOL can also slow
    #    down if too many structures are provided and a fast machine may
    #    generate structures too quickly for PyMOL to read, the
    #    "Buffer clean up" message
    # uncomment the line below to use the PyMOL_Observer
    protocols.moves.AddPyMOLObserver(test_pose, True)

    # 13. perform protein-protein docking
    counter = 0    # for pretty output to PyMOL
    while not jd.job_complete:
        # a. set necessary variables for this trajectory
        # -reset the test pose to original (centroid) structure
        test_pose.assign(pose)
        # -change the pose name, for pretty output to PyMOL
        counter += 1
        test_pose.pdb_info().name(job_output + '_' + str(counter))

        # b. perturb the structure for this trajectory
        perturb.apply(test_pose)

        # c. perform docking
        dock_prot.apply(test_pose)
        #### alternate application of the DockingProtocol pieces
        #docking_low.apply(test_pose)
        #docking_high.apply(test_pose)

        # d. output the decoy structure
        to_fullatom.apply(test_pose)    # ensure the output is fullatom
        # to PyMOL
        test_pose.pdb_info().name(job_output + '_' + str( counter ) + '_fa')
        # to a PDB file
        jd.output_decoy(test_pose)
        
        # Manually break after first iteration?
        break

################################################################################
# INTERPRETING RESULTS

"""
The (Py)JobDistributor will output the lowest scoring pose for each trajectory
(as a PDB file), recording the score in <job_output>.fasc. Generally,
the decoy generated with the lowest score contains the best prediction
for the protein conformation. PDB files produced from docking will contain
both docking partners in their predicted conformation. When inspecting these
PDB files (or the PyMOL_Observer output) be aware that PyMOL can introduce or
predict bonds that do not exist, particularly for close atoms. This rarely
occurs when using the PyMOLMover.keep_history feature (since PyRosetta will
sample some conformation space that has clashes).

The PyMOL_Observer will output a series of structures directly produced by the
DockingProtocol. Unfortunately, this includes intermediate structures during the
application of SwitchResidueTypeSetMovers (very pretty to watch...but somewhat
useless). A LARGE number of structures are output to PyMOL and your
machine may have difficulty loading all of these structures. If this occurs,
try changing the PyMOL_Observer keep_history to False or running the
protocol without the PyMOL_Observer.

"""

################################################################################
# COMMANDLINE COMPATIBILITY

# everything below is added to provide commandline usage,
#   the available options are specified below
# this method:
#    1. defines the available options
#    2. loads in the commandline or default values
#    3. calls sample_docking with these values

# parser object for managing input options
# all defaults are for the example using "test_in.pdb" with reduced
#    cycles/jobs to provide results quickly
parser = optparse.OptionParser()
parser.add_option('--pdb_filename', dest = 'pdb_filename',
    default = '../test/data/test_dock.pdb',    # default example PDB
    help = 'the PDB file containing the proteins to dock')
# for more information on "partners", see sample_docking step 2.
parser.add_option('--partners', dest = 'partners',
    default = 'E_I',    # default for the example test_dock.pdb
    help = 'the relative chain partners for docking')
# perturbation options
parser.add_option('--translation', dest = 'translation',
    default = '3.0',
    help = 'magnitude of the random translation applied (in Angstroms)')
parser.add_option('--rotation', dest = 'rotation',
    default = '8.0',
    help = 'magnitude of the random rotation applied (in degrees)')
# PyJobDistributor options
parser.add_option( '--jobs', dest='jobs' ,
    default = '1',    # default to single trajectory for speed
    help = 'the number of jobs (trajectories) to perform')
parser.add_option( '--job_output', dest = 'job_output',
    default = 'dock_output',    # if a specific output name is desired
    help = 'the name preceding all output, output PDB files and .fasc')
(options,args) = parser.parse_args()

# PDB file option
pdb_filename = options.pdb_filename
partners = options.partners
# perturbation options
translation = float(options.translation)
rotation = float(options.rotation)
# JobDistributor options
jobs = int(options.jobs)
job_output = options.job_output

sample_docking(pdb_filename, partners, translation, rotation,
    jobs, job_output)

################################################################################
# ALTERNATE SCENARIOS

################
# A Real Example
"""
All of the default variables and parameters used above are specific to
the example with "test_in.pdb", which is supposed to be simple,
straightforward, and speedy. Here is a more practical example:

The influenza protein NS1 is implicated to increase virulence. Suppose you are
curious about the molecular mechanism of binding for this protein and decide
to investigate using PyRosetta.

1. Download a copy of RCSB PDB file 3RT3 (remove waters and any other HETATM)
2. Make a directory containing:
        -the PDB file for 3RT3 (cleaned of HETATMs and waters)
            lets name it "3RT3.clean.pdb" here
        -this sample script (technically not required, but otherwise the
            commands in 3. would change since docking.py wouldn't be here)
3. Run the script from the commandline with appropriate arguments:

>python docking.py --pdb_filename 3RT3.clean.pdb --partners B_C --jobs 400 --job_output 3RT3_docking_output --translation 3 --rotation 8 --PyMOLMover_ip off

        -The partners option, "B_C" is PDB specific, if you change the chain
            IDs in 3RT3, make sure this matches
        -400 trajectories is low, sampling docking conformations is difficult,
            typically thousands of (800-1000) trajectories are attempted
        -The common values for translation and rotation are 3 (Angstroms) and
            8 (degrees) respectively, since initial configurations are
            randomized, the rotation should not significantly alter the sampling
            space, if the translation is too big or small,it won't provide a
            useful move since minimization will move the structures to a notably
            different conformation (either due a clash or sliding into contact)
        -This script features the PyMOL_Observer, and consequently a commandline
            option ("off") to avoid using it, Monte Carlo simulations are not
            expected to produce kinetically meaningful results and as such,
            viewing the intermediates is only useful when understanding a
            protocol and rarely produces insight beyond the final output

4. Wait for output, this will take a while (performing 400 trajectories
        of the full DockingProtocol is intensive)
5. Analyze the results (see INTERPRETING RESULTS above)

Note: this is a direct port of the Rosetta docking protocol and provides example
syntax for using this method within Python. This script also provides example
code for PyRosetta. A priori, a large number of jobs (~1000 or more) is required
to achieve useful results. The best protocols are somewhat protein-specific and
the scoring and sampling methods here may be customized to produce better
results on your protein.

"""

##############################
# Changing Docking Sampling
"""
Monte Carlo sampling of docking is abstractly straightforward however the
specific values required to yield good predictions (rotation, translation,
scoring weights) are elusive. The process involves iterative steps of rigid-
body docking (shape complementarity, determining which shapes can fit) and
structure prediction (how the new interface context changes the protein
conformation). Although the DockingProtocol is nicely packaged as a single
application, object effectively performs rounds of perturbation (translation,
rotation, distance minimization) and optimization (sidechain packing, structural
refinement).

Many other approaches to docking exist however PyRosetta exposes a myriad of
tools for constructing novel Monte Carlo algorithms for docking prediction.

Please try alternate sampling methods to better understand how these
algorithms perform and to find what moves best suite your problem.

"""

#############################
# Changing Docking Scoring
"""
The Rosetta docking and docking_min score functions have weights optimized for
performance. Other scoring terms or weights may provide useful biasing for
docking. Literature on this docking method suggests that the implicit
solvation scores drive the formation of protein-protein interfaces in Rosetta
while scores representing van der Waals repulsive forces prevent structures
from collapsing into one another.

Please try alternate scoring functions or unique selection methods to better
understand which scoring terms contribute to performance and to find what
scoring best suites your problem.

"""

###############
# Other Docking
"""
Binding interaction drive MANY biological processes. While protein-protein
docking addresses many possibilities, there are other compounds which a user
may wish to investigate. Rosetta (and thus the core methods exposed in
PyRosetta) provides protocols for predicting docking with other chemical
entities such as DNA, RNA (see dna_docking.py), and small ligands (see
ligand_docking.py). Rosetta provides symmetry tools for predicting symmetric
protein complexes (another kind of protein-protein interaction) but these
methods are not currently supported in PyRosetta. Peptide (or other flexible
molecule) docking does not currently have a fine tuned protocol however
PyRosetta provides all the tools necessary for such algorithms. When binding
a highly flexible molecule, the sampling space becomes much larger,
compounding the already large search space of docking and conformation
prediction. Such algorithms are anticipated to require extremely intensive
sampling to achieve useful results. Antibody docking is supported by Rosetta
is likewise confounded by a very large search space. Small molecule docking
covers a large variety of biologically relevant interactions including
(potentially) the placement of structurally important water molecules.

Problem                    Rosetta tool        In PyRosetta?
protein-protein docking    DockingProtocol     yes
small molecule docking     LigandDockProtocol  components, use DockMCMProtocol
nucleic acid interface     DNAInterfaceFinder  components, use DockMCMProtocol
symmetric complexes        SymDockProtocol     no
flexible polymer docking   -                   no
antibody docking           AntibodyModeler     no
water placement            -                   tools, no direct protocol

"""

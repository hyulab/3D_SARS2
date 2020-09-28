# the Interface object is required by this script
from rosetta.protocols.scoring import Interface
from rosetta import *
from pyrosetta import *
from pyrosetta.rosetta import protocols

init()  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!

# returns the "interaction energy" for the pose with a given mutation
def interface_ddG(pose, muts, movable_jumps, scorefxn='', cutoff=8.0, out_filename=''):
    # 1. create a reference copy of the pose
    wt = Pose()    # the "wild-type"
    wt.assign(pose)
    
    # 2. setup a specific default ScoreFunction
    if not scorefxn:
        # this is a modified version of the scoring function discussed in
        #    PNAS 2002 (22)14116-21, without environment dependent hbonding
        scorefxn = ScoreFunction()
        scorefxn.set_weight(fa_atr, 0.44)
        scorefxn.set_weight(fa_rep, 0.07)
        scorefxn.set_weight(fa_sol, 1.0)
        scorefxn.set_weight(hbond_bb_sc, 0.5)
        scorefxn.set_weight(hbond_sc, 1.0)
    
    # 3. create a copy of the pose for mutation
    mutant = Pose()
    mutant.assign(pose)
    
    # 4. mutate the desired residue
    # the pack_radius argument of mutate_residue (see below) is redundant
    #    for this application since the area around the mutation is already
    #    repacked
    print "\n\n\n\n\n\n\nAPPLYING MUTATIONS\n\n\n\n\n\n"
    for mutant_position, mutant_aa in muts:
        mutant = mutate_residue(mutant, mutant_position, mutant_aa, 0.0, scorefxn)
    
    print "\n\n\n\n\n\n\nCALCULATING ENERGIES\n\n\n\n\n\n"
    # 5. calculate the "interaction energy"
    # the method calc_interaction_energy is exposed in PyRosetta however it
    #    does not alter the protein conformation after translation and may miss
    #    significant interactions
    # an alternate method for manually separating and scoring is provided called
    #    calc_binding_energy (see Interaction Energy vs. Binding Energy below)
    wt_score, wt_complex, wt_packed_complex, wt_packed_separate = calc_binding_energy(wt, scorefxn, [x[0] for x in muts], cutoff)
    
    mut_score, mut_complex, mut_packed_complex, mut_packed_separate = calc_binding_energy(mutant, scorefxn, [x[0] for x in muts], cutoff)
    #### the method calc_interaction_energy separates an input pose by
    ####    500 Angstroms along the jump defined in a Vector1 of jump numbers
    ####    for movable jumps, a ScoreFunction must also be provided
    #### if setup_foldtree has not been applied, calc_interaction_energy may be
    ####    wrong (since the jumps may be wrong)
    #wt_score = calc_interaction_energy(wt, scorefxn, movable_jumps)
    #mut_score = calc_interaction_energy(mutant, scorefxn, movable_jumps)
    ddg = mut_score - wt_score

    # -write the mutant structure to a PDB file
    print "\n\n\n\n\n\n\nWRITING OUTPUT\n\n\n\n\n\n"
    if out_filename:
        wt_packed_complex.dump_pdb(out_filename.replace(".pdb", "WT_Complex.pdb"))
        wt_packed_separate.dump_pdb(out_filename.replace(".pdb", "WT_Separate.pdb"))
        mut_packed_complex.dump_pdb(out_filename.replace(".pdb", "Mut_Complex.pdb"))
        mut_packed_separate.dump_pdb(out_filename.replace(".pdb", "Mut_Separate.pdb"))
    
    return ddg, wt_score, mut_score, wt_complex, mut_complex
# FUNCTION END

# a different version of mutate_residue is provided in PyRosetta v2.0 and
#    earlier that does not optionally repack nearby residues

# replaces the residue at  <resid>  in  <pose>  with  <new_res>  with repacking
def mutate_residue(pose, mutant_position, mutant_aa,
        pack_radius = 0.0, pack_scorefxn = '' ):
    """
    Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
        and repack any residues within  <pack_radius>  Angstroms of the mutating
        residue's center (nbr_atom) using  <pack_scorefxn>
    note: <mutant_aa>  is the single letter name for the desired ResidueType

    example:
        mutate_residue(pose, 30, A)
    See also:
        Pose
        PackRotamersMover
        MutateResidue
        pose_from_sequence
    """
    #### a MutateResidue Mover exists similar to this except it does not pack
    ####    the area around the mutant residue (no pack_radius feature)
    #mutator = MutateResidue(mutant_position, mutant_aa)
    #mutator.apply(test_pose)

    if pose.is_fullatom() == False:
        IOError( 'mutate_residue only works with fullatom poses' )

    test_pose = Pose()
    test_pose.assign(pose)

    # create a standard scorefxn by default
    if not pack_scorefxn:
        pack_scorefxn = get_fa_scorefxn() #  create_score_function('standard')

    task = standard_packer_task(test_pose)

    # the Vector1 of booleans (a specific object) is needed for specifying the
    #    mutation, this demonstrates another more direct method of setting
    #    PackerTask options for design
    aa_bool = rosetta.utility.vector1_bool()
    # PyRosetta uses several ways of tracking amino acids (ResidueTypes)
    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    # aa_from_oneletter returns the integer representation of an amino acid
    #    from its one letter code
    # convert mutant_aa to its integer representation
    mutant_aa = core.chemical.aa_from_oneletter_code(mutant_aa)

    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design
    # to do this, construct a Vector1 of booleans indicating which amino acid
    #    (by its numerical designation, see above) to allow
    for i in range(1, 21):
        # in Python, logical expression are evaluated with priority, thus the
        #    line below appends to aa_bool the truth (True or False) of the
        #    statement i == mutant_aa
        aa_bool.append( i == int(mutant_aa) )

    # modify the mutating residue's assignment in the PackerTask using the
    #    Vector1 of booleans across the proteogenic amino acids
    task.nonconst_residue_task(mutant_position
        ).restrict_absent_canonical_aas(aa_bool)

    # prevent residues from packing by setting the per-residue "options" of
    #    the PackerTask
    center = pose.residue(mutant_position).nbr_atom_xyz()
    for i in range(1, pose.total_residue() + 1):
        # only pack the mutating residue and any within the pack_radius
        if not i == mutant_position or center.distance_squared(
                test_pose.residue(i).nbr_atom_xyz()) > pack_radius**2:
            task.nonconst_residue_task(i).prevent_repacking()

    # apply the mutation and pack nearby residues
    packer = protocols.minimization_packing.PackRotamersMover(pack_scorefxn, task)
    packer.apply(test_pose)

    return test_pose
# FUNCTION END

# there is a significant difference between interaction energy and binding
#    energy (see Interaction Energy vs. Binding Energy below)
def calc_binding_energy(pose, scorefxn, mut_poses, cutoff = 8.0):
    # create a copy of the pose for manipulation
    test_pose = Pose()
    test_pose.assign(pose)

    # setup packer options
    # the sidechain conformations of residues "near the interface", defined as
    #    within  <cutoff>  Angstroms of an interface residue, may change and
    #    must be repacked, if all residues are repacked, aberrant sidechain
    #    conformations near the interface, but independent of complex
    #    interactions, will be repacked for the mutant and wild-type structures
    #    preventing them from adding noise to the score difference
    # this method of setting up a PackerTask is different from packer_task.py
    tf = standard_task_factory()    # create a TaskFactory
    tf.push_back(core.pack.task.operation.RestrictToRepacking())    # restrict it to repacking

    # this object contains repacking options, instead of turning the residues
    #    "On" or "Off" directly, this will create an object for these options
    #    and assign it to the TaskFactory
    prevent_repacking = core.pack.task.operation.PreventRepacking()

    # the "center" (nbr_atom) of the mutant residue, for distance calculation
    centers = [test_pose.residue(x).nbr_atom_xyz() for x in mut_poses]
    for i in range(1, test_pose.total_residue() + 1):
        # the .distance_squared method is (a little) lighter than .norm
        # if the residue is further than <cutoff> Angstroms away, do not repack
        flag = True
        for center in centers:            
            if center.distance_squared(
                    test_pose.residue(i).nbr_atom_xyz()) <= cutoff**2:
                flag = False
                break
        if(flag):
            prevent_repacking.include_residue(i)
        else:
            print "\n\n\nRESIDUE NEAR INTERFACE\n\n\n"

    # apply these settings to the TaskFactory
    tf.push_back(prevent_repacking)

    # setup a PackRotamersMover
    packer = protocols.minimization_packing.PackRotamersMover(scorefxn)
    packer.task_factory(tf)

    #### create a Mover for performing translation
    #### RigidBodyTransMover is SUPPOSED to translate docking partners of a
    ####    pose based on an axis and magnitude
    #### test it using the PyMOLMover, it does not perform a simple translation
    ####    I also observed a "Hbond Tripped" error when packing after applying
    ####    the Mover, it appears to store inf and NaN values into hbonds
    #transmover = RigidBodyTransMover()
    # calc_interaction_energy separates the chains by 500.0 Angstroms,
    #    so does this Mover
    # if using this Mover, the step_size MUST be a float
    # if this setting is left to default, it will move the proteins
    #    VERY far apart
    #transmover.step_size( 5.0 )

    # repack the test_pose
    packer.apply(test_pose)
    packed_complex = Pose()
    packed_complex.assign(test_pose)

    # score this structure
    before = scorefxn(test_pose)

    # separate the docking partners
    #### since RigidBodyTransMover DOES NOT WORK, it is not used
    #transmover.apply(test_pose)

    # here are two methods for applying a translation onto a pose structure
    # both require an xyzVector
    xyz = rosetta.numeric.xyzVector_double_t()    # a Vector for coordinates
    xyz.x = 500.0    # arbitrary separation magnitude, in the x direction
    xyz.y = 0.0    #...I didn't have this and it defaulted to 1e251...?
    xyz.z = 0.0    #...btw thats like 1e225 light years,
                   #    over 5e245 yrs at Warp Factor 9.999 (thanks M. Pacella)

    #### here is a hacky method for translating the downstream partner of a
    #    pose protein-protein complex (must by two-body!)
    chain2starts = len(pose.chain_sequence(1)) + 1
    for r in range(chain2starts, test_pose.total_residue() + 1):
        for a in range(1, test_pose.residue(r).natoms() + 1):
            test_pose.residue(r).set_xyz(a,
                test_pose.residue(r).xyz(a) + xyz)

    # here is an elegant way to do it, it assumes that jump number 1
    #    defines the docking partners "connectivity"
    # the pose.jump method returns a jump object CREATED from the pose jump
    #    data, the pose itself does not own a Jump object, thus you can use
    #    Jump methods, such as pose.jump(1).set_translation, however the object
    #    has not been properly constructed for manipulation, thus performing
    #    a change does not cause any problems, but is not permanently applied
    #translate = test_pose.jump( 1 )    # copy this information explicitly
    # adjust its translation via vector addition
    #translate.set_translation( translate.get_translation() + xyz )
    #test_pose.set_jump( 1 , translate )
    # as explained above, this call will NOT work
    #test_pose.jump(1).set_translation( test_pose.get_translation() + xyz )

    # repack the test_pose after separation
    packer.apply(test_pose)

    # return the change in score
    return before - scorefxn(test_pose), before, packed_complex, test_pose
# FUNCTION END

import sys
import pandas as pd

try:
    _, p1, p2, pdb, rank, muts, trials, interface_cutoff, out_base = sys.argv
except:
    raise
    _, p1, p2, pdb, rank, muts, trials, interface_cutoff = sys.argv
    out_base = "Data/ddG_Mutated_Structures/"

muts = [(int(x.split("_")[0]), x.split("_")[1]) for x in muts.split(",")]
rank = int(rank)
trials = int(trials)
interface_cutoff = float(interface_cutoff)


# 1. create a pose from the desired PDB file
pose = Pose()
pose_from_file(pose, pdb)

muts = [(pose.pdb_info().pdb2pose("A", x[0]), x[1]) for x in muts]

# 2. setup the docking FoldTree and other related parameters
dock_jump = 1
movable_jumps = Vector1([dock_jump])
protocols.docking.setup_foldtree(pose, "A_B", movable_jumps)

# 3. create ScoreFuncions for the Interface and "ddG" calculations
# the pose's Energies objects MUST be updated for the Interface object to
#    work normally
scorefxn = get_fa_scorefxn() #  create_score_function('standard')
scorefxn(pose)    # needed for proper Interface calculation

# setup a "ddG" ScoreFunction, custom weights
ddG_scorefxn = ScoreFunction()
ddG_scorefxn.set_weight(core.scoring.fa_atr, 0.44)
ddG_scorefxn.set_weight(core.scoring.fa_rep, 0.07)
ddG_scorefxn.set_weight(core.scoring.fa_sol, 1.0)
ddG_scorefxn.set_weight(core.scoring.hbond_bb_sc, 0.5)
ddG_scorefxn.set_weight(core.scoring.hbond_sc, 1.0)

# 6. Mutate Protein and Calculate ddG
summary = []
for trial in range(trials):
    # store the ddG values in a dictionary
    ddG_mutants = {}
    
    filename = "{0}{1}_{2}_{3}_SARS_Muts_{4}.pdb".format(out_base, p1, p2, rank, trial+1)
    ddg, wt_score, mut_score, wt_complex, mut_complex = interface_ddG(pose, muts, movable_jumps, ddG_scorefxn, interface_cutoff, filename)
    
    summary.append([p1, p2, rank, trial, wt_complex, wt_score, mut_complex, mut_score, ddg, filename])
summary = pd.DataFrame(summary, columns=["P1", "P2", "Docking_Rank", "ddG_Trial", "WT_Score", "WT_dG", "Mut_Score", "Mut_dG", "ddG", "pdbfile"])
if(not os.path.exists("{0}Summary_Logs/".format(out_base))):
   os.mkdir("{0}Summary_Logs/".format(out_base))
summary.to_csv("{0}Summary_Logs/{1}_{2}_{3}".format(out_base, p1, p2, rank), sep="\t", index=None)
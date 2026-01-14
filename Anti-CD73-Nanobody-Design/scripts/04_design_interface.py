import pyrosetta
from pyrosetta import rosetta
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    InitializeFromCommandline,
    IncludeCurrent,
    OperateOnResidueSubset,
    RestrictToRepackingRLT,
    PreventRepackingRLT
)
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector,
    NeighborhoodResidueSelector,
    NotResidueSelector,
    AndResidueSelector
)
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover, MinMover

# 1. Init
pyrosetta.init("-mute all")
print("1. Loading docked complex...")
pose = pyrosetta.pose_from_pdb("final_docked.pdb")
scorefxn = pyrosetta.create_score_function("ref2015")

# 2. Define the Interface (Where the action happens)
# We need to tell Rosetta WHERE to mutate.
print("2. Defining Interface for Design...")

# Select Chain A (CD73 - Target) and Chain B (Nanobody)
chainA = ChainSelector("A")
chainB = ChainSelector("B")

# Select residues on the Nanobody (B) that are close (8A) to the Target (A)
# These are the residues we want to DESIGN (Mutate)
interface_on_nanobody = AndResidueSelector(
    chainB, 
    NeighborhoodResidueSelector(chainA, 8.0, False)
)

# Select residues on the Target (A) that are close to the Nanobody
# These we want to REPACK (Wiggle, but don't mutate target!)
interface_on_target = AndResidueSelector(
    chainA, 
    NeighborhoodResidueSelector(chainB, 8.0, False)
)

# Everything else? Freeze it to save time.
not_interface = NotResidueSelector(
    rosetta.core.select.residue_selector.OrResidueSelector(interface_on_nanobody, interface_on_target)
)

# 3. Setup the Packer Task (The "Brain" of the mutation)
tf = TaskFactory()
tf.push_back(InitializeFromCommandline())
tf.push_back(IncludeCurrent()) # Always consider the original amino acid

# RULE 1: Don't touch the non-interface parts
tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), not_interface))

# RULE 2: Only wiggle (repack) the target interface, don't mutate it
tf.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), interface_on_target))

# (Implicit Rule 3: The interface_on_nanobody is left alone, so it defaults to DESIGN/MUTATE)

# 4. Run the Design
print("3. Running Sequence Design (Mutating Nanobody Interface)...")
packer = PackRotamersMover(scorefxn)
packer.task_factory(tf)
packer.apply(pose)

# 5. Minimize (Relax the new mutations)
print("4. Minimizing energy of new sequence...")
minmover = MinMover()
minmover.movemap_factory(pyrosetta.rosetta.core.select.movemap.MoveMapFactory())
minmover.score_function(scorefxn)
minmover.apply(pose)

# 6. Save and Report
final_score = scorefxn(pose)
print("-" * 40)
print(f"Original Score: -1357.66 REU")
print(f"Designed Score: {final_score:.2f} REU")

pose.dump_pdb("designed_nanobody.pdb")
print("Saved 'designed_nanobody.pdb'")

# 7. Print the Mutations (What changed?)
# We compare the sequence of the scaffold vs the new design
# (A simple sequence print for manual checking)
print("\nNew Nanobody Sequence (Chain B):")
print(pose.chain_sequence(2))

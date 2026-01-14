import pyrosetta
from pyrosetta import rosetta

# 1. Initialize with specific options for docking
# -ex1 -ex2aro: Allow sidechains to rotate more freely (helps fix clashes)
pyrosetta.init("-mute all -ex1 -ex2aro")

print("1. Loading Safe Complex...")
pose = pyrosetta.pose_from_pdb("start_complex.pdb")

# 2. Define the 'FoldTree' (Kinematics)
# This tells Rosetta that Chain A (Target) is the anchor, and Chain B (Nanobody) is movable.
# The "1" means they are connected by 'Jump 1'.
rosetta.protocols.docking.setup_foldtree(pose, "A_B", rosetta.utility.vector1_int(1))

# 3. Setup the Smart Docking Protocol
# This object automatically performs:
# - Low-Resolution Search (Centroid) to resolve the 74,000 clash
# - High-Resolution Refinement (Full Atom) to optimize binding
dock = rosetta.protocols.docking.DockingProtocol()
dock.set_movable_jumps(rosetta.utility.vector1_int(1)) 

# 4. Create the Physics Rulebook
scorefxn = pyrosetta.create_score_function("ref2015")

# 5. Run the Simulation
print("2. Running Docking Simulation...")
print("   (This acts like 'Soft Clay' first to fix the overlap, then hardens.)")
print("   Please wait (approx 3-5 minutes)...")

dock.apply(pose)

# 6. Report Results
final_score = scorefxn(pose)
print("-" * 30)
print(f"Final Docked Score: {final_score:.2f} REU")

if final_score < -100:
    print("STATUS: SUCCESS. Excellent binding energy.")
elif final_score < 0:
    print("STATUS: GOOD. Stable complex formed.")
else:
    print("STATUS: POOR. The proteins might just be touching without binding.")

pose.dump_pdb("final_docked.pdb")
print("Saved 'final_docked.pdb'")

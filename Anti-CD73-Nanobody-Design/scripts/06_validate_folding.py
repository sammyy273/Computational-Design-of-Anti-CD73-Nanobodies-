import pyrosetta
from pyrosetta import rosetta
import requests
import time

pyrosetta.init("-mute all")

def solve_rmsd(pose1, pose2):
    # Align them on the backbone atoms
    # RMSD (Root Mean Square Deviation) measures the average distance between atoms.
    # Lower is better. < 2.0 Angstroms = "It folds correctly!"
    return rosetta.core.scoring.CA_rmsd(pose1, pose2)

print("1. Preparing Sequence...")
# Load your design to get the sequence
design_pose = pyrosetta.pose_from_pdb("designed_nanobody.pdb")
# Extract Chain B (The Nanobody)
nanobody_sequence = design_pose.chain_sequence(2)
print(f"   Sequence: {nanobody_sequence[:20]}...")

print("2. Sending to ESMFold AI (Meta)...")
# We use the public API to fold the sequence
url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
response = requests.post(url, data=nanobody_sequence, verify=False)

if response.status_code == 200:
    print("   -> Prediction received!")
    with open("esm_prediction.pdb", "w") as f:
        f.write(response.text)
else:
    print(f"   [ERROR] API request failed: {response.status_code}")
    exit()

print("3. Comparing Design vs. AI Prediction...")
# Load the AI prediction
esm_pose = pyrosetta.pose_from_pdb("esm_prediction.pdb")

# Load JUST the nanobody from your design (Split by chain)
# Your design file has Target + Nanobody. We need just the Nanobody to compare.
split_design = design_pose.split_by_chain()
design_nanobody = split_design[2] # Chain 2 is the nanobody

# Calculate RMSD
rmsd = solve_rmsd(design_nanobody, esm_pose)

print("-" * 40)
print(f"RMSD (Error Metric): {rmsd:.3f} Å")
print("-" * 40)

if rmsd < 2.0:
    print("RESULT: PASS. The AI confirms your sequence folds exactly as designed.")
elif rmsd < 4.0:
    print("RESULT: ACCEPTABLE. The general fold is correct, loops might vary.")
else:
    print("RESULT: FAIL. The sequence folds into a different shape than expected.")

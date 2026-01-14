import pyrosetta
from pyrosetta import rosetta

# 1. Init
pyrosetta.init("-mute all")

def get_sequence(pdb_file, chain_id):
    pose = pyrosetta.pose_from_pdb(pdb_file)
    # Get the sequence for the specific chain (1=Target, 2=Nanobody)
    # In our setup, Nanobody is always the second chain moved into place.
    # But let's be safe and grab the chain by ID if possible, or just grab Chain 1 from scaffold.
    # Since 'scaffold.pdb' is JUST the nanobody, it is Chain 1.
    # In 'designed_nanobody.pdb', the nanobody is Chain B (index 2).
    
    if "scaffold" in pdb_file:
        return pose.chain_sequence(1)
    else:
        return pose.chain_sequence(2)

print("1. Comparing sequences...")

# Load the sequences
seq_original = get_sequence("scaffold.pdb", 1)
seq_design   = get_sequence("designed_nanobody.pdb", 2)

print(f"Original Length: {len(seq_original)}")
print(f"Designed Length: {len(seq_design)}")

if len(seq_original) != len(seq_design):
    print("[WARNING] Lengths differ! Alignment might be tricky.")
    # Simple truncation to match lengths for comparison
    min_len = min(len(seq_original), len(seq_design))
    seq_original = seq_original[:min_len]
    seq_design = seq_design[:min_len]

# Compare residue by residue
mutations = []
print("\n--- MUTATION REPORT ---")
print(f"{'Residue':<10} {'Original':<10} {'Designed':<10} {'Type'}")
print("-" * 45)

for i in range(len(seq_original)):
    aa_orig = seq_original[i]
    aa_new  = seq_design[i]
    
    if aa_orig != aa_new:
        residue_num = i + 1
        # Classify mutation
        mut_type = "Surface" # Default
        # CDR3 is roughly 95-105
        if 95 <= residue_num <= 105:
            mut_type = "CDR3 (Binding Loop)"
        
        mutations.append(f"{aa_orig}{residue_num}{aa_new}")
        print(f"{residue_num:<10} {aa_orig:<10} {aa_new:<10} {mut_type}")

print("-" * 45)
print(f"Total Mutations: {len(mutations)}")
print(f"Mutation String: {', '.join(mutations)}")

# Write to file for your report
with open("mutation_report.txt", "w") as f:
    f.write("Original: " + seq_original + "\n")
    f.write("Designed: " + seq_design + "\n")
    f.write("Mutations: " + ", ".join(mutations) + "\n")
    
print("\nSaved details to 'mutation_report.txt'.")

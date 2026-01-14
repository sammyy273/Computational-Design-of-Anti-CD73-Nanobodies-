import pyrosetta
from pyrosetta import rosetta

# Initialize PyRosetta (mute noisy output)
pyrosetta.init("-mute all")

def extract_chain_clean(input_pdb, chain_letter, output_file):
    print(f"Reading {input_pdb}...")
    
    # 1. Load the full PDB
    # This loads everything: Protein, DNA, Water, Ions
    full_pose = pyrosetta.pose_from_pdb(input_pdb)
    
    # 2. Split the complex into a list of separate poses
    # This automatically separates chains and usually separates waters/ligands
    split_poses = full_pose.split_by_chain()
    
    target_pose = None
    
    # 3. Find the specific chain we want
    for i in range(1, len(split_poses) + 1):
        p = split_poses[i]
        
        # Get the PDB Chain ID (e.g., 'A', 'B') of the first residue
        if p.total_residue() > 0:
            current_chain = p.pdb_info().chain(1)
            
            if current_chain == chain_letter:
                target_pose = p
                break
    
    # 4. Save it if found
    if target_pose:
        print(f"   -> Found Chain {chain_letter}. Residues: {target_pose.total_residue()}")
        target_pose.dump_pdb(output_file)
        print(f"   -> Saved to {output_file}")
    else:
        print(f"   [ERROR] Chain {chain_letter} not found in {input_pdb}!")


print("-" * 40)
# 1. Get CD73 (Chain A) from 4H1S
extract_chain_clean("4H1S.pdb", "A", "target.pdb")

print("-" * 40)
# 2. Get Nanobody (Chain B) from 5F1O
extract_chain_clean("5F1O.pdb", "B", "scaffold.pdb")

print("-" * 40)
print("Done. Check your folder for target.pdb and scaffold.pdb")

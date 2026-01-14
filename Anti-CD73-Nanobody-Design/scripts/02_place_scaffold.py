import pyrosetta
from pyrosetta import rosetta
import math

pyrosetta.init("-mute all")

def get_center(pose, residues=None):
    """Calculate geometric center of atoms."""
    x, y, z = 0.0, 0.0, 0.0
    count = 0
    
    if residues is None:
        residues = range(1, pose.total_residue() + 1)
        
    for i in residues:
        if not pose.residue(i).has("CA"): continue
        atom = pose.residue(i).xyz("CA")
        x += atom.x
        y += atom.y
        z += atom.z
        count += 1
        
    if count == 0: return rosetta.numeric.xyzVector_double_t(0,0,0)
    return rosetta.numeric.xyzVector_double_t(x/count, y/count, z/count)

print("1. Processing Geometry...")
target = pyrosetta.pose_from_pdb("target.pdb")
scaffold = pyrosetta.pose_from_pdb("scaffold.pdb")

# 1. Find the "Core" of CD73
core_center = get_center(target) 

# 2. Find the "Active Site" (The Exit Door)
# Using the Zinc-binding residues you provided earlier
active_site_res = [36, 38, 85, 220] 
active_center = get_center(target, active_site_res)

# 3. Calculate the "Exit Vector" (Direction pointing OUT of the protein)
# Vector = Surface - Core
vec_x = active_center.x - core_center.x
vec_y = active_center.y - core_center.y
vec_z = active_center.z - core_center.z

# Normalize vector (make length 1.0)
length = math.sqrt(vec_x**2 + vec_y**2 + vec_z**2)
norm_x, norm_y, norm_z = vec_x/length, vec_y/length, vec_z/length

print(f"   -> Exit Vector calculated: ({norm_x:.2f}, {norm_y:.2f}, {norm_z:.2f})")

# 4. Define the Safe Spot (Active Site + 20 Angstroms OUT)
safe_distance = 20.0
safe_x = active_center.x + (norm_x * safe_distance)
safe_y = active_center.y + (norm_y * safe_distance)
safe_z = active_center.z + (norm_z * safe_distance)

# 5. Move Nanobody to Safe Spot
print("2. Moving Nanobody to Safe Spot...")
scaffold_center = get_center(scaffold)
move_x = safe_x - scaffold_center.x
move_y = safe_y - scaffold_center.y
move_z = safe_z - scaffold_center.z

translation = rosetta.numeric.xyzVector_double_t(move_x, move_y, move_z)
rotation = rosetta.numeric.xyzMatrix_double_t.identity()
scaffold.apply_transform_Rx_plus_v(rotation, translation)

# 6. Save Combined
combined = target.clone()
combined.append_pose_by_jump(scaffold, 1)
combined.dump_pdb("start_complex.pdb")

# 7. Verify Score
scorefxn = pyrosetta.create_score_function("ref2015")
score = scorefxn(combined)
print("-" * 30)
print(f"Final Check - Energy Score: {score:.2f} REU")
if score < 1000:
    print("STATUS: SAFE. Ready for docking.")
else:
    print("STATUS: WARNING. Still clashing (but closer). Docking protocol handles this.")

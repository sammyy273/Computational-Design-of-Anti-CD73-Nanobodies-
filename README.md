# Computational Design of Anti-CD73 Nanobodies
_A computational pipeline to design high-affinity nanobody antagonists against the cancer target CD73._

Solid tumors use the enzyme CD73 to suppress the immune system. Conventional antibody drugs are often too large to penetrate dense tumor tissue. This project utilizes PyRosetta (Rosetta Macromolecular Modeling Suite) to computationally design a Single-Domain Antibody (Nanobody) that is 1/10th the size of a standard antibody and binds to the CD73 active site.

# Key Results
| **Metric**     |   **Result**    |                  **Meaning**                                     |
| -------------- | --------------- | -----------------------------------------------------------------|
| Binding Energy | -1362 REU       | Strong, stable complex formation.                                |
| Mutations      | S7D, Q13H, S21L | Specific surface optimizations for solubility & packing.         |
| Validation     | 1.97 Å RMSD     | AI (ESMFold) confirms the design folds correctly (98% accuracy). |

# Structural Validation
Comparison of the Physics-based design (Cyan) vs. the Deep Learning prediction (Salmon).

_(Note: The low RMSD indicates high structural confidence in the designed sequence.)_

# Methodology & Pipeline

**_The project follows a reproducible 6-step computational workflow:_**

**Structure Prep:** Isolation of Human CD73 (PDB: 4H1S) and a Llama VHH scaffold (PDB: 5F1O).

**Geometric Alignment:** Vector-based positioning of the nanobody 20Å from the active site to prevent steric clashes.

**Docking:** Monte Carlo-based docking simulation (DockingProtocol) using Centroid and Full-Atom scoring.

**Affinity Maturation:** Sequence design (PackRotamersMover) to optimize the binding interface.

**Analysis:** Mutation profiling and energy scoring (Ref2015).

**Validation:** Orthogonal validation using ESMFold (Meta AI) to predict the 3D structure of the de novo sequence.

# Installation & Usage
# _Prerequisites_
​Python 3.10 (Required for PyRosetta compatibility)

​PyRosetta

# Setup
**1. Clone the repository**
git clone [https://github.com/sammyy273/Anti-CD73-Nanobody-Design.git](https://github.com/sammyy273/Anti-CD73-Nanobody-Design.git)

cd Anti-CD73-Nanobody-Design

**OR**

**2. Install it:**

_Download the correct file (For 64-bit; **Run Windows Power shell as adminstrator**):_

wget https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.MinSizeRel.python310.ubuntu.wheel/pyrosetta-2025.45+release.d79cb06334-cp310-cp310-linux_x86_64.whl

**_Install it:_**

pip install pyrosetta-2025.45+release.d79cb06334-cp310-cp310-linux_x86_64.whl

**_Test it using:_**

python -c "import pyrosetta; pyrosetta.init()"

import pyrosetta

from pyrosetta.toolbox import pose_from_rcsb

**_Initialize PyRosetta_**

pyrosetta.init()

# Run the Pipeline
**_​The scripts are numbered sequentially for reproducibility:_**

python scripts/01_prepare_pdbs.py      # Clean raw PDBs

python scripts/02_place_scaffold.py    # Align Nanobody to Active Site

python scripts/03_run_docking.py       # Run Docking Simulation

python scripts/04_design_interface.py  # Mutate/Design Sequence

python scripts/06_validate_folding.py  # Check structure with AI

# License
​This project is for academic/educational purposes. PyRosetta requires a separate license.



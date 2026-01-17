QE-KIT (v0.1.0)

An Interactive Workflow Automation Tool for Quantum ESPRESSO

QE-KIT is a Python-based ecosystem designed to bridge the gap between structural modeling and high-level physics properties. It automates the tedious aspects of Quantum ESPRESSO (QE) simulations, from CIF file conversion and "Smart" pseudopotential mapping to structural refinement and phonon post-processing.
🛠️ Installation & Setup
1. Prerequisites

    Python 3.8+

    Quantum ESPRESSO Utilities: Ensure cell2ibrav.x and cif2cell are installed and available in your system $PATH.

    QE-POTCAR Library: A structured directory of .UPF files (see Library Structure below).

2. Deployment

Clone the repository and run the automated setup script. This will create a virtual environment, install dependencies, and register the qekit command globally.
Bash

git clone https://github.com/yourusername/QE-KIT.git
cd QE-KIT
chmod +x setup_env.sh
./setup_env.sh

During setup, you will be prompted to enter the absolute path to your QE-POTCAR library. This path persists across sessions.
📂 QE-POTCAR Library Structure

The tool utilizes a "Smart Scanner" to traverse your pseudopotential library. To function correctly, your library must follow this hierarchy:
Plaintext

QE-POTCAR/
├── LDA/
│   ├── full/          # Full Relativistic
│   │   ├── PAW/
│   │   └── USPP/      # Ultrasoft
│   └── scalar/        # Scalar Relativistic
│       ├── PAW/
│       └── USPP/
├── PBE/
│   ├── full/
│   └── scalar/
└── NC-mt/             # Norm-Conserving (Required for Phonon/Optical)

🚀 The 3-Zone Workflow
Zone 1: Pre-Processing

    Module 100 (Structure2In): Converts .cif files to a clean scf.dat template.

    Module 101 (Pseudo Select): Automatically scans your library and maps the correct .UPF files to your elements based on the selected functional (PBE, LDA, or NC).

    Module 102 (K-Point Gen): Generates uniform grids for SCF or high-symmetry paths for Band/Phonon dispersions.

Zone 2: Input Generation

Surgically extracts geometry and species metadata from scf.dat to generate ready-to-run input files:

    200-204: SCF, NSCF, Bands, and Optical inputs.

    201 (VC-Relax): Structural optimization input.

    205 (Phonon ph.x): Automated setup for DFPT calculations, including mass and dynamical matrix naming.

Zone 3: Post-Processing & Refinement

    Module 300 (Structural Refinement): 1. Parses relax.out for converged coordinates. 2. Runs cell2ibrav.x to find the exact symmetry parameters (ibrav and celldm). 3. Updates scf.in surgically, removing the old CELL_PARAMETERS and ATOMIC_SPECIES to prevent QE errors.

    Module 301 (Phonon Processing): 1. Scans ph.in for metadata. 2. Generates q2r.in, matdyn.in (for DOS), and matdyndisp.in (for Bands). 3. Automatically appends the K-path from your scf.dat.

📝 Example Usage

    Open a terminal in your project folder.

    Type qekit.

    Select 100 to load your CIF.

    Select 101 to pick your library flavor.

    Generate your vc-relax.in and run it.

    Once finished, run qekit and select 300 to refine your structure into a standard ibrav format for property runs.

Farhan Noor - Developer & Maintainer

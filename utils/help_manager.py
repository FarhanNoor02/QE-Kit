import sys

def display_help():
    """Prints a professional, categorized help menu to the terminal."""
    help_text = r"""
QE-Kit: A Pre- & Post-Processing Suite for Quantum ESPRESSO
===========================================================
Version: 0.3.4 | Lead Dev: Farhan Noor (Univ. of Dhaka)

USAGE:
    qekit           : Launch interactive menu (Default)
    qekit --help    : Show this help documentation
    qekit --[code]  : Run a specific module in headless mode

ZONE 1: STRUCTURAL DISCOVERY
    --100 : Convert .cif files to PWscf input templates.
    --101 : Interactive Pseudopotential selection. Selects from the QE-POTCAR 
            library, which must be configured during initial setup.
    --102 : Generate K-point grid. Available options: 
            (A) Uniform (Automatic grid)
            (B) Homogeneous (Explicit crystal grid; e.g., for optical spectra)
            (C) Band (High-symmetry path)

ZONE 2: INPUT ARCHITECT
    --200 : Generate SCF input files. Users will specify the ecutwfc and ecutrho 
            thresholds, read directly from the pseudopotential files.
    --206 : Generate Phonon (ph.in) input files (Automated).
    --207: Zone 2: Generate fs.in for Fermi surface extraction from NSCF output.
    --208: Zone 2: Generate thermo_control for elastic/thermodynamic properties.

ZONE 3: POST-PROCESSOR
    --300 : Symmetry refinement (Updates scf.in with optimized lattice).
    --301 : Phonon post-processing (Automated). Generates files for use with 
            q2r.x, matdyn.x, and plotband.x.
    --303 : Isotropic optical constant derivation (n, k, R, alpha).
    --305: Thermo_PW result extraction

ZONE 4: DATA VISUALIZATION
    --401 : Run HT-Phonon Pipeline (Bash Driver).
    --404 : Interactive Brillouin Zone Visualizer (Matplotlib).
    --405 : Calculated XRD Pattern Plotter (Pymatgen).
    --406 : 3D Fermi Surface Visualizer (scikit-image)

DOCUMENTATION & SUPPORT:
    GitHub: https://github.com/FarhanNoor02/QE-Kit
    """
    print(help_text)
    sys.exit(0)

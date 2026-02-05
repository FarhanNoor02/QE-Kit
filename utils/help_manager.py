import sys

def display_help():
    """Prints a professional, categorized help menu to the terminal."""
    help_text = """
QE-Kit: A Pre- & Post-Processing Suite for Quantum ESPRESSO
===========================================================
Version: 0.3.3 | Lead Dev: Farhan Noor

USAGE:
    qekit           : Launch interactive menu (Default)
    qekit --help    : Show this help documentation
    qekit --[code]  : Run a specific module in headless mode

ZONE 1: STRUCTURAL DISCOVERY
    --100 : Convert .cif files to PWscf input templates.
    --101 : Interactive Pseudopotential selection. Selects from the QE-POTCAR library which is to
            be installed during initial setup of qekit
    --102 : Generate K-point grid. 
            Available options: A) Uniform (Automatic grid),(B) Homogeneous (Explicit crystal grid;
            list all Kpoints- for example: used for optical spectra), (C) Band (High-symmetry path)

ZONE 2: INPUT ARCHITECT
    --200 : Generate SCF input files. Users will specify the ecut and ewfc thresholds; to be read
            from the pseudopotential files
    --206 : Generate Phonon (ph.in) input files (Automated).

ZONE 3: PHYSICS & ANALYSIS
    --300 : Symmetry refinement (Updates scf.in with optimized lattice).
    --301 : Phonon post-processing (Automated). Generates files for use with q2r.x, matdyn.x and
            plotband.x
    --303 : Isotropic optical constant derivation (n, k, R, alpha).

DOCUMENTATION & SUPPORT:
    GitHub: https://github.com/FarhanNoor02/QE-Kit
    """
    print(help_text)
    sys.exit(0)

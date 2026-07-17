# utils/citations.py
import os

def display_citations():
    # Clear the screen for a clean reading experience
    os.system('clear' if os.name == 'posix' else 'cls')
    
    citations_text = """
======================================================================
                 QE-KIT PACKAGE DEPENDENCIES & CITATIONS
======================================================================
If you use QE-Kit in your research, please ensure you cite the underlying 
engines and libraries that power these calculations.

[1] Pymatgen (Python Materials Genomics)
    Ong, S. P., et al.(2013).Python Materials Genomics (pymatgen): A Robust, Open-Source Python Library for
    Materials Analysis. Computational Materials Science, 68, 314-319. doi:10.1016/j.commatsci.2012.10.028

[2] ASE (Atomic Simulation Environment)
    Larsen, A. H., et al. (2017). The atomic simulation environment-a Python library for working with atoms.
    J Phys Condens Matter. 2017 Jul 12;29(27):273002. doi: 10.1088/1361-648X/aa680e.

[3] Spglib (Symmetry analysis & cell refinement)
   A. Togo, I. Tanaka, "Spglib: a software library for crystal symmetry search", arXiv:1808.01590 (2018)

[4] Seekpath (High-symmetry K-path generation)
    Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka.
    Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017) .

[5] cif2cell (CIF parser)
    Torbjorn Bjorkman, "CIF2Cell: Generating geometries for electronic structure programs",
    Computer Physics Communications 182, 1183-1186 (2011) doi: 10.1016/j.cpc.2011.01.013

======================================================================
    """
    
    print(citations_text)
    input("\nPress [ENTER] to return to the Main Menu...")

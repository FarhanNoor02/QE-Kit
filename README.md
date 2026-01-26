# QE-Kit (v0.2.2) 
### A Pre- & Post-Processing Suite for Quantum ESPRESSO

**Author:** Farhan Noor  
**Affiliation:** Department of Physics, University of Dhaka  
**Status:** Active Development (v0.2.2)

QE-Kit is a Python-powered Command Line Interface (CLI) designed to automate the end-to-end research workflow for Density Functional Theory (DFT) calculations. Developed from over two years of experience optimizing Quantum ESPRESSO workflows, QE-Kit bridges the gap between initial structural discovery and publication-quality results.
---

## 🚀 Key Features

### Zone 1: Pre-Processing & Discovery
* **Structure Transformation:** Intelligent CIF/XYZ to QE Input conversion using ASE and Spglib.
* **Symmetry Analysis:** Automated detection of space groups and high-symmetry K-paths.
* **QE-PotLib Integration:** A managed library system for optimized pseudopotential selection.

### Zone 2: Input Generation
* **Streamlined Setup:** Rapid generation of SCF, NSCF, VC-Relax, and Phonon input files.
* **Modular Logic:** Dedicated modules for Electronic (Bands/DOS) and Optoelectronic (Epsilon) calculation parameters.

### Zone 3: Post-Processing & Analysis
* **Symmetry Refinement:** Automated `cell2ibrav` conversion for optimized structures.
* **Orbital Summation:** Advanced PDOS processing with automated orbital-wise summation.
* **Optoelectronics Suite:** Calculation of isotropic averaged optical constants ($n, k, R, \alpha$) from raw dielectric functions.

---

## 🛠 Installation

QE-Kit is best installed as a global tool using `pipx`:

```bash
pipx install git+https://github.com/FarhanNoor02/QE-Kit.git

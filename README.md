# QE-Kit (v0.3.4) 
### A Pre- & Post-Processing Suite for Quantum ESPRESSO

**Author:** Farhan Noor  
**Affiliation:** Department of Physics, University of Dhaka  
**Status:** Active Development (v0.3.4)

QE-Kit is a Python-powered Command Line Interface (CLI) designed to automate the end-to-end research workflow for Density Functional Theory (DFT) calculations. Developed from over two years of experience optimizing Quantum ESPRESSO workflows, QE-Kit bridges the gap between initial structural discovery and publication-quality results.

---

## 🚀 Architectural Overview

QE-Kit is organized into four distinct **Zones**, representing the logical progression of a materials science research project.



### 🌐 Zone 1: Structural Discovery
Focuses on the transition from experimental crystallography to computational models.
* **Module 100 (Structure to Input):** Robust conversion of `.cif` or `.xyz` files into standard `pw.x` input templates using `ASE`. 
* **Module 101 (Pseudo-Selector):** Interactive selection from your local pseudopotential library. Includes a **Kinetic Energy Advisor** that parses potential files to suggest optimized cutoffs.
* **Module 102 (K-Path Generation):** Utilizes `SeeK-path` to automatically detect high-symmetry points in the Brillouin zone for standardized band structure plots.

### 🏗️ Zone 2: Input File Architect
Rapidly generates specialized input files with optimized defaults for the Quantum ESPRESSO suite.
* **200 - 202 (Electronic Core):** Streamlined setup for SCF, VC-Relax, and NSCF calculations.
* **203 - 204 (Property Mapping):** Templates for Density of States (DOS/PDOS) and Band Structure.
* **205 (Optical Suite):** Specific setup for `epsilon.x` to study complex dielectric responses.
* **206 (Phonon Suite):** Automated generation of `ph.in` for lattice dynamics and vibrational analysis.

### 🧪 Zone 3: Physics & Analysis
The "Brain" of the suite, where raw data is converted into physical insights.
* **Module 300 (Refine Symmetry):** A critical post-relaxation tool. It parses `vc-relax.out`, refines the cell using `Spglib`, and updates the `scf.in` with the exact `ibrav` and atomic coordinates.
* **Module 301 (Phonon Dispersion):** Automates the post-processing of dynamical matrices into a readable dispersion format.
* **Module 302 (PDOS Summation):** Sums orbital-wise contributions (e.g., total $d$-orbital contribution for a transition metal) from multiple `.pdos_atm` files.
* **Module 303 (Optical Constant Derivation):** Calculates isotropic averaged refractive index ($n$), extinction coefficient ($k$), reflectivity ($R$), and absorption ($\alpha$) from raw dielectric function tensors.



### ⚡ Zone 4: Workflow & Visualization
The automation and plotting layer for high-throughput research.
* **HT-Phonon Driver:** A supervised Bash pipeline (`ht_phonon.sh`) that orchestrates the entire phonon workflow—from structural input to final dispersion—without manual intervention.
* **Gnuplot Visualization Suite:**
    * `plot-overlay.sh`: Overlays different data sets (e.g., Spin-up vs Spin-down or Total DOS vs PDOS).
    * `plot-subplots.sh`: Multi-pane vertical stacking for optical spectra to ensure scale-independent analysis across different physical constants.

---

## 🛠 Installation

QE-Kit is best installed as a global tool using `pipx` to maintain an isolated environment:

```bash
# Install directly from the repository
pipx install git+[https://github.com/FarhanNoor02/QE-Kit.git](https://github.com/FarhanNoor02/QE-Kit.git)
```
## System Requirements

    1. Python: 3.8 or higher.

    2. Quantum ESPRESSO: Core binaries (pw.x, ph.x, epsilon.x, etc.) should be in your PATH.

    3. QE-Potcar library: https://github.com/FarhanNoor02/QEPotcar; the pseudopotential library is created via web scraping and 
       it colects and organises all the pseudopotentials used in QE. Path to this library is to be given during installation.

    4. Gnuplot: Required for the visualization suite in Zone 4.

## 📖 Usage
Interactive Mode

Simply type the command to launch the stylish CLI menu:

``` bash

qekit
```
### Headless/Automation Mode
For use in HPC scripts or remote clusters, use command flags:
Bash

qekit --help    # View scientific documentation and command list
qekit --300     # Run automated symmetry refinement
qekit --206     # Generate phonon inputs automatically

## 📈 Roadmap

    v0.4.0: 

    1. Upcoming: DFT+U (Hubbard) support with linear response (hp.x) automation.

    2. Upcoming: Spin-Orbit Coupling (SOC) and Wannier90 integration.

    3. Upcoming: Electron-Phonon suite for EPW code connectivity.

## To upgrade to the latest version (v0.3.4)
```bash
pipx upgrade qekit
```
## CITATION:
Farhan Noor, F. QE-Kit [Computer software]. https://github.com/FarhanNoor02/QE-Kit

# QE-Kit <img width="75" height="75" alt="uploaded-img" src="https://github.com/user-attachments/assets/08e97847-3685-413d-b630-2cc56c243a61" />

**A Pre- & Post-Processing Suite for Quantum ESPRESSO**
*Developed by: Farhan Noor, University of Dhaka*

QE-Kit is a comprehensive, modular Python and Bash toolkit designed to automate the workflow of Density Functional Theory (DFT) calculations using Quantum ESPRESSO. It bridges the gap between raw crystallographic data and publication-ready physical properties, featuring an interactive CLI, smart pseudopotential management, and robust XML-based calculation state verification.

---

## 🚀 Features

- **Interactive & Headless Modes:** Navigate via a beautiful, interactive CLI using `questionary`, or run modules headlessly in HPC batch scripts (e.g., `./qekit.py --305`).
- **State-Aware Generation:** Modules parse QE XML schemas to ensure prerequisite calculations were successful before generating downstream inputs.
- **Smart Pseudopotential Linking:** Automatically detects local `.UPF` files or fetches them from your global library based on strict element mapping.
- **Advanced Workflows:** Built-in support for standard DFT, Electron-Phonon Coupling (EPC) for superconductivity, and thermodynamic/elastic property extraction via `thermo_pw`.

---

## 📂 Architecture & Modules

The suite is divided into four primary "Zones" and a foundational "Utils" layer.

### 🟢 ZONE 1 | Structural Discovery (Pre-Processing)
*Focuses on preparing the initial crystalline structure and matching it with correct pseudopotentials.*

*   **`100: structure2in`**: Converts raw crystallographic data (like `.cif` files) into a foundational `scf.dat` template. Defines the `ibrav`, `CELL_PARAMETERS`, `ATOMIC_SPECIES`, and `ATOMIC_POSITIONS` blocks.
*   **`101: pseudo_select`**: A "smart" pseudopotential linker. Reads `scf.dat` and checks if required `.UPF` files exist locally. If found, it natively assigns `pseudo_dir='./'`. If missing, it prompts the user to select functional/relativistic options from the global `QE-POTCAR` library, rigorously maps filenames to elements, and injects the absolute path.
*   **`102: kpath_gen`**: Generates high-symmetry K-point paths for the Brillouin zone, necessary for band structure and phonon dispersion calculations.

### 🔵 ZONE 2 | Input File Architect
*Focuses on reading the `scf.dat` template and generating specific, execution-ready `.in` files for Quantum ESPRESSO binaries.*

*   **`200: scf_gen`**: Generates the ground-state Self-Consistent Field (`scf.in`) input file.
*   **`201: vcrelax_gen`**: Generates the variable-cell relaxation input (`vc-relax.in`) for structural optimization.
*   **`202: nscf_gen`**: Generates the Non-Self-Consistent Field input (`nscf.in`), setting up dense uniform grids required for DOS or Fermi surface runs.
*   **`203: pdos_gen`**: Generates the `projwfc.in` file to calculate the Projected Density of States.
*   **`204: bands_gen`**: Generates the `bands.in` file to calculate electronic band dispersion along the k-path.
*   **`205: optical_gen`**: Generates inputs for optical properties (e.g., for `epsilon.x`).
*   **`206: phonon_gen`**: Generates the input file for Density Functional Perturbation Theory (`ph.x`) calculations.
*   **`207: fs_gen`**: Generates `fs.in` for the `fs.x` executable. Utilizes `check_calc_status` to ensure the parent NSCF calculation completed successfully before writing the file.
*   **`208: thermopw_gen`**: Configures the `thermo_control` file for the `thermo_pw.x` driver. Verifies the prior SCF run status and captures user input for external pressure.

### 🟣 ZONE 3 | Post-Processor
*Focuses on extracting, calculating, and formatting raw output data into physical properties.*

*   **`300: optimized`**: Refines symmetry (using `cell2ibrav`) after a `vc-relax` run, updating lattice parameters and atomic positions for subsequent ground-state runs.
*   **`301: phonon_proc`**: Processes raw phonon output data into plottable dispersion curves.
*   **`302: pdos_proc`**: Sums and processes Projected Density of States data for plotting.
*   **`303: optical_proc`**: Derives optical constants (absorption, reflectivity, refractive index) from QE outputs.
*   **`304: epc_processor`**: Analyzes Electron-Phonon Coupling (EPC) data to calculate superconducting parameters like $T_c$.
*   **`305: thermopw_proc`**: Scans `scf.out` from a `thermo_pw` run. Extracts Voigt-Reuss-Hill averages, Debye temperature, applied strain tensors, and the symmetrized elastic stiffness matrix (in GPa), formatting everything into a clean `result.txt` file.

### 🟠 ZONE 4 | Data Visualization
*Focuses on rendering graphical representations of the processed data.*

*   **`401: ht_phonon.sh`**: Bash driver for high-throughput automation of phonon workflows.
*   **`402: plot-overlay.sh`**: Gnuplot utility for overlaying multiple datasets (e.g., comparing band structures).
*   **`403: plot-subplots.sh`**: Gnuplot utility for generating side-by-side subplot panels (e.g., Bands + DOS).
*   **`404: bz_plotter`**: Uses Matplotlib to render an interactive 3D visualization of the Brillouin Zone.
*   **`405: xrd_plotter`**: Uses Pymatgen to calculate and plot the theoretical X-Ray Diffraction (XRD) pattern based on the relaxed crystal structure.

### 🛠️ The `utils/` Library (Backend Engines)
*These modules operate invisibly in the background, providing core logic to the main Zones.*

*   **`config_manager`**: Manages user-level settings (like saving/fetching the global `PSEUDO_LIB_PATH`).
*   **`help_manager`**: Manages the CLI documentation and headless argument flags.
*   **`qe_xml_parser`**: A robust diagnostic tool that parses Quantum ESPRESSO's `data-file-schema.xml`. Powers the calculation verification to prevent downstream modules from failing blindly.
*   **`check_strain`**: Scans `thermo_pw` logs to extract the precise $3 \times 3$ strain tensors applied to the crystal during elastic calculations.
*   **`symmetric_mat`**: Reads raw elastic compliances, forces the $6 \times 6$ matrix into strict symmetry ($C_{ij} = C_{ji}$), converts units from kbar to GPa, and formats it for output.

---

## 💻 Usage

**Interactive Mode:**
Simply launch the orchestrator to bring up the menu interface:
```bash
   qekit.py
```
## 📈 Roadmap

    v0.4.0: 

    1. Upcoming: DFT+U (Hubbard) support with linear response (hp.x) automation.

    2. Upcoming: Spin-Orbit Coupling (SOC) and Wannier90 integration.

    3. Upcoming: Electron-Phonon suite for EPW code connectivity.

## To upgrade to the latest version
```bash
pipx upgrade qekit
```
## CITATION:
Farhan Noor, F. QE-Kit [Computer software]. https://github.com/FarhanNoor02/QE-Kit

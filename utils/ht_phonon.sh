#!/bin/bash

# QE-Kit High-Throughput Phonon Driver (v0.3.2)
# Fully Automated Refinement and Post-Processing logic

echo "========================================================="
echo "   QE-KIT: HIGH-THROUGHPUT PHONON AUTOMATION SUITE      "
echo "========================================================="

# 1. Detection
CIF_FILE=$(ls *.cif | head -n 1)
if [ -z "$CIF_FILE" ]; then
    echo "[!] Error: No .cif file found in this directory."
    exit 1
fi
echo "[+] Detected structure: $CIF_FILE"

# 2. Initialization
read -p "[?] Enter prefix for calculations: " PREFIX
read -p "[?] Enter number of processors (J): " J

# 3. Zone 1: Pre-Processing (Interactive)
# We keep this interactive to ensure pseudopotentials and K-points are correct.
echo -e "\n--- Step 1: Running QE-Kit Zone 1 (Structural Setup) ---"
qekit

# 4. VC-RELAX Execution
if [ -f "relax.in" ]; then
    echo "[+] Starting VC-RELAX on $J cores..."
    mpirun -np $J pw.x < relax.in > relax.out
    grep -q "JOB DONE" relax.out && echo "[✔] VC-RELAX complete." || { echo "[!] VC-RELAX failed."; exit 1; }
else
    echo "[!] Error: relax.in not found."
    exit 1
fi

# 5. Symmetry Refinement (NOW AUTOMATED)
echo -e "\n--- Step 2: Automated Cell Refinement (Module 300) ---"
qekit --300

# 6. SCF Execution
if [ -f "scf.in" ]; then
    echo "[+] Starting SCF calculation..."
    mpirun -np $J pw.x < scf.in > scf.out
    grep -q "JOB DONE" scf.out && echo "[✔] SCF complete." || { echo "[!] SCF failed."; exit 1; }
else
    echo "[!] Error: scf.in not found."
    exit 1
fi

# 7. Phonon Input Generation (NOW AUTOMATED)
echo -e "\n--- Step 3: Automated ph.in Generation (Module 206) ---"
qekit --206

# 8. Phonon Execution
if [ -f "ph.in" ]; then
    echo "[+] Starting Phonon calculation (This may take time)..."
    mpirun -np $J ph.x < ph.in > ph.out
    grep -q "JOB DONE" ph.out && echo "[✔] PH complete." || { echo "[!] PH failed."; exit 1; }
else
    echo "[!] Error: ph.in not found."
    exit 1
fi

# 9. Post-Processing (NOW AUTOMATED)
echo -e "\n--- Step 4: Automated Phonon Post-Processing (Module 301) ---"
qekit --301

echo "========================================================="
echo "   PHONON WORKFLOW COMPLETE FOR $PREFIX"
echo "========================================================="

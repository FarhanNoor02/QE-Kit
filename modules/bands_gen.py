# modules/bands_gen.py
import os
import questionary
from utils import qe_xml_parser

def run_bands_gen():
    dat_file = "scf.dat"
    
    print("\n--- [204] Bands Input Generation ---")
    
    if not os.path.exists(dat_file):
        print(f"[!] Error: {dat_file} not found. Please complete Zone 1 first.")
        return

    prefix = questionary.text("Enter calculation prefix:", default="ncn").ask()
    ecutwfc = questionary.text("ecutwfc (Ry):", default="50").ask()
    ecutrho = questionary.text("ecutrho (Ry):", default="400").ask()

    with open(dat_file, 'r') as f:
        lines = f.readlines()

    # 1. Extract Pseudo Dir
    pseudo_dir = "./"
    for line in lines:
        if line.startswith("# PSEUDO_DIR_PATH:"):
            pseudo_dir = line.split(":")[1].strip()

    # 2. Extract &SYSTEM variables (filtering out the ones we override)
    system_vars = ""
    in_scf_system = False
    for line in lines:
        if "&SYSTEM" in line.upper(): 
            in_scf_system = True
            continue
        if line.strip() == "/" and in_scf_system: 
            in_scf_system = False
            break # Stop after the first &SYSTEM block
        if in_scf_system:
            if not any(x in line for x in ["ecutwfc", "ecutrho", "occupations", "smearing", "degauss"]):
                system_vars += line

    # 3. SURGICAL EXTRACTION OF CARDS ONLY
    # We only want the specific blocks, ignoring redundant &SYSTEM blocks or file paths
    relevant_cards = ""
    target_cards = ["ATOMIC_SPECIES", "ATOMIC_POSITIONS", "CELL_PARAMETERS", "K_POINTS"]
    capture = False
    
    for line in lines:
        # Check if the line starts one of our target blocks
        strip_line = line.strip().upper()
        if any(strip_line.startswith(card) for card in target_cards):
            capture = True
        
        # Stop capturing if we hit a namelist start or a closing slash
        elif strip_line.startswith("&") or strip_line == "/":
            capture = False
            
        if capture:
            relevant_cards += line

    # 4. K-Path Validation
    if "{crystal_b}" not in relevant_cards:
        print("\n" + "!"*70)
        print("[!] ERROR: High-symmetry K-path {crystal_b} not found in scf.dat.")
        print("[!] ACTION: Re-run Zone 1 (Module 102) with Option 'C'.")
        print("!"*70 + "\n")
        return

    # 5. Build Final Template
    bands_template = f"""&CONTROL
  calculation = 'bands'
  prefix = '{prefix}'
  pseudo_dir = '{pseudo_dir}'
  outdir = './outdir/'
  verbosity = 'high'
  tprnfor = .true.
  tstress = .true.
/
&SYSTEM
{system_vars.strip()}
  occupations = 'smearing'
  smearing = 'gauss'
  degauss = 1.0d-2
  ecutwfc = {ecutwfc}
  ecutrho = {ecutrho}
/
&ELECTRONS
  conv_thr = 1.0d-8
  mixing_beta = 0.5d0
/

{relevant_cards.strip()}
"""

    with open("bands.in", "w") as f:
        f.write(bands_template)

    print("\n[+] Success! Clean 'bands.in' generated.")

if __name__ == "__main__":
    run_bands_gen()

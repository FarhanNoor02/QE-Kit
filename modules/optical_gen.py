# modules/optical_gen.py
import os
import questionary
from utils import qe_xml_parser

def run_optical_gen():
    dat_file = "scf.dat"
    print("\n--- [205] Optical (Dielectric) Generation ---")

    # 1. Verification Step
    prefix = questionary.text("Enter calculation prefix:", default="aco").ask()
    calc_type, status = qe_xml_parser.check_calc_status(prefix)

    # XML check for correct pseudos
    with open(dat_file, 'r') as f:
        lines = f.readlines()
    
    content = "".join(lines)
    if "pbe-mt_fhi.UPF" not in content:
        print("\n" + "!"*75)
        print("[!] WARNING: NC-mt pseudopotentials not detected in scf.dat.")
        print("[!] Optical calculations require NC-mt potentials for accurate response.")
        print("[!] ACTION: Please re-run Module 101 selecting 'NC-mt' first.")
        print("!"*75 + "\n")
        if not questionary.confirm("Proceed with potentially incorrect pseudos?").ask():
            return

    # 2. User Physics Prompts
    etot_conv = questionary.text("Total energy conv threshold (etot_conv_thr):", default="1.0d-6").ask()
    forc_conv = questionary.text("Force conv threshold (forc_conv_thr):", default="1.0d-6").ask()
    electron_conv = questionary.text("conv_thr (electrons):", default="1.0d-8").ask()
    nbnd = questionary.text("Number of bands (nbnd):", default="60").ask()
    
    # Epsilon parameters
    intersmear = questionary.text("Interband broadening (intersmear, eV):", default="0.01").ask()
    intrasmear = questionary.text("Intraband broadening (intrasmear, eV):", default="0.001").ask()
    wmax = questionary.text("Max photon energy (wmax, eV):", default="20.0").ask()
    nw = questionary.text("Number of frequency points (nw):", default="500").ask()

    # 3. Surgical Extraction of Data from scf.dat
    pseudo_dir = "./"
    system_vars = ""
    relevant_cards = ""
    
    # Search for pseudo_dir
    for line in lines:
        if line.startswith("# PSEUDO_DIR_PATH:"):
            pseudo_dir = line.split(":")[1].strip()

    # Extract clean SYSTEM variables (only from the first block)
    in_system = False
    for line in lines:
        if "&SYSTEM" in line.upper(): in_system = True; continue
        if line.strip() == "/" and in_system: in_system = False; break
        if in_system:
            if not any(x in line for x in ["ecutwfc", "ecutrho", "occupations", "smearing", "degauss", "nbnd", "nosym", "noinv"]):
                system_vars += line

    # Extract ONLY relevant Cards
    target_cards = ["ATOMIC_SPECIES", "ATOMIC_POSITIONS", "CELL_PARAMETERS", "K_POINTS"]
    capture = False
    for line in lines:
        strip_line = line.strip().upper()
        if any(strip_line.startswith(card) for card in target_cards):
            capture = True
        elif strip_line.startswith("&") or strip_line == "/":
            capture = False
        
        if capture:
            relevant_cards += line

    # 4. Build scf_opt.in
    scf_opt_template = f"""&CONTROL
  calculation = 'scf'
  prefix = '{prefix}'
  pseudo_dir = '{pseudo_dir}'
  outdir = './outdir/'
  verbosity = 'high'
  tprnfor = .true.
  tstress = .true.
  etot_conv_thr = {etot_conv}
  forc_conv_thr = {forc_conv}
/
&SYSTEM
{system_vars.strip()}
  occupations = 'smearing'
  smearing = 'gauss'
  degauss = 0.02
  ecutwfc = 50
  ecutrho = 500
  nbnd = {nbnd}
  nosym = .true.
  noinv = .true.
/
&ELECTRONS
  conv_thr = {electron_conv}
  mixing_beta = 0.5d0
/

{relevant_cards.strip()}
"""

    # 5. Build epsilon.in
    epsilon_template = f"""&inputpp
    outdir = './outdir/'
    prefix = '{prefix}'
    calculation = 'eps'
/
&energy_grid
    smeartype = 'gauss'
    intersmear = {intersmear}d0
    intrasmear = {intrasmear}d0
    wmax = {wmax}d0
    wmin = 0.01d0
    nw = {nw}
    shift = 0d0
/
"""

    # 6. Writing files
    with open("scf_opt.in", "w") as f:
        f.write(scf_opt_template)
    with open("epsilon.in", "w") as f:
        f.write(epsilon_template)

    print("\n" + "="*60)
    print("[+] SUCCESS: Generated 'scf_opt.in' and 'epsilon.in'.")
    print("[!] ACTION: Ensure you have a UNIFORM K-MESH (Module 102 -> Option A).")
    print("="*60)

if __name__ == "__main__":
    run_optical_gen()

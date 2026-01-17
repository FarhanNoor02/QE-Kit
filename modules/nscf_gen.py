# modules/nscf_gen.py
import os
import questionary
from utils import qe_xml_parser  # Using the new utility tool

def run_nscf_gen():
    dat_file = "scf.dat"
    
    print("\n--- [202] NSCF Input Generation ---")
    
    # 1. Verification Step using the Utility Tool
    prefix = questionary.text("Enter calculation prefix:", default="ncn").ask()
    calc_type, status = qe_xml_parser.check_calc_status(prefix)

    # Logic Check: NSCF requires a finished SCF or VC-RELAX
    if status == "Not Found":
        print(f"\n[!] Error: No XML data found for prefix '{prefix}'.")
        print("    SCF calculation pending or 'outdir' is missing.")
        return
    
    if status == "Failed":
        print(f"\n[!] Error: The last calculation ('{calc_type}') failed.")
        print("    Please fix the SCF run before generating NSCF input.")
        return

    if calc_type not in ["scf", "vc-relax"]:
        print(f"\n[!] Warning: Found '{calc_type.upper()}' results, but NSCF usually follows SCF.")
        if not questionary.confirm("Do you want to proceed?").ask():
            return

    # 2. Extract Pseudo Path from scf.dat
    if not os.path.exists(dat_file):
        print(f"\n[!] Error: '{dat_file}' not found.")
        return

    pseudo_dir = "./"
    with open(dat_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("# PSEUDO_DIR_PATH:"):
            pseudo_dir = line.split(":")[1].strip()

    # 3. Prompts for User Input
    electron_conv = questionary.text("Electron conv threshold (conv_thr):", default="1.0d-8").ask()
    ecutwfc = questionary.text("ecutwfc (Ry):", default="50").ask()
    ecutrho = questionary.text("ecutrho (Ry):", default="375").ask()
    
    # 4. Extraction logic
    system_params = ""
    cards_content = ""
    in_system = False
    
    for line in lines:
        if line.startswith("&SYSTEM"):
            in_system = True
            continue
        if line.strip() == "/" and in_system:
            in_system = False
            continue
        
        if in_system:
            if not any(k in line for k in ["ecutwfc", "ecutrho", "occupations", "smearing", "degauss"]):
                system_params += line
        elif not line.startswith("#"): 
            cards_content += line

    # 5. Build the final nscf.in
    nscf_template = f"""&CONTROL
  calculation = 'nscf'
  prefix = '{prefix}'
  pseudo_dir = '{pseudo_dir}'
  outdir = './outdir/'
  verbosity = 'high'
/
&SYSTEM
{system_params.strip()}
  occupations = 'smearing'
  smearing = 'gauss'
  degauss = 1.0d-2
  ecutwfc = {ecutwfc}
  ecutrho = {ecutrho}
/
&ELECTRONS
  conv_thr = {electron_conv}
  mixing_beta = 0.5d0
/

{cards_content.strip()}
"""

    with open("nscf.in", "w") as f:
        f.write(nscf_template)

    print(f"\n[+] Success! 'nscf.in' generated for prefix '{prefix}'.")

if __name__ == "__main__":
    run_nscf_gen()

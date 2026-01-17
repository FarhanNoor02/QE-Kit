# modules/scf_gen.py
import os
import questionary

def run_scf_gen():
    dat_file = "scf.dat"
    
    if not os.path.exists(dat_file):
        print(f"\n[!] Error: '{dat_file}' not found. Run Pre-Processing first.")
        return

    # 1. Read scf.dat to extract Pseudo Path and Structural Blocks
    pseudo_dir = "./" # fallback
    with open(dat_file, 'r') as f:
        lines = f.readlines()

    # Search for the path tag we created in Module 101
    for line in lines:
        if line.startswith("# PSEUDO_DIR_PATH:"):
            pseudo_dir = line.split(":")[1].strip()

    print("\n--- [200] SCF Input Generation ---")
    print(f"Detected Pseudo Directory: {pseudo_dir}")

    # 2. Prompts for User Input
    prefix = questionary.text("Enter prefix (e.g. 'ncn'):", default="ncn").ask()
    etot_conv = questionary.text("Total energy conv threshold (etot_conv_thr):", default="1.0d-6").ask()
    forc_conv = questionary.text("Force conv threshold (forc_conv_thr):", default="1.0d-6").ask()
    electron_conv = questionary.text("SCF electron conv threshold (conv_thr):", default="1.0d-8").ask()
    ecutwfc = questionary.text("ecutwfc (Ry):", default="50").ask()
    ecutrho = questionary.text("ecutrho (Ry):", default="375").ask()

    # 3. Extraction logic for Structure and K-points
    # We need to extract: ibrav, nat, ntyp, celldm, and the cards
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
            # We filter out tags we are overriding (ecutwfc, etc)
            if not any(k in line for k in ["ecutwfc", "ecutrho", "occupations", "smearing", "degauss"]):
                system_params += line
        elif not line.startswith("#"): # Skip our metadata comments
            cards_content += line

    # 4. Build the final scf.in
    scf_template = f"""&CONTROL
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
{system_params.strip()}
  occupations = 'smearing'
  smearing = 'gauss'
  degauss = 1.0d-2
  ecutwfc = {ecutwfc}
  ecutrho = {ecutrho}
  la2F = .true.
/
&ELECTRONS
  conv_thr = {electron_conv}
  mixing_beta = 0.5d0
/

{cards_content.strip()}
"""

    with open("scf.in", "w") as f:
        f.write(scf_template)

    print("\n[+] Success! 'scf.in' has been generated.")

if __name__ == "__main__":
    run_scf_gen()

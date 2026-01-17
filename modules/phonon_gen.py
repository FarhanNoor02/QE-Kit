# modules/phonon_gen.py
import os
import questionary
from utils import qe_xml_parser

def run_phonon_gen():
    dat_file = "scf.dat"
    print("\n--- [206] Phonon Input Generation (ph.x) ---")

    # 1. Pre-flight Check
    prefix = questionary.text("Enter calculation prefix:", default="ncn").ask()
    calc_type, status = qe_xml_parser.check_calc_status(prefix)

    if status != "Success":
        print(f"\n[!] Error: No successful SCF found for prefix '{prefix}'.")
        return

    # 2. SURGICAL EXTRACTION OF MASSES (Order-Sensitive)
    if not os.path.exists(dat_file):
        print(f"[!] Error: {dat_file} missing.")
        return

    with open(dat_file, 'r') as f:
        lines = f.readlines()

    ordered_masses = []
    in_species = False
    
    for line in lines:
        upper_line = line.upper().strip()
        if upper_line.startswith("ATOMIC_SPECIES"):
            in_species = True
            continue
        
        if in_species:
            # If we hit an empty line or the next Card, stop
            if not upper_line or any(upper_line.startswith(c) for c in ["ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS"]):
                in_species = False
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                # Store (Element, Mass) to ensure we know exactly which is which
                ordered_masses.append({"element": parts[0], "mass": parts[1]})

    # 3. User Input for Ph Physics
    tr2_ph = questionary.text("Phonon tr2_ph:", default="1.0d-14").ask()
    q_grid = questionary.text("q-grid (nq1 nq2 nq3):", default="3 3 3").ask()
    fildyn = questionary.text("fildyn:", default=f"{prefix}.dyn").ask()
    
    nq1, nq2, nq3 = q_grid.split() if len(q_grid.split()) == 3 else ("3", "3", "3")

    # 4. Generate amass tags IN ORDER
    amass_block = ""
    print("\n[v] Mapping Atomic Masses from scf.dat:")
    for i, item in enumerate(ordered_masses):
        amass_block += f"  amass({i+1}) = {item['mass']}\n"
        print(f"    amass({i+1}) -> {item['element']}: {item['mass']}")

    # 5. Final Template Assembly
    phonon_template = f"""&inputph
  outdir = './outdir/'
  prefix = '{prefix}'
  tr2_ph = {tr2_ph}
{amass_block.rstrip()}
  nq1 = {nq1}
  nq2 = {nq2}
  nq3 = {nq3}
  fildyn = '{fildyn}'
  fildvscf = 'dvscf'
  ldisp = .true.
  trans = .true.
  electron_phonon = 'interpolated'
  el_ph_sigma = 0.02
  el_ph_nsigma = 10
  nmix_ph = 8
/
"""

    with open("ph.in", "w") as f:
        f.write(phonon_template)

    print("\n" + "="*60)
    print("[+] SUCCESS: 'ph.in' generated with strict amass ordering.")
    print("="*60)

if __name__ == "__main__":
    run_phonon_gen()

import os
import questionary
import re

def run_301_phonon_processing():
    print("\n--- [301] Phonon Post-Processing Tool ---")
    
    ph_in = "ph.in"
    scf_dat = "scf.dat"

    # 1. Verification Steps
    if not os.path.exists(ph_in):
        print("[!] Error: ph.in not found. Ensure the phonon calculation was set up.")
        return

    if not os.path.exists(scf_dat):
        print("[!] Error: scf.dat not found."); return

    # 2. Check for K-Path in scf.dat
    with open(scf_dat, 'r') as f:
        scf_content = f.read()

    if "K_POINTS {crystal_b}" not in scf_content:
        print("\n[!] Error: K-path not found in scf.dat.")
        print("    Repeat pre-processing with Kpath for band first (Module 102).")
        return

    # 3. Scan ph.in for prefix, fildyn, and amass values
    with open(ph_in, 'r') as f:
        ph_text = f.read()

    # Extraction using Regex
    prefix_match = re.search(r"prefix\s*=\s*['\"]([^'\"]+)['\"]", ph_text)
    fildyn_match = re.search(r"fildyn\s*=\s*['\"]([^'\"]+)['\"]", ph_text)
    amass_matches = re.findall(r"amass\((\d+)\)\s*=\s*([\d\.]+)", ph_text)

    if not prefix_match or not fildyn_match:
        print("[!] Error: Could not extract prefix or fildyn from ph.in.")
        return

    prefix = prefix_match.group(1)
    fildyn = fildyn_match.group(1)
    amass_dict = {idx: val for idx, val in amass_matches}

    # 4. User Prompts for DOS and Grids
    nk1 = questionary.text("Enter nk1 for Phonon DOS:", default="10").ask()
    nk2 = questionary.text("Enter nk2 for Phonon DOS:", default="10").ask()
    nk3 = questionary.text("Enter nk3 for Phonon DOS:", default="10").ask()
    ndos = questionary.text("Enter ndos:", default="50").ask()

    # 5. Generate q2r.in
    q2r_content = f"""&input
  zasr='crystal', 
  fildyn='{fildyn}', 
  flfrc='{prefix}.fc', 
  la2F=.true.
/
"""
    with open("q2r.in", "w") as f:
        f.write(q2r_content)
    print("[+] Generated q2r.in")

    # 6. Generate matdyn.in (Phonon DOS)
    amass_str = "\n".join([f"  amass({i})={v}" for i, v in amass_dict.items()])
    matdyn_dos = f"""&input
  asr='crystal',
{amass_str}
  flfrc='{prefix}.fc',
  flfrq='{prefix}.freq',
  la2F=.true.,
  dos=.true.,
  fldos='phonon.dos',
  nk1={nk1}, nk2={nk2}, nk3={nk3}, ndos={ndos}
/
"""
    with open("matdyn.in", "w") as f:
        f.write(matdyn_dos)
    print("[+] Generated matdyn.in")

    # 7. Generate matdyndisp.in (Phonon Bands/Dispersion)
    # Extract K-Path block from scf.dat
    kpath_lines = []
    found_k = False
    lines = scf_content.splitlines()
    for i, line in enumerate(lines):
        if "K_POINTS {crystal_b}" in line:
            found_k = True
            # The next line is the number of points
            num_points = lines[i+1].strip()
            kpath_lines.append(num_points)
            # Add the actual points until a blank line or namelist
            for j in range(i+2, len(lines)):
                if not lines[j].strip() or lines[j].strip().startswith(("&", "/")):
                    break
                kpath_lines.append(lines[j])
            break

    matdyn_disp = f"""&input
  asr='crystal',
{amass_str}
  flfrc='{prefix}.fc',
  flfrq='{prefix}.freq',
  q_in_band_form=.true.,
  q_in_cryst_coord=.true.
/
"""
    with open("matdyndisp.in", "w") as f:
        f.write(matdyn_disp)
        f.write("\n".join(kpath_lines) + "\n")
    
    print("[+] Generated matdyndisp.in with extracted K-path")
    print("\n[v] Post-processing files ready. Sequence: q2r.x -> matdyn.x (DOS) -> matdyn.x (Dispersion)")


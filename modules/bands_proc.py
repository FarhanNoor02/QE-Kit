# modules/bands_proc.py
import os
import re
import subprocess
from utils import qe_xml_parser

def run_306_bands_processing(automated=False):
    print("\n--- [306] Band Structure Post-Processing (bands.x) ---")
    
    bands_in = "bands.in"
    nscf_out = "nscf.out"
    bandpp_in = "bandpp.in"
    bandpp_out = "bandpp.out"

    if not os.path.exists(bands_in):
        print(f"[!] Error: '{bands_in}' not found in the current directory.")
        return

    # 1. Extract prefix and outdir from bands.in FIRST
    # We need these variables to feed into the XML parser
    prefix = "pwscf"
    outdir = "./outdir/"
    
    with open(bands_in, 'r') as f:
        bands_content = f.read()
        
    prefix_match = re.search(r"prefix\s*=\s*['\"]([^'\"]+)['\"]", bands_content, re.IGNORECASE)
    outdir_match = re.search(r"outdir\s*=\s*['\"]([^'\"]+)['\"]", bands_content, re.IGNORECASE)
    
    if prefix_match: prefix = prefix_match.group(1)
    if outdir_match: outdir = outdir_match.group(1)

    # 2. Strict XML Health Verification
    print(f"[*] Verifying run health for prefix '{prefix}' in '{outdir}'...")
    calc_type, status = qe_xml_parser.check_calc_status(prefix, outdir)
    
    if status == "Not Found":
        print(f"[!] CRITICAL: XML schema not found.")
        print("    Ensure pw.x (calculation='bands') completed successfully.")
        return
    elif status == "Corrupted":
        print(f"[!] CRITICAL: XML schema is corrupted.")
        return
    elif status == "Failed":
        print(f"[!] CRITICAL: Previous calculation crashed (QE Exit Status != 0).")
        print("    Please review your pw.x output files for errors.")
        return
    
    if calc_type != "bands":
        print(f"[!] Warning: The XML reports the last run was '{calc_type}', not 'bands'.")
        print("    This could mean the save directory was overwritten. Proceeding with caution...\n")
    else:
        print("[+] XML schema verified successfully. (Status: Success)\n")

    # 3. Prepare bandpp.in for bands.x
    bandpp_content = f"""&bands
  outdir='{outdir}'
  prefix='{prefix}'
  filband='{prefix}.band.dat'
/
"""
    with open(bandpp_in, 'w') as f:
        f.write(bandpp_content)
    print(f"[*] Generated '{bandpp_in}'.")

    # 4. Execute bands.x 
    print("[*] Executing bands.x post-processing...")
    try:
        cmd = f"bands.x < {bandpp_in} > {bandpp_out}"
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        print(f"[!] Error executing bands.x. Check '{bandpp_out}' for detailed errors.")
        return

    if not os.path.exists(bandpp_out):
        print(f"[!] Error: '{bandpp_out}' was not generated.")
        return

    # 5. Parse bandpp.out for X-Coordinates
    x_coords = []
    with open(bandpp_out, 'r') as f:
        for line in f:
            if "high-symmetry point:" in line and "x coordinate" in line:
                match = re.search(r"x coordinate\s+([0-9\.]+)", line)
                if match:
                    x_coords.append(match.group(1))

    # 6. Parse bands.in for K-point Labels
    labels = []
    lines = bands_content.splitlines()
    in_kpoints = False
    kpts_count = 0
    
    for line in lines:
        if "K_POINTS" in line.upper():
            in_kpoints = True
            continue
        if in_kpoints:
            if kpts_count == 0 and len(line.split()) == 1 and line.strip().isdigit():
                kpts_count = int(line.strip())
                continue
                
            label_match = re.search(r"!\s*([A-Za-z0-9\_]+)", line)
            if label_match:
                labels.append(label_match.group(1).upper())
            elif line.strip() != "":
                labels.append("?") 

            if len(labels) == kpts_count and kpts_count > 0:
                break

    # 7. Extract Fermi Energy from nscf.out
    e_fermi = "Not found"
    if os.path.exists(nscf_out):
        with open(nscf_out, 'r') as f:
            for line in f:
                if "the Fermi energy is" in line:
                    match = re.search(r"the Fermi energy is\s+([0-9\.\-]+)\s+ev", line)
                    if match:
                        e_fermi = match.group(1)
                    break
    else:
        print(f"[!] '{nscf_out}' not found. Cannot extract Fermi energy.")

    # 8. Print Console Summary
    print("\n" + "="*55)
    print(" 📊 BAND STRUCTURE DATA EXTRACTED")
    print("="*55)
    print(f"[*] Fermi Energy (E_F): {e_fermi} eV\n")
    print("[*] High-Symmetry K-Point Mapping:")
    print("-" * 55)
    
    kpath_str = []
    for i in range(min(len(labels), len(x_coords))):
        print(f"    {labels[i]:<10} ->   x = {x_coords[i]}")
        kpath_str.append(labels[i])
        
    print("-" * 55)
    print(f"[*] K-Path sequence: {' - '.join(kpath_str)}")
    print("="*55 + "\n")
    print(f"[+] Plottable dataset written to: '{prefix}.band.dat.gnu'")

if __name__ == "__main__":
    run_306_bands_processing()

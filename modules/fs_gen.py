import os
import re
from utils.qe_xml_parser import check_calc_status

def extract_prefix():
    """Scans local input files to accurately extract the user's defined prefix."""
    for file in ["scf.in", "nscf.in", "bands.in"]:
        if os.path.exists(file):
            with open(file, 'r') as f:
                for line in f:
                    match = re.search(r"prefix\s*=\s*['\"]([^'\"]+)['\"]", line, re.IGNORECASE)
                    if match:
                        return match.group(1)
    return "pwscf"

def run_207_fs_gen():
    """Generates the fs.in input file for Quantum ESPRESSO's fs.x executable."""
    print("\n--- [Zone 2] Fermi Surface Input Generator (fs.x) ---")
    
    prefix = extract_prefix()
    print(f"[*] Extracted calculation prefix: '{prefix}'")
    
    # Utilize the XML parser to verify the exact state of the calculation
    print("[*] Verifying calculation state via XML schema...")
    calc_type, status = check_calc_status(prefix)
    
    # 1. Validation Logic
    if status == "Success":
        if calc_type == "nscf":
            print(f"  [✓] Verified successful '{calc_type}' calculation. Ready for fs.x.")
        elif calc_type == "scf":
            print(f"  [!] Warning: Last successful run was '{calc_type}'.")
            print("      fs.x requires a dense, unshifted NSCF grid. Please run NSCF first.")
        else:
            print(f"  [!] Warning: Last successful run was '{calc_type}'. Expected 'nscf'.")
    elif status == "Failed":
        print(f"  [!] CRITICAL: The previous '{calc_type}' calculation failed (exit_status != 0).")
        print("      Do not run fs.x until the DFT calculation is resolved.")
    else:
        print("  [!] Warning: Could not locate or parse XML schema. Proceeding blindly...")

    # 2. Build Input
    fs_content = f"""&fermi
  outdir = './outdir'
  prefix = '{prefix}'
/
"""
    
    with open("fs.in", 'w') as f:
        f.write(fs_content)
        
    print("\n[+] Successfully generated 'fs.in'!")
    print("    Run this using: fs.x < fs.in > fs.out\n")

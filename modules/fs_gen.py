# modules/fs_gen.py
import os
import re
from utils.qe_xml_parser import check_calc_status

def extract_prefix():
    for file in ["scf.in", "nscf.in", "bands.in"]:
        if os.path.exists(file):
            with open(file, 'r') as f:
                match = re.search(r"prefix\s*=\s*['\"]([^'\"]+)['\"]", f.read(), re.IGNORECASE)
                if match: return match.group(1)
    return "pwscf"

def run_207_fs_gen():
    """Generates the fs.in input file."""
    print("\n--- [Zone 2] Fermi Surface Input Generator (fs.x) ---")
    prefix = extract_prefix()
    
    # Validation Logic
    _, status = check_calc_status(prefix)
    if status == "Success":
        print(f"  [✓] Verified calculation state for prefix: '{prefix}'")
    else:
        print(f"  [!] Warning: Calculation state for '{prefix}' is '{status}'. Proceed with caution.")

    # Build Input
    fs_content = f"&fermi\n  outdir = './outdir'\n  prefix = '{prefix}'\n/\n"
    with open("fs.in", "w") as f:
        f.write(fs_content)
    print("[+] Generated 'fs.in'. Run with: mpirun -np <N> fs.x < fs.in > fs.out")

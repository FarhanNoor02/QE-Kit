# modules/thermopw_gen.py
import os
import re
from utils.qe_xml_parser import check_calc_status

def run_208_thermopw_gen():
    """Generates the thermo_control input file for thermo_pw.x."""
    print("\n--- [Zone 2] Thermo_PW Input Generator ---")

    if not os.path.exists("scf.in"):
        print("[!] Error: 'scf.in' not found.")
        return

    # Verify SCF run status using your XML parser utility
    with open("scf.in", 'r') as f:
        prefix_match = re.search(r"prefix\s*=\s*['\"]([^'\"]+)['\"]", f.read(), re.IGNORECASE)
        prefix = prefix_match.group(1) if prefix_match else "pwscf"

    _, status = check_calc_status(prefix)
    if status != "Success":
        print(f"[!] Warning: No successful SCF run verified for prefix '{prefix}'.")
    else:
        print(f"[+] Verified successful SCF run for '{prefix}'.")
    
    # Get pressure input
    try:
        pressure = float(input("[?] Enter pressure in kilobar [default: 0.0]: ") or "0.0")
    except ValueError:
        print("[!] Invalid input. Defaulting to 0.0 kbar.")
        pressure = 0.0
    
    # Write the thermo_control file
    with open("thermo_control", 'w') as f:
        f.write(f"&INPUT_THERMO\n  what = 'mur_lc_elastic_constants',\n  pressure = {pressure:.2f}\n/\n")
    
    print("[+] 'thermo_control' generated successfully.")

import os
import re
import numpy as np

# Physical constants for derivations
H_PLANCK = 4.135667696e-15  # eV·s
C_LIGHT = 2.99792458e10     # cm/s
PI = np.pi

def load_eps_data_averaged(filename):
    """
    Reads QE optical files (epsr, epsi, eels), handles malformed formatting 
    via Regex, and returns the Isotropic Average of the x, y, and z components.
    """
    energy_vals = []
    avg_vals = []
    
    if not os.path.exists(filename):
        return None, None
        
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines, comments, and lines with overflow asterisks
            if not line or line.startswith('#') or '*' in line:
                continue
                
            # Regex fix: insert space before minus signs (QE formatting bug)
            clean_line = re.sub(r'(?<!^)(?<![eE])-', ' -', line)
            
            try:
                parts = clean_line.split()
                # Col 0: Energy, Col 1,2,3: x,y,z components
                e_val = float(parts[0])
                val_x, val_y, val_z = float(parts[1]), float(parts[2]), float(parts[3])
                
                # Compute Isotropic Average ( (x+y+z)/3 )
                avg_res = (val_x + val_y + val_z) / 3.0
                
                energy_vals.append(e_val)
                avg_vals.append(avg_res)
            except (ValueError, IndexError):
                continue
                
    return np.array(energy_vals), np.array(avg_vals)

def run_303_optical_processing():
    print("\n--- [303] Optical Properties: Isotropic Averaging & Derived Constants ---")
    
    prefix = input("Enter the calculation prefix (e.g., nc): ").strip()
    if not prefix:
        print("[!] Prefix cannot be empty."); return

    # Define standard filenames based on prefix
    file_r = f"epsr_{prefix}.dat"
    file_i = f"epsi_{prefix}.dat"
    file_eels = f"eels_{prefix}.dat"

    # 1. Load and Average Dielectric Parts
    print(f" [+] Reading and Averaging {file_r} and {file_i}...")
    energy, eps1 = load_eps_data_averaged(file_r)
    _, eps2 = load_eps_data_averaged(file_i)

    if eps1 is None or eps2 is None:
        print(f"[!] Error: Mandatory files {file_r} or {file_i} not found.")
        return

    # 2. Load and Average EELS if present
    eels_exists = os.path.exists(file_eels)
    _, eels_avg = load_eps_data_averaged(file_eels) if eels_exists else (None, None)

    # 3. Synchronize array lengths
    # We use the minimum length to avoid index errors if files differ by a line
    list_lengths = [len(eps1), len(eps2)]
    if eels_exists: list_lengths.append(len(eels_avg))
    
    min_len = min(list_lengths)
    energy = energy[:min_len]
    eps1 = eps1[:min_len]
    eps2 = eps2[:min_len]

    print(" [+] Calculating n, k, R, and alpha (Absorption)...")
    
    # Mathematical Derivations
    eps_mag = np.sqrt(eps1**2 + eps2**2)
    n = np.sqrt((eps_mag + eps1) / 2)
    k = np.sqrt((eps_mag - eps1) / 2)
    R = ((n - 1)**2 + k**2) / ((n + 1)**2 + k**2)
    alpha = (4 * PI * k * energy) / (H_PLANCK * C_LIGHT)

    # 4. Consolidate and Save
    output_fn = f"optical_final_{prefix}.dat"
    
    if eels_exists:
        eels_avg = eels_avg[:min_len]
        data_stack = np.column_stack((energy, eps1, eps2, n, k, R, alpha, eels_avg))
        header = "Energy(eV)  Eps_Real    Eps_Imag    n           k           R           alpha(cm-1) EELS_Avg"
    else:
        data_stack = np.column_stack((energy, eps1, eps2, n, k, R, alpha))
        header = "Energy(eV)  Eps_Real    Eps_Imag    n           k           R           alpha(cm-1)"

    np.savetxt(output_fn, data_stack, header=header, fmt='%12.6f')
    
    print(f"\n[v] SUCCESS: Processed data saved to '{output_fn}'")
    if eels_exists:
        print(" [+] Note: EELS data was detected and successfully averaged.")

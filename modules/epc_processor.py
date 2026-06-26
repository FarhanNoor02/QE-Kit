import os
import glob
import numpy as np
from scipy.integrate import simpson, cumulative_trapezoid

def parse_a2fdos_file(filepath):
    """Parses the a2F.dos file, ignoring headers and trailing text."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            # Strip leading spaces to properly catch the '#' comment lines
            if line.strip().startswith('#'):
                continue
            
            parts = line.split()
            if not parts:
                continue
                
            try:
                # Check if the line contains numeric data
                float(parts[0])
                data.append([float(x) for x in parts])
            except ValueError:
                # Safely break when hitting the "lambda = ..." footer
                break
                
    return np.array(data)

def run_304_epc_processor():
    """Calculates lambda, w_ln, and Tc from all a2F.dos files in the directory."""
    print("\n--- [Zone 3] Electron-Phonon Coupling & Tc Estimator ---")
    
    a2f_files = sorted(glob.glob("a2F.dos*"))
    
    if not a2f_files:
        print("[!] Error: No 'a2F.dos*' files found in the current directory.")
        return

    try:
        mu_star = float(input("[?] Enter estimated Coulomb pseudopotential (μ*) [default: 0.10]: ") or "0.10")
    except ValueError:
        mu_star = 0.10
        print("[*] Invalid input. Defaulting to μ* = 0.10")

    # --- Unit Conversion Constants ---
    RY_TO_EV = 13.6056980659
    RY_TO_K = 157887.3
    RY_TO_CM1 = 109737.31  # User-defined conversion factor

    print("\n" + "="*85)
    print(f"{'Filename':<15} | {'λ':<8} | {'ω_ln (cm^-1)':<12} | {'ω_ln (K)':<10} | {'Tc (K)':<10}")
    print("="*85)

    for file in a2f_files:
        data = parse_a2fdos_file(file)
        
        if data.size == 0:
            print(f"{file:<15} | [!] Empty or unreadable data")
            continue
            
        omega_ry = data[:, 0]
        a2f = data[:, 1]
        
        # Avoid division by zero at omega = 0 by filtering out near-zero frequencies
        valid_idx = omega_ry > 1e-6
        w_ry = omega_ry[valid_idx]
        a2f_val = a2f[valid_idx]
        
        # 1. Calculate Total Lambda
        integrand_lambda = a2f_val / w_ry
        lambda_tot = 2.0 * simpson(integrand_lambda, x=w_ry)
        
        # 2. Calculate Logarithmic Average Frequency (w_ln)
        integrand_ln = (a2f_val / w_ry) * np.log(w_ry)
        w_ln_ry = np.exp((2.0 / lambda_tot) * simpson(integrand_ln, x=w_ry))
        
        w_ln_K = w_ln_ry * RY_TO_K
        w_ln_cm1 = w_ln_ry * RY_TO_CM1 # Converted to wavenumbers
        
        # 3. Calculate Tc (Allen-Dynes Equation)
        if lambda_tot > mu_star:
            numerator = -1.04 * (1.0 + lambda_tot)
            denominator = lambda_tot - mu_star * (1.0 + 0.62 * lambda_tot)
            tc_k = (w_ln_K / 1.2) * np.exp(numerator / denominator)
        else:
            tc_k = 0.0
            
        print(f"{file:<15} | {lambda_tot:<8.4f} | {w_ln_cm1:<12.2f} | {w_ln_K:<10.2f} | {tc_k:<10.3f}")

        # 4. Calculate and export frequency-dependent lambda(omega)
        w_ev = w_ry * RY_TO_EV
        w_cm1_arr = w_ry * RY_TO_CM1
        lambda_omega = 2.0 * cumulative_trapezoid(integrand_lambda, w_ry, initial=0)
        
        suffix = file.replace('a2F.dos', '')
        if suffix.startswith('.'):
            suffix = suffix[1:]
        out_name = f"lambda_omega_{suffix}.dat" if suffix else "lambda_omega.dat"
        
        # Stack 3 columns: eV, cm^-1, lambda
        output_data = np.column_stack((w_ev, w_cm1_arr, lambda_omega))
        np.savetxt(out_name, output_data, header='omega(eV) omega(cm^-1) lambda(omega)', fmt='%.6f')

    print("="*85)
    print("[+] Frequency-dependent λ(ω) exported as .dat files for Zone 4 plotting.\n")

if __name__ == "__main__":
    # Optional: allows the script to be run directly for testing
    run_304_epc_processor()

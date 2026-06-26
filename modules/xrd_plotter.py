import glob
import os
import matplotlib.pyplot as plt

def run_405_xrd_plotter():
    """Reads a CIF file, calculates the XRD pattern, and exports a clean plot."""
    print("\n--- [Zone 4] Calculated XRD Pattern Visualizer ---")
    
    try:
        from pymatgen.core import Structure
        from pymatgen.analysis.diffraction.xrd import XRDCalculator
    except ImportError:
        print("[!] Error: 'pymatgen' is not installed. Please add it to your dependencies.")
        return

    # 1. Automatically detect the first .cif file in the current directory
    cif_files = glob.glob("*.cif")

    if not cif_files:
        print("[!] Error: No .cif files found in the current directory.")
        return

    detected_cif = cif_files[0]
    print(f"[*] Found CIF file: '{detected_cif}' -> Processing...\n")

    # 2. Load the structure
    try:
        structure = Structure.from_file(detected_cif)
        formula = structure.composition.reduced_formula
        print(f"[*] Loaded structure formula: {formula}")
    except Exception as e:
        print(f"[!] Error reading {detected_cif}: {e}")
        return

    # 3. Initialize XRD calculator (Cu-Kalpha radiation)
    print("[*] Calculating diffraction pattern (Cu-Kα)...")
    xrd_calc = XRDCalculator(wavelength="CuKa")
    pattern = xrd_calc.get_pattern(structure)

    # 4. Extract and output ALL data to the console
    print("\n" + "="*55)
    print(f"{'2-Theta (deg)':<15} | {'Intensity (%)':<15} | {'HKL Miller Indices'}")
    print("="*55)

    for i in range(len(pattern.x)):
        two_theta = pattern.x[i]
        intensity = pattern.y[i]
        
        # Extract Miller indices from pymatgen data structure
        hkl_list = pattern.hkls[i]
        hkl_str = ", ".join([str(d['hkl']) for d in hkl_list])
        
        print(f"{two_theta:<15.4f} | {intensity:<15.2f} | ({hkl_str})")

    print("="*55 + "\n")

    # 5. Generate and save the completely clean plot
    print("[*] Generating and saving clean XRD plot...")
    plt.figure(figsize=(10, 5))

    plt.vlines(pattern.x, 0, pattern.y, colors='black', linewidth=1.5)

    plt.title(f"Calculated XRD Pattern of {formula}", fontsize=14, fontweight='bold', pad=15)
    plt.xlabel(r"$2\theta\ (^\circ)$", fontsize=12)
    plt.ylabel("Intensity (a.u.)", fontsize=12)

    plt.xlim(10, 90)
    plt.ylim(0, 105)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.grid(axis='x', linestyle=':', alpha=0.6)

    output_image = f"{formula}_clean_xrd.png"
    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[+] Success! Clean plot saved as: {output_image}\n")

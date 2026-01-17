# modules/structure2in.py
import subprocess
import os
import glob

def run_structure2in():
    """
    Finds a CIF file in the current directory and converts it to scf.dat
    using cif2cell for Quantum ESPRESSO.
    """
    # 1. Search for CIF files in the current directory
    cif_files = glob.glob("*.cif")
    
    if not cif_files:
        print("Error: No .cif files found in the current directory.")
        return False
    
    # For now, we take the first one found, or ask if multiple exist
    if len(cif_files) > 1:
        print("Multiple CIF files found:")
        for i, f in enumerate(cif_files):
            print(f" {i}) {f}")
        idx = int(input("Select file index: "))
        selected_cif = cif_files[idx]
    else:
        selected_cif = cif_files[0]

    print(f"--- Processing: {selected_cif} ---")

    # 2. Define the output filename
    output_file = "scf.dat"

    # 3. Call cif2cell via subprocess
    # -p quantum-espresso: specifies the output format
    # --setup-all: generates a more complete input (pseudo names, etc.)
    # -o: specifies the output filename
    try:
        command = [
            "cif2cell", 
            selected_cif, 
            "-p", "quantum-espresso", 
            "-o", output_file
        ]
        
        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"Successfully created {output_file}")
            return True
        else:
            print(f"Error running cif2cell: {result.stderr}")
            return False

    except FileNotFoundError:
        print("Error: 'cif2cell' is not installed or not in your PATH.")
        print("Install it with: pip install cif2cell")
        return False

# This allows you to test the module individually
if __name__ == "__main__":
    run_structure2in()

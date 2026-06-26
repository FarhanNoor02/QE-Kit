import glob
from ase.io import read

def run_404_bz_plotter():
    """Reads a CIF file in the directory and plots the Brillouin Zone."""
    print("\n--- [Zone 4] Brillouin Zone Visualizer ---")
    
    cif_files = glob.glob("*.cif")
    
    if not cif_files:
        print("[!] Error: No .cif files found in the current directory.")
        return
        
    if len(cif_files) > 1:
        print(f"[*] Found multiple CIF files: {cif_files}")
        print(f"[*] Proceeding with the first one: '{cif_files[0]}'\n")
        
    cif_path = cif_files[0]
    
    try:
        atoms = read(cif_path)
    except Exception as e:
        print(f"[!] Error reading {cif_path}: {e}")
        return
        
    lattice = atoms.cell.get_bravais_lattice()
    
    print(f"========================================")
    print(f"  Processing: {cif_path}")
    print(f"  Lattice:    {lattice.description()}")
    print(f"========================================")
    print("[*] Launching interactive 3D plot... (Close the plot window to return to menu)")
    
    # Generate and display the interactive 3D Brillouin Zone plot
    lattice.plot_bz(show=True, interactive=True)

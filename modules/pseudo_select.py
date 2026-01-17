# modules/pseudo_select.py
import os
import questionary
from utils import config_manager

def run_pseudo_select():
    scf_dat = "scf.dat"
    
    # 1. Fetch Global Path
    PSEUDO_LIB_PATH = config_manager.get_pseudo_dir()
    
    if not PSEUDO_LIB_PATH:
        print("\n[!] Error: Parent Pseudo Directory not configured in Settings.")
        return

    if not os.path.exists(scf_dat):
        print(f"\n[!] Error: '{scf_dat}' not found. Run Module 100 first.")
        return

    print("\n--- [101] Smart Pseudopotential Selection ---")

    # 2. Navigate the QE-POTCAR Hierarchy
    functional = questionary.select(
        "Select Exchange Functional:",
        choices=[
            {"name": "PBE (GGA)", "value": "PBE"},
            {"name": "LDA (Local Density)", "value": "LDA"},
            {"name": "NC-mt (Norm-Conserving)", "value": "NC-mt"}
        ]
    ).ask()

    # NC-mt is top-level in your tree
    if functional == "NC-mt":
        target_dir = os.path.join(PSEUDO_LIB_PATH, "NC-mt")
    else:
        calc_type = questionary.select(
            "Select Calculation Type:",
            choices=[
                {"name": "Scalar Relativistic", "value": "scalar"},
                {"name": "Full Relativistic", "value": "full"}
            ]
        ).ask()

        pp_type = questionary.select(
            "Select Formalism:",
            choices=[
                {"name": "PAW (Projector Augmented Wave)", "value": "PAW"},
                {"name": "USPP (Ultrasoft)", "value": "USPP"}
            ]
        ).ask()
        
        target_dir = os.path.join(PSEUDO_LIB_PATH, functional, calc_type, pp_type)
    
    if not os.path.exists(target_dir):
        print(f"\n[!] Error: Directory not found: {target_dir}")
        return

    # 3. Read and Update scf.dat with Smart File Matching
    with open(scf_dat, 'r') as f:
        lines = f.readlines()

    new_lines = []
    # Inject the directory metadata for other modules to use
    new_lines.append(f"# PSEUDO_DIR_PATH: {target_dir}\n")
    
    # Get list of all available UPFs in the target directory
    available_upfs = [f for f in os.listdir(target_dir) if f.upper().endswith(".UPF")]
    
    in_species_block = False
    
    print(f"\n[v] Scanning library: {target_dir}")

    for line in lines:
        if line.startswith("# PSEUDO_DIR_PATH:"): continue
            
        if line.strip().startswith("ATOMIC_SPECIES"):
            in_species_block = True
            new_lines.append(line)
            continue
        
        if in_species_block:
            parts = line.split()
            # Exit block detection
            if not parts or any(card in line for card in ["ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS"]):
                in_species_block = False
                new_lines.append(line)
            else:
                element = parts[0]
                mass = parts[1]
                
                # SMART MATCHING: Find the file starting with 'Element.'
                # Example: 'Al.' matches 'Al.pz-n-rrkjus_psl.0.1.UPF'
                match = None
                for upf in available_upfs:
                    if upf.startswith(f"{element}."):
                        match = upf
                        break
                
                if match:
                    new_lines.append(f"  {element:<5} {mass:<10} {match}\n")
                    print(f"    Mapped: {element} -> {match}")
                else:
                    print(f"    [!] Warning: No pseudo found for {element} in this folder.")
                    new_lines.append(line) # Keep original if not found
        else:
            new_lines.append(line)

    # 4. Save Updates
    with open(scf_dat, 'w') as f:
        f.writelines(new_lines)

    print(f"\n[+] Success! {scf_dat} updated with specific QE-POTCAR filenames.")

if __name__ == "__main__":
    run_pseudo_select()

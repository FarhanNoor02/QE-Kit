# modules/pseudo_select.py
import os
import questionary
from utils import config_manager

def run_pseudo_select():
    target_file = "scf.dat"
    
    if not os.path.exists(target_file):
        print(f"\n[!] Error: '{target_file}' not found. Run Module 100 first.")
        return

    print(f"\n--- [101] Smart Pseudopotential Selection ---")

    # 1. Identify required elements from scf.dat
    with open(target_file, 'r') as f:
        lines = f.readlines()
        
    elements = []
    in_species_block = False
    for line in lines:
        if line.strip().startswith("ATOMIC_SPECIES"):
            in_species_block = True
            continue
        if in_species_block:
            # Stop parsing elements if we hit the next block
            if not line.strip() or any(card in line for card in ["ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS"]):
                in_species_block = False
                continue
            parts = line.split()
            if len(parts) >= 2:
                elements.append(parts[0]) # Extract the element symbol (e.g., 'Nb', 'C')
                
    if not elements:
        print("[!] Error: No elements found under ATOMIC_SPECIES in scf.dat.")
        return

    print(f"[*] Elements detected: {', '.join(elements)}")

    # 2. Scan current directory for .UPF files
    local_upfs = [f for f in os.listdir('.') if f.upper().endswith(".UPF")]
    local_matches = {}
    
    for el in elements:
        # Match element.UPF or element_UPF (e.g., 'Nb.' matches 'Nb.pbe-spn-rrkjus_psl.1.0.0.UPF')
        match = next((upf for upf in local_upfs if upf.startswith(f"{el}.") or upf.startswith(f"{el}_")), None)
        if match:
            local_matches[el] = match

    # 3. Decision Logic: Local vs Global Library
    if len(local_matches) == len(elements):
        print("\n[+] All required pseudopotentials found in the current directory.")
        pseudo_path = "./"
        matches_to_use = local_matches
    else:
        missing = [el for el in elements if el not in local_matches]
        print(f"\n[!] Missing local UPFs for: {', '.join(missing)}")
        print("[*] Accessing global library...")
        
        PSEUDO_LIB_PATH = config_manager.get_pseudo_dir()
        if not PSEUDO_LIB_PATH:
            print("\n[!] Error: Parent Pseudo Directory not configured in Settings.")
            return

        functional = questionary.select(
            "Select Exchange Functional:",
            choices=[
                {"name": "PBE (GGA)", "value": "PBE"},
                {"name": "LDA (Local Density)", "value": "LDA"},
                {"name": "NC-mt (Norm-Conserving)", "value": "NC-mt"}
            ]
        ).ask()

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
            
        print(f"\n[v] Scanning library: {target_dir}")
        lib_upfs = [f for f in os.listdir(target_dir) if f.upper().endswith(".UPF")]
        
        lib_matches = {}
        for el in elements:
            match = next((upf for upf in lib_upfs if upf.startswith(f"{el}.") or upf.startswith(f"{el}_")), None)
            if match:
                lib_matches[el] = match
            else:
                print(f"    [!] Warning: No pseudo found for {el} in library.")
                
        pseudo_path = target_dir
        matches_to_use = lib_matches

    # 4. Write updates back to scf.dat
    new_lines = []
    # Inject the PSEUDO_DIR_PATH metadata exactly once at the top
    new_lines.append(f"# PSEUDO_DIR_PATH: {pseudo_path}\n")
    
    in_species_block = False
    for line in lines:
        if line.startswith("# PSEUDO_DIR_PATH:"): 
            continue # Remove old path metadata
            
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
                continue
                
            element = parts[0]
            mass = parts[1]
            
            # Map the actual UPF filename instead of _PSEUDO
            if element in matches_to_use:
                upf_name = matches_to_use[element]
                new_lines.append(f"  {element:<4} {mass:<10} {upf_name}\n")
                print(f"  Mapped: {element} -> {upf_name}")
            else:
                new_lines.append(line) # Fallback to original if not found
        else:
            new_lines.append(line)

    with open(target_file, 'w') as f:
        f.writelines(new_lines)
        
    print(f"\n[+] Success! '{target_file}' updated.")
    print(f"[*] pseudo_dir tag set to: '{pseudo_path}' (scf_gen.py will read this)")

if __name__ == "__main__":
    run_pseudo_select()

# modules/pseudo_select.py
import os
import re
import questionary
from utils import config_manager

def extract_required_upfs(filepath):
    """Parses ATOMIC_SPECIES block to extract required UPF filenames."""
    upfs = []
    in_species_block = False
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip().startswith("ATOMIC_SPECIES"):
                in_species_block = True
                continue
            
            if in_species_block:
                # Stop if we hit another block header
                if not line.strip() or any(card in line for card in ["ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS"]):
                    in_species_block = False
                    continue
                
                parts = line.split()
                if len(parts) >= 3:
                    upfs.append(parts[2])
    return upfs

def update_pseudo_dir(filepath, new_dir_path):
    """Updates or injects the pseudo_dir tag in the &CONTROL block."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Regex to find and replace existing pseudo_dir
    pattern = r"(pseudo_dir\s*=\s*)(['\"])(.*?)\2"
    
    if re.search(pattern, content, re.IGNORECASE):
        new_content = re.sub(pattern, rf"\g<1>'{new_dir_path}'", content, flags=re.IGNORECASE)
    else:
        # Inject if missing (assumes &CONTROL exists)
        new_content = re.sub(r"(&CONTROL\s*\n)", rf"\1  pseudo_dir='{new_dir_path}',\n", content, flags=re.IGNORECASE)
        
    with open(filepath, 'w') as f:
        f.write(new_content)

def run_pseudo_select():
    # Prefer scf.in if it exists, otherwise use template scf.dat
    target_file = "scf.in" if os.path.exists("scf.in") else "scf.dat"
    
    if not os.path.exists(target_file):
        print(f"\n[!] Error: Neither 'scf.in' nor 'scf.dat' found. Run Module 100 first.")
        return

    print(f"\n--- [101] Smart Pseudopotential Selection ({target_file}) ---")

    # 1. Local Verification Phase (The "See" Logic)
    required_upfs = extract_required_upfs(target_file)
    
    if not required_upfs:
        print("[!] Warning: No pseudopotentials found in ATOMIC_SPECIES block. Proceeding to library...")
    else:
        print("[*] Scanning local directory for specified pseudopotentials...")
        all_local = True
        for upf in required_upfs:
            if os.path.exists(upf):
                print(f"  [✓] Found locally: {upf}")
            else:
                print(f"  [!] Missing locally: {upf}")
                all_local = False
        
        if all_local:
            print("\n[+] All required pseudopotentials exist in the current directory.")
            print(f"[*] Updating {target_file} to use pseudo_dir='./'")
            update_pseudo_dir(target_file, "./")
            return

    # 2. Global Library Fetch Phase (Triggered only if local files are missing)
    print("\n[!] Local pseudopotentials incomplete. Accessing global library...")
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

    # 3. Update target_file ATOMIC_SPECIES with smart mapping
    with open(target_file, 'r') as f:
        lines = f.readlines()

    new_lines = []
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
            if not parts or any(card in line for card in ["ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS"]):
                in_species_block = False
                new_lines.append(line)
            else:
                element = parts[0]
                mass = parts[1]
                
                # Smart Matching: Look for file starting with 'Element.'
                match = next((upf for upf in available_upfs if upf.startswith(f"{element}.")), None)
                
                if match:
                    new_lines.append(f"  {element:<5} {mass:<10} {match}\n")
                    print(f"    Mapped: {element} -> {match}")
                else:
                    print(f"    [!] Warning: No pseudo found for {element} in library.")
                    new_lines.append(line)
        else:
            new_lines.append(line)

    with open(target_file, 'w') as f:
        f.writelines(new_lines)
        
    # 4. Update pseudo_dir tag to the absolute path of the selected library
    update_pseudo_dir(target_file, target_dir)

    print(f"\n[+] Success! {target_file} updated with specific library filenames and path.")

if __name__ == "__main__":
    run_pseudo_select()

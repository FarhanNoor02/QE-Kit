import os
import glob
import re
import subprocess

def run_302_pdos_processing():
    print("\n--- [302] PDOS Post-Processing: sumpdos.x Automation ---")
    
    # 1. Get Prefix from User
    prefix = input("Enter the calculation prefix (e.g., nfb): ").strip()
    if not prefix:
        print("[!] Prefix cannot be empty."); return

    # 2. Identify PDOS files using the QE naming convention:
    # prefix.pdos.pdos_atm#ID(Atom)_wfc#ID(orbital)
    pattern = f"{prefix}.pdos.pdos_atm#*"
    pdos_files = glob.glob(pattern)

    if not pdos_files:
        print(f"[!] No PDOS files found matching pattern: {pattern}")
        return

    # 3. Parse files to detect unique Atoms and their Orbitals
    # Regex to extract Atom name and Orbital letter from the QE filename
    # Example: ...atm#1(Nb)_wfc#5(d) -> Group 1: Nb, Group 2: d
    parser = re.compile(r"atm#\d+\((\w+)\)_wfc#\d+\((\w)\)")
    
    atoms_data = {} # Structure: {'Nb': {'s', 'p', 'd'}, 'Fe': {'s', 'p', 'd'}}

    for f in pdos_files:
        match = parser.search(f)
        if match:
            atom = match.group(1)
            orbital = match.group(2)
            
            if atom not in atoms_data:
                atoms_data[atom] = set()
            atoms_data[atom].add(orbital)

    if not atoms_data:
        print("[!] Could not parse atom/orbital information from filenames.")
        return

    # 4. Print detected summary
    print("\n[v] Detected Species and Orbitals:")
    for atom, orbitals in atoms_data.items():
        print(f"    - {atom}: {', '.join(sorted(orbitals))}")

    # 5. Execute sumpdos.x commands
    print("\n--- Running sumpdos.x ---")
    
    for atom, orbitals in atoms_data.items():
        # A. Sum by Orbital for this Atom (e.g., Nb_s.dat, Nb_p.dat)
        for orb in orbitals:
            output_file = f"{atom}_{orb}.dat"
            # Command: sumpdos.x *prefix.pdos.pdos_atm#*(Atom)_wfc#*(orb)* > output.dat
            # We use the wildcard inside the command string for shell expansion
            cmd = f"sumpdos.x *({atom})*({orb}) > {output_file}"
            
            print(f" [+] Generating: {output_file}")
            subprocess.run(cmd, shell=True)

        # B. Total Sum for this Atom (e.g., Nb.dat)
        total_output = f"{atom}.dat"
        cmd_total = f"sumpdos.x *({atom})* > {total_output}"
        
        print(f" [+] Generating Total: {total_output}")
        subprocess.run(cmd_total, shell=True)

    print("\n[v] PDOS processing complete. Files generated in the current directory.")

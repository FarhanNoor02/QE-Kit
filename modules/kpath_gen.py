# modules/kpath_gen.py
import os
import seekpath
import numpy as np
import spglib
import questionary

def run_kpath_gen():
    dat_file = "scf.dat"
    if not os.path.exists(dat_file):
        print(f"Error: {dat_file} not found. Run Module 100 first.")
        return

    # 1. User Selection for K-Point Type
    k_type = questionary.select(
        "Select K-Point Type:",
        qmark="-->",
        choices=[
            {"name": "A) Uniform (Automatic grid)", "value": "A"},
            {"name": "B) Homogeneous (Explicit crystal grid)", "value": "B"},
            {"name": "C) Band (High-symmetry path)", "value": "C"}
        ]
    ).ask()

    k_output = ""

    if k_type == "A":
        grid = questionary.text("Enter grid (e.g., 4 4 4):", default="4 4 4").ask()
        k_output = f"K_POINTS {{automatic}}\n {grid} 0 0 0\n"

    elif k_type == "B":
        grid_str = questionary.text("Enter homogeneous grid (e.g., 6 6 6):", default="6 6 6").ask()
        try:
            nx, ny, nz = map(int, grid_str.split())
            total_points = nx * ny * nz
            weight = 1.0 / total_points
            k_output = f"K_POINTS crystal\n{total_points}\n"
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        x, y, z = i/nx, j/ny, k/nz
                        k_output += f"  {x:12.8f} {y:12.8f} {z:12.8f}  {weight:12.8e}\n"
            print(f"[+] Generated {total_points} explicit points.")
        except ValueError:
            print("Invalid grid format.")
            return

    elif k_type == "C":
        structure = parse_dat_for_seekpath(dat_file)
        if structure:
            try:
                # FIX: Pass symprec directly to seekpath and feed it the raw structure.
                # Seekpath handles the primitive reduction internally safely.
                result = seekpath.get_path(structure, with_time_reversal=True, symprec=1e-3)
                k_output = format_seekpath_to_qe(result)
                
                orig_atoms = len(structure[2])
                prim_atoms = len(result['primitive_types'])
                
                print(f"[+] Space Group identified: {result['spacegroup_international']} ({result['spacegroup_number']})")
                if orig_atoms != prim_atoms:
                    print(f"[+] Primitive cell found. Atoms reduced from {orig_atoms} to {prim_atoms}")
                    print("[!] WARNING: Seekpath generated K-points based on the primitive cell.")
                    print("    Ensure your scf.dat matches this primitive structure for band calculations.")
            except Exception as e:
                print(f"[!] Seekpath failed: {e}")
                print("[*] Try running '300: Refine Symmetry' to clean up coordinates first.")
                return

    # 2. WIPE AND REPLACE LOGIC
    if k_output:
        update_scf_dat_with_kpoints(dat_file, k_output)
    else:
        print("\n[!] No K-Points were added.")


def update_scf_dat_with_kpoints(filename, new_kpoints):
    """Reads scf.dat, removes any existing K_POINTS block, and adds the new one."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    final_lines = []
    skip_mode = False

    for line in lines:
        if "K_POINTS" in line:
            skip_mode = True
            continue
        
        if skip_mode:
            if line.strip() == "" or any(card in line for card in ["ATOMIC_", "CELL_"]):
                skip_mode = False
                if any(card in line for card in ["ATOMIC_", "CELL_"]):
                    final_lines.append(line)
            continue
        
        final_lines.append(line)

    with open(filename, 'w') as f:
        f.writelines(final_lines)
        f.write("\n" + new_kpoints)
    
    print(f"\n[+] Success! {filename} updated. Old K-points removed and replaced.")


def parse_dat_for_seekpath(filename):
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        cell, positions, numbers = [], [], []
        species_map, next_num = {}, 1
        in_cell, in_pos = False, False
        for line in lines:
            if "CELL_PARAMETERS" in line: in_cell = True; in_pos = False; continue
            if "ATOMIC_POSITIONS" in line: in_pos, in_cell = True, False; continue
            if in_cell and len(cell) < 3:
                parts = line.split()
                if len(parts) == 3: cell.append([float(x) for x in parts])
            elif in_pos:
                parts = line.split()
                if len(parts) < 4: continue
                if parts[0] not in species_map:
                    species_map[parts[0]] = next_num
                    next_num += 1
                numbers.append(species_map[parts[0]])
                positions.append([float(x) for x in parts[1:4]])
        
        raw_structure = (np.array(cell), np.array(positions), numbers)
        
        # FIX: Return raw_structure directly. 
        # Pre-processing with spglib.find_primitive caused numerical collisions inside seekpath.
        return raw_structure

    except Exception as e:
        print(f"Parsing error: {e}"); return None


def format_seekpath_to_qe(result):
    path = result['path']
    coords = result['point_coords']
    output = f"K_POINTS {{crystal_b}}\n{len(path) + 1}\n"
    first_label = path[0][0]
    p = coords[first_label]
    output += f"  {p[0]:.4f} {p[1]:.4f} {p[2]:.4f} 20 ! {first_label}\n"
    for segment in path:
        label = segment[1]
        p = coords[label]
        output += f"  {p[0]:.4f} {p[1]:.4f} {p[2]:.4f} 20 ! {label}\n"
    return output

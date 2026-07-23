import os
import seekpath
import numpy as np
import spglib
import questionary
import re

def run_kpath_gen():
    dat_file = "scf.dat"
    if not os.path.exists(dat_file):
        print(f"[!] Error: {dat_file} not found. Run Module 100 first.")
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
            print("[!] Invalid grid format.")
            return

    elif k_type == "C":
        structure = parse_dat_for_seekpath(dat_file)
        if structure:
            try:
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
                print(f"[!] Seekpath failed:")
                print(f"    {e}\n")
                print("[*] Possible causes:")
                print("    • malformed CELL_PARAMETERS")
                print("    • incorrect ATOMIC_POSITIONS")
                print("    • atom count inconsistent with nat")
                print("    • symmetry tolerance too strict")
                print("    • structure is not a valid primitive/conventional cell")
                return

    # 2. WIPE AND REPLACE LOGIC
    if k_output:
        update_scf_dat_with_kpoints(dat_file, k_output)
    else:
        print("\n[!] No K-Points were added.")


def update_scf_dat_with_kpoints(filename, new_kpoints):
    """Safely replaces the K_POINTS block without deleting trailing cards."""
    with open(filename) as f:
        lines = f.readlines()

    output = []
    skipping = False

    qe_cards = (
        "ATOMIC_SPECIES",
        "ATOMIC_POSITIONS",
        "CELL_PARAMETERS",
        "K_POINTS",
        "CONSTRAINTS",
        "OCCUPATIONS",
        "HUBBARD",
        "&",
    )

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("K_POINTS"):
            skipping = True
            continue

        if skipping:
            if stripped == "" or any(stripped.startswith(card) for card in qe_cards):
                skipping = False
                if stripped != "":
                    output.append(line)
            continue

        output.append(line)

    with open(filename, "w") as f:
        f.writelines(output)
        if not output[-1].endswith("\n"):
            f.write("\n")
        f.write("\n")
        f.write(new_kpoints)

    print(f"\n[+] Success! {filename} updated.")


def parse_dat_for_seekpath(filename):
    """
    Parse scf.dat into the (cell, positions, numbers) tuple expected by Seekpath.
    Robust against trailing K-point grids and handles both ibrav and CELL_PARAMETERS.
    """
    try:
        with open(filename, "r") as f:
            lines = f.readlines()

        cell = []
        positions = []
        numbers = []
        species_map = {}
        next_num = 1

        nat = None
        ibrav = 0
        c_over_a = 1.0

        in_system = False
        in_cell = False
        in_pos = False

        qe_cards = (
            "K_POINTS",
            "CELL_PARAMETERS",
            "ATOMIC_SPECIES",
            "CONSTRAINTS",
            "OCCUPATIONS",
            "HUBBARD",
        )

        for line in lines:
            stripped = line.strip()

            # ---------------------------
            # Read nat and ibrav from &SYSTEM
            # ---------------------------
            if stripped.startswith("&SYSTEM"):
                in_system = True
                continue

            if in_system:
                if stripped == "/":
                    in_system = False
                    continue

                # Extract variables safely (handles inline commas)
                if "nat" in stripped:
                    match = re.search(r"nat\s*=\s*(\d+)", stripped)
                    if match: nat = int(match.group(1))
                if "ibrav" in stripped:
                    match = re.search(r"ibrav\s*=\s*(-?\d+)", stripped)
                    if match: ibrav = int(match.group(1))
                if "celldm(3)" in stripped:
                    match = re.search(r"celldm\(3\)\s*=\s*([\d\.]+)", stripped)
                    if match: c_over_a = float(match.group(1))

            # ---------------------------
            # CELL_PARAMETERS
            # ---------------------------
            if stripped.startswith("CELL_PARAMETERS"):
                in_cell = True
                in_pos = False
                continue

            if in_cell:
                if len(cell) < 3:
                    parts = stripped.split()
                    if len(parts) != 3:
                        raise ValueError("CELL_PARAMETERS must contain exactly 3 vectors.")
                    cell.append([float(x) for x in parts])
                    if len(cell) == 3:
                        in_cell = False
                continue

            # ---------------------------
            # ATOMIC_POSITIONS
            # ---------------------------
            if stripped.startswith("ATOMIC_POSITIONS"):
                in_pos = True
                continue

            if in_pos:
                if (
                    stripped == ""
                    or stripped.startswith("&")
                    or any(stripped.startswith(card) for card in qe_cards)
                ):
                    in_pos = False
                    continue  # Stop parsing atoms safely

                parts = stripped.split()
                if len(parts) < 4:
                    continue

                if parts[0] not in species_map:
                    species_map[parts[0]] = next_num
                    next_num += 1

                numbers.append(species_map[parts[0]])
                positions.append([float(x) for x in parts[1:4]])

        # ---------------------------
        # Handle ibrav cell generation if CELL_PARAMETERS was missing
        # ---------------------------
        if len(cell) == 0 and ibrav > 0:
            if ibrav == 1:   cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            elif ibrav == 2: cell = [[-0.5, 0.0, 0.5], [0.0, 0.5, 0.5], [-0.5, 0.5, 0.0]]
            elif ibrav == 3: cell = [[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5]]
            elif ibrav == 4: cell = [[1.0, 0.0, 0.0], [-0.5, 0.86602540378, 0.0], [0.0, 0.0, c_over_a]]
            elif ibrav == 6: cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, c_over_a]]
            else:
                raise ValueError(f"ibrav={ibrav} detected but not supported without explicit CELL_PARAMETERS.")

        # ---------------------------
        # Sanity checks
        # ---------------------------
        if len(cell) != 3:
            raise ValueError("Failed to read or construct 3x3 CELL_PARAMETERS.")

        if nat is not None and len(positions) != nat:
            raise ValueError(f"Expected {nat} atoms but parsed {len(positions)} atoms.")

        return (np.array(cell, dtype=float), np.array(positions, dtype=float), numbers)

    except Exception as e:
        print(f"[!] Parsing error: {e}")
        return None


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

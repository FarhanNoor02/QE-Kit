import os
import subprocess
import re

def run_300_structure_refinement():
    print("\n--- [300] Structural Refinement: CELL -> IBRAV/CELLDM ---")
    
    relax_out = "relax.out"
    scf_in = "scf.in"
    # Note: We use scf.dat as the 'template' for species/kpoints if available, 
    # but for this workflow, we surgically update scf.in directly.
    
    if not os.path.exists(relax_out) or not os.path.exists(scf_in):
        print("[!] Error: Necessary files (relax.out or scf.in) missing."); return

    # 1. Extract the Final Coordinates block from relax.out
    with open(relax_out, 'r') as f:
        lines = f.readlines()

    final_block = []
    found = False
    for line in lines:
        if "Begin final coordinates" in line:
            found = True
            final_block = [] 
            continue
        if found:
            if "End final coordinates" in line:
                found = False
                break
            final_block.append(line)

    if not final_block:
        print("[!] Final coordinates not found. Check convergence in relax.out."); return

    # 2. Extract Data from final_block (Positions and Cell)
    pos_lines = []
    cell_matrix_lines = []
    is_pos = False
    is_cell = False

    for line in final_block:
        if "ATOMIC_POSITIONS" in line:
            pos_lines.append(line)
            is_pos = True; is_cell = False
        elif "CELL_PARAMETERS" in line:
            is_cell = True; is_pos = False
        elif is_pos:
            pos_lines.append(line)
        elif is_cell:
            if len(line.split()) == 3:
                cell_matrix_lines.append(line)

    # 3. Find 'alat' for cell2ibrav.x
    alat = None
    for line in reversed(lines):
        if "lattice parameter (alat)" in line:
            alat = line.split()[4]
            break

    # 4. Invoke cell2ibrav.x
    with open("c2i.in", "w") as f:
        f.writelines(cell_matrix_lines)
        f.write(f"{alat}\n")

    try:
        process = subprocess.run(["cell2ibrav.x"], input=open("c2i.in").read(), 
                                 capture_output=True, text=True)
        c2i_output = process.stdout
    except FileNotFoundError:
        print("[!] Error: cell2ibrav.x not in PATH."); return

    # 5. Parse ibrav and ALL celldm
    new_vals = {}
    ibrav_match = re.search(r'ibrav\s*=\s*(\d+)', c2i_output)
    if ibrav_match: new_vals['ibrav'] = ibrav_match.group(1)
    celldm_matches = re.findall(r'celldm\((\d+)\)\s*=\s*([\d\.]+)', c2i_output)
    for idx, val in celldm_matches:
        new_vals[f'celldm({idx})'] = val

    # 6. Surgical scf.in rewrite
    with open(scf_in, 'r') as f:
        old_scf = f.readlines()

    new_scf = []
    in_system = False
    positions_injected = False
    skip_mode = None 
    
    lattice_vars = re.compile(r'^\s*(ibrav|A|B|C|cosbc|cosac|cosab|celldm\(\d\))\s*=', re.IGNORECASE)

    for line in old_scf:
        clean_line = line.strip().upper()

        # Handle &SYSTEM namelist
        if "&SYSTEM" in clean_line:
            in_system = True
            new_scf.append(line)
            # Inject new symmetry parameters
            param_line = f"  ibrav = {new_vals.get('ibrav', '0')}"
            for k, v in sorted(new_vals.items()):
                if 'celldm' in k: param_line += f", {k} = {v}"
            new_scf.append(param_line + ",\n")
            continue
        
        if in_system:
            if "/" in clean_line: in_system = False
            if lattice_vars.match(line): continue

        # Handle Block Removal and Replacement
        if "ATOMIC_SPECIES" in clean_line:
            skip_mode = 'species'
            continue
        if "CELL_PARAMETERS" in clean_line:
            skip_mode = 'cell'
            continue
        if "ATOMIC_POSITIONS" in clean_line:
            if not positions_injected:
                new_scf.extend(pos_lines)
                positions_injected = True
            skip_mode = 'positions'
            continue
        
        # Reset skip_mode when hitting K_POINTS or next namelist
        if any(x in clean_line for x in ["K_POINTS", "&ELECTRONS", "&IONS", "&CONTROL"]):
            skip_mode = None

        # Line exclusion logic
        if skip_mode == 'species' and len(line.split()) >= 3: continue
        if skip_mode == 'cell' and len(line.split()) == 3: continue
        if skip_mode == 'positions' and len(line.split()) >= 4: continue
        
        if skip_mode is None:
            new_scf.append(line)

    with open("scf_refined.in", "w") as f:
        f.writelines(new_scf)
    
    print(f"\n[+] Refinement Complete: 'scf_refined.in' generated.")
    print(f"    - Standardized to ibrav {new_vals.get('ibrav')}")
    print(f"    - Updated ATOMIC_POSITIONS from relax.out")
    print(f"    - Preserved K_POINTS and electronic settings.")

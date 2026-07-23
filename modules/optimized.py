# modules/optimized.py
import os
import subprocess
import re

def run_300_structure_refinement():
    print("\n--- [300] Structural Refinement: CELL -> IBRAV/CELLDM ---")
    
    relax_out = "relax.out"
    scf_in = "scf.in"
    dat_file = "scf.dat"
    
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

    lattice_vars = re.compile(r'^\s*(ibrav|A|B|C|cosbc|cosac|cosab|celldm\(\d\))\s*=', re.IGNORECASE)

    # 6. Surgical scf.in rewrite
    with open(scf_in, 'r') as f:
        old_scf = f.readlines()

    new_scf = []
    in_system = False
    positions_injected = False
    skip_mode = None 

    for line in old_scf:
        clean_line = line.strip().upper()

        if "&SYSTEM" in clean_line:
            in_system = True
            new_scf.append(line)
            param_line = f"  ibrav = {new_vals.get('ibrav', '0')}"
            for k, v in sorted(new_vals.items()):
                if 'celldm' in k: param_line += f", {k} = {v}"
            new_scf.append(param_line + ",\n")
            continue
        
        if in_system:
            if "/" in clean_line: in_system = False
            if lattice_vars.match(line): continue

        if "ATOMIC_SPECIES" in clean_line:
            skip_mode = None
            new_scf.append(line)
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
        
        if any(x in clean_line for x in ["K_POINTS", "&ELECTRONS", "&IONS", "&CONTROL", "ATOMIC_SPECIES"]):
            skip_mode = None

        if skip_mode == 'cell' and len(line.split()) == 3: continue
        if skip_mode == 'positions' and len(line.split()) >= 4: continue
        
        if skip_mode is None:
            new_scf.append(line)

    with open("scf_refined.in", "w") as f:
        f.writelines(new_scf)

    # 7. Sync the Master Template (scf.dat) for downstream Zone 2 compatibility
    if os.path.exists(dat_file):
        with open(dat_file, 'r') as f:
            old_dat = f.readlines()

        new_dat = []
        in_system_dat = False
        positions_injected_dat = False
        skip_mode_dat = None

        for line in old_dat:
            clean_line = line.strip().upper()

            # Preserve # PSEUDO_DIR_PATH headers
            if line.strip().startswith("#"):
                new_dat.append(line)
                continue

            if "&SYSTEM" in clean_line:
                in_system_dat = True
                new_dat.append(line)
                param_line = f"  ibrav = {new_vals.get('ibrav', '0')}"
                for k, v in sorted(new_vals.items()):
                    if 'celldm' in k: param_line += f", {k} = {v}"
                new_dat.append(param_line + ",\n")
                continue
            
            if in_system_dat:
                if "/" in clean_line: in_system_dat = False
                if lattice_vars.match(line): continue

            if "ATOMIC_SPECIES" in clean_line:
                skip_mode_dat = None
                new_dat.append(line)
                continue

            if "CELL_PARAMETERS" in clean_line:
                skip_mode_dat = 'cell'
                continue
                
            if "ATOMIC_POSITIONS" in clean_line:
                if not positions_injected_dat:
                    new_dat.extend(pos_lines)
                    positions_injected_dat = True
                skip_mode_dat = 'positions'
                continue
            
            if any(x in clean_line for x in ["K_POINTS", "ATOMIC_SPECIES", "&"]):
                if not clean_line.startswith("&SYSTEM"):
                    skip_mode_dat = None

            if skip_mode_dat == 'cell' and len(line.split()) == 3: continue
            if skip_mode_dat == 'positions' and len(line.split()) >= 4: continue
            
            if skip_mode_dat is None:
                new_dat.append(line)

        with open(dat_file, "w") as f:
            f.writelines(new_dat)
    
    print(f"\n[+] Refinement Complete: 'scf_refined.in' generated.")
    print(f"    - Standardized to ibrav {new_vals.get('ibrav')}")
    print(f"    - Updated ATOMIC_POSITIONS from relax.out")
    print(f"    - Master template 'scf.dat' synced with refined structure.")

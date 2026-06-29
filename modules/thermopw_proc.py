# modules/thermopw_proc.py
import os
import re
from utils import check_strain, symmetric_mat

def run_305_thermopw_proc():
    print("\n--- [Zone 3] Thermo_PW Post-Processor ---")
    scf_out = "scf.out"
    result_file = "result.txt"

    if not os.path.exists(scf_out):
        print("[!] Error: scf.out not found.")
        return

    with open(result.txt, "w") as res:
        res.write("--- Thermo_PW Analysis Results ---\n\n")

        # 1. Parse VRH and Debye
        with open(scf_out, "r") as f:
            content = f.read()
            vrh = re.search(r"Voigt-Reuss-Hill average of the two approximations:(.*?)\n\s*\n", content, re.DOTALL)
            if vrh:
                res.write("Voigt-Reuss-Hill Averages:\n" + vrh.group(1).strip() + "\n\n")
            
            debye = re.search(r"Debye temperature =\s*([\d\.]+)\s*K", content)
            if debye:
                res.write(f"Debye Temperature: {debye.group(1)} K\n\n")

        # 2. Add Strain Tensors
        res.write("--- Strain Tensors ---\n")
        for idx, s in enumerate(check_strain.parse_strains(scf_out), 1):
            res.write(f"Strain {idx}:\n{s}\n\n")

        # 3. Add Symmetric Elastic Matrix
        res.write("--- Symmetric Elastic Constant Matrix ---\n")
        # Ensure path points to the correct file in elastic_constants
        if os.path.exists("elastic_constants/output_el_cons.dat.g1"):
            mat_str = symmetric_mat.get_elastic_matrix_string()
            res.write(mat_str)
        else:
            res.write("Elastic constants file not found in elastic_constants/")

    print(f"[+] Post-processing complete. Summary saved to {result_file}")

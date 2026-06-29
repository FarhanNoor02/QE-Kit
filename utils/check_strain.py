import numpy as np
import os

def parse_strains(filename="scf.out"):
    strains = []
    if not os.path.exists(filename): return strains
    with open(filename, "r") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if "Applying the following strain" in lines[i]:
            strain_matrix = []
            for j in range(1, 4):
                values = lines[i+j].strip().strip("()").split(",")
                strain_matrix.append([float(v) for v in values])
            strains.append(np.array(strain_matrix))
            i += 4
        else:
            i += 1
    return strains

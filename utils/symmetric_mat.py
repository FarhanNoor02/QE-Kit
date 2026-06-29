import numpy as np

def get_elastic_matrix_string(filename="elastic_constants/output_el_cons.dat.g1"):
    with open(filename, "r") as f:
        numbers = [float(x) for x in f.read().split()]
    
    C_raw = np.array(numbers[:36]).reshape((6, 6))
    C_sym = 0.5 * (C_raw + C_raw.T) / 10.0 # Convert kbar to GPa
    
    output = "Elastic Stiffness Matrix (GPa):\n"
    for row in C_sym:
        output += "  ".join(f"{val:12.7f}" for val in row) + "\n"
    return output

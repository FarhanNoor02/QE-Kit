# modules/pdos_gen.py
import os
import questionary
from utils import qe_xml_parser

def run_pdos_gen():
    print("\n--- [203] Projected DOS (PDOS) Generation ---")
    
    # 1. Verification Step
    prefix = questionary.text("Enter calculation prefix:", default="ncn").ask()
    calc_type, status = qe_xml_parser.check_calc_status(prefix)

    if status == "Not Found":
        print(f"\n[!] Error: No XML data found for prefix '{prefix}'.")
        print("    You must run an NSCF calculation before generating PDOS.")
        return
    
    if calc_type != "nscf":
        print(f"\n[!] Warning: Last calculation was '{calc_type.upper()}', not 'NSCF'.")
        print("    PDOS requires a fixed charge density run (NSCF).")
        if not questionary.confirm("Do you want to proceed anyway?").ask():
            return

    # 2. PDOS Specific User Inputs
    filpdos = questionary.text("Enter output PDOS filename (filpdos):", default=f"{prefix}.pdos").ask()
    degauss = questionary.text("Smearing width (degauss, Ry):", default="0.01").ask()
    delta_e = questionary.text("Energy step (DeltaE, eV):", default="0.005").ask()

    # 3. Assemble pdos.in for projwfc.x
    pdos_template = f"""&projwfc
    prefix = '{prefix}',
    outdir = './outdir/',
    ngauss       = 0,
    degauss      = {degauss},
    DeltaE       = {delta_e},
    filpdos = '{filpdos}'
/
"""

    # 4. Write the file
    with open("pdos.in", "w") as f:
        f.write(pdos_template)

    print("\n" + "-"*50)
    print("Projected DOS input file ready to be run with projwfc.x")
    print("-" * 50)

if __name__ == "__main__":
    run_pdos_gen()

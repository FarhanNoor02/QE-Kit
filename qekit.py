import os
import questionary  
import sys
from questionary import Style, Separator
from modules import (
    structure2in, pseudo_select, kpath_gen, scf_gen, vcrelax_gen, 
    nscf_gen, pdos_gen, bands_gen, optical_gen, phonon_gen, 
    optimized, phonon_proc, optical_proc, pdos_proc
)

from utils import config_manager

# --- Stylish Layout ---
custom_style = Style([
    ('separator', 'fg:#ffcc00 bold'),       # Gold for Zone Dividers
    ('qmark', 'fg:#00ff00 bold'),          # Green for prompt markers
    ('question', 'bold'),                  # Bold text for questions
    ('pointer', 'fg:#00ff00 bold'),        # Green selection arrow
    ('highlighted', 'fg:#00ff00 bold'),    # Green text when highlighted
    ('selected', 'fg:#ccff00'),            # Neon Yellow for final choice
    ('instruction', 'fg:#888888 italic'),  # Dim help text
])

def show_banner():
    # os.system('clear' if os.name == 'posix' else 'cls') # Optional: Clear terminal on start
    banner = r"""
    ╔════════════════════════════════════════════════════════════════╗
    ║   ____    _____      _      ______  _______                    ║
    ║  / __ \  |____|     | |/ /  |_  _| |__  __|                    ║
    ║ | |  | | | |__      | ' /    | |      | |                      ║
    ║ | |  | | |  __| --  | <      | |      | |                      ║
    ║ | |__| | | |___     | . \   _| |_     | |                      ║
    ║  \ __\\  |_____|    |_|\_\ |_____|    |_|                      ║
    ║ A Pre- & Post-Processing Suite for Quantum ESPRESSO            ║
    ║                                                                ║
    ║ [ By: Farhan Noor, University of Dhaka ]                       ║
    ╚════════════════════════════════════════════════════════════════╝
    """
    print(banner)

def pre_processing_menu():
    while True:
        choice = questionary.select(
            "--- ZONE 1 | STRUCTURAL DISCOVERY ---",
            style=custom_style,
            choices=[
                "100: Structure to Input",
                "101: Pseudopotential Selection",
                "102: K-Point Path Generation",
                Separator(""),
                "<- Return to Main Menu"
            ]
        ).ask()

        if choice == "100: Structure to Input": structure2in.run_structure2in()
        elif choice == "101: Pseudopotential Selection": pseudo_select.run_pseudo_select()
        elif choice == "102: K-Point Path Generation": kpath_gen.run_kpath_gen()
        elif choice == "<- Return to Main Menu": break

def input_generation_menu():
    while True:
        choice = questionary.select(
            "--- ZONE 2 | INPUT FILE ARCHITECT ---",
            style=custom_style,
            choices=[
                "200: SCF Setup", "201: VC-Relax Setup", "202: NSCF Setup",
                "203: DOS/PDOS Setup", "204: Bands Setup", "205: Optical Setup",
                "206: Phonon Setup", Separator(""), "<- Return to Main Menu"
            ]
        ).ask()
        
        if "200" in choice: scf_gen.run_scf_gen()
        elif "201" in choice: vcrelax_gen.run_vcrelax_gen()
        elif "202" in choice: nscf_gen.run_nscf_gen()
        elif "203" in choice: pdos_gen.run_pdos_gen()
        elif "204" in choice: bands_gen.run_bands_gen()
        elif "205" in choice: optical_gen.run_optical_gen()
        elif "206" in choice: phonon_gen.run_phonon_gen()
        elif "<- Return to Main Menu" in choice: break

def post_processing_menu():
    while True:
        choice = questionary.select(
            "--- ZONE 3 | PHYSICS & ANALYSIS ---",
            style=custom_style,
            choices=[
                "300: Refine Symmetry (cell2ibrav)", "301: Phonon Dispersion Setup",
                "302: PDOS Summation", "303: Optical Constant Derivation",
                Separator(""), "<- Return to Main Menu"
            ]
        ).ask()
        
        if "300" in choice: optimized.run_300_structure_refinement()
        elif "301" in choice: phonon_proc.run_301_phonon_processing()
        elif "302" in choice: pdos_proc.run_302_pdos_summation()
        elif "303" in choice: optical_proc.run_303_optical_processing()
        elif "<- Return to Main Menu" in choice: break

def automation_visualization_menu():
    utils_dir = os.path.join(os.path.dirname(__file__), "utils")
    while True:
        choice = questionary.select(
            "--- ZONE 4 | WORKFLOW & VISUALIZATION ---",
            style=custom_style,
            choices=[
                "401: Run HT-Phonon Pipeline (Bash Driver)",
                "402: Overlaid Plotter (Gnuplot)",
                "403: Subplot/Spectral Plotter (Gnuplot)",
                Separator(""),
                "<- Return to Main Menu"
            ]
        ).ask()

        if "401" in choice: os.system(f"bash {os.path.join(utils_dir, 'ht_phonon.sh')}")
        elif "402" in choice: os.system(f"bash {os.path.join(utils_dir, 'plot-overlay.sh')}")
        elif "403" in choice: os.system(f"bash {os.path.join(utils_dir, 'plot-subplots.sh')}")
        elif "<- Return to Main Menu" in choice: break

def main():
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg == "--300": optimized.run_300_structure_refinement()
        elif arg == "--206": phonon_gen.run_phonon_gen(automated=True)
        elif arg == "--301": phonon_proc.run_301_phonon_processing()
        sys.exit(0)

    while True:
        show_banner()
        category = questionary.select(
            "Main Menu » Select Workstream",
            style=custom_style,
            choices=[
                Separator("========================================="),
                " ▶ ZONE 1 | Structural Discovery",
                " ▶ ZONE 2 | Input File Architect",
                " ▶ ZONE 3 | Physics & Analysis",
                " ▶ ZONE 4 | Workflow & Visualization",
                Separator("========================================="),
                " ⚙ Settings",
                " ✖ Exit Application"
            ]
        ).ask()

        if "ZONE 1" in category: pre_processing_menu()
        elif "ZONE 2" in category: input_generation_menu()
        elif "ZONE 3" in category: post_processing_menu()
        elif "ZONE 4" in category: automation_visualization_menu()
        elif "Settings" in category:
            current_path = config_manager.get_pseudo_dir()
            print(f"\n[Current Pseudo Directory]: {current_path}")
            new_path = questionary.text("Enter new absolute path:").ask()
            if new_path and os.path.isdir(new_path):
                config_manager.save_pseudo_dir(new_path)
                print("[+] Path updated.")
        elif "Exit" in category:
            print("\nExiting QE-Kit. Happy Computing!")
            sys.exit()

if __name__ == "__main__":
    main()

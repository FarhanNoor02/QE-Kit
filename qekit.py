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

# --- Professional Terminal Styling ---
# Using Green (Success), Gold (Headers), and Neon Yellow (Selections)
custom_style = Style([
    ('separator', 'fg:#ffcc00 bold'),       # Gold for Zone Headers
    ('qmark', 'fg:#00ff00 bold'),          # Green for prompt markers
    ('question', 'bold'),                  # Bold text for questions
    ('pointer', 'fg:#00ff00 bold'),        # Green selection arrow
    ('highlighted', 'fg:#00ff00 bold'),    # Green text when highlighted
    ('selected', 'fg:#ccff00'),            # Yellow for selected text
    ('instruction', 'fg:#888888 italic'),  # Grey for help text
])

def show_banner():
    # os.system('clear' if os.name == 'posix' else 'cls') # Optional: Clear terminal on start
    banner = r"""
    ╔════════════════════════════════════════════════════════════════╗
    ║   ____    _____      _      _____  _______                     ║
    ║  / __ \  |____|     | |/ /  |_  _| |_  __|                     ║
    ║ | |  | | | |__      | ' /    | |     | |                       ║
    ║ | |  | | |  __| --  | <      | |     | |                       ║
    ║ | |__| | | |___     | . \   _| |_    | |                       ║
    ║  \ __\\  |_____|    |_|\_\ |_____|   |_|                       ║
    ║ A Pre- & Post-Processing Suite for Quantum ESPRESSO            ║
    ║   [ Build v0.3 ] | [ University of Dhaka ]                     ║
    ║   [ Lead Dev: Farhan Noor ]                                    ║
    ╚════════════════════════════════════════════════════════════════╝
    """
    print(banner)

def pre_processing_menu():
    while True:
        choice = questionary.select(
            "--- ZONE 1: PRE-PROCESSING ---",
            style=custom_style,
            choices=[
                "100: Structure to Input (CIF/XYZ)",
                "101: Pseudopotential Selection (QE-PotLib)",
                "102: K-Point Path Selection",
                Separator(""),
                "<- Return to Main Menu"
            ]
        ).ask()

        if choice == "100: Structure to Input (CIF/XYZ)":
            structure2in.run_structure2in()
        elif choice == "101: Pseudopotential Selection (QE-PotLib)":
            pseudo_select.run_pseudo_select()
        elif choice == "102: K-Point Path Selection":
            kpath_gen.run_kpath_gen()
        elif choice == "<- Return to Main Menu":
            break

def input_generation_menu():
    while True:
        choice = questionary.select(
            "--- ZONE 2: INPUT GENERATION ---",
            style=custom_style,
            qmark="-->",
            choices=[
                Separator(" [ Primary Calculations ] "),
                "200: Self-Consistent Field (SCF)",
                "201: Variable-Cell Relaxation (VC-relax)",
                "202: Non-Self-Consistent Field (NSCF)",
                Separator(" [ Property Calculations ] "),
                "203: Density of States (DOS)",
                "204: Band Structure",
                "205: Optical Properties (Epsilon)",
                "206: Phonon Calculations",
                Separator(""),
                "<- Return to Main Menu"
            ]
        ).ask()
        
        if choice == "200: Self-Consistent Field (SCF)":
            scf_gen.run_scf_gen()
        elif choice == "201: Variable-Cell Relaxation (VC-relax)":
            vcrelax_gen.run_vcrelax_gen()
        elif choice == "202: Non-Self-Consistent Field (NSCF)":
            nscf_gen.run_nscf_gen()
        elif choice == "203: Density of States (DOS)":
            pdos_gen.run_pdos_gen()
        elif choice == "204: Band Structure":
            bands_gen.run_bands_gen()
        elif choice == "205: Optical Properties (Epsilon)":
            optical_gen.run_optical_gen()
        elif choice == "206: Phonon Calculations":
            phonon_gen.run_phonon_gen()
        elif choice == "<- Return to Main Menu":
            break

def post_processing_menu():
    while True:
        choice = questionary.select(
            "--- ZONE 3: POST-PROCESSING & ANALYSIS ---",
            style=custom_style,
            choices=[
                Separator(" [ Symmetry & Structure ] "),
                "300: Symmetry Refinement (cell2ibrav)",
                Separator(" [ Electronic & Vibrational ] "),
                "301: Phonon Dispersion Setup (q2r/matdyn)",
                "302: PDOS Orbital Summation",
                Separator(" [ Optoelectronics ] "),
                "303: Optical Constants (n, k, R, alpha)",
                Separator(""),
                "<- Return to Main Menu"
            ]
        ).ask()
        
        if choice == "300: Symmetry Refinement (cell2ibrav)":
            optimized.run_300_structure_refinement()
        elif choice == "301: Phonon Dispersion Setup (q2r/matdyn)":
            phonon_proc.run_301_phonon_processing()
        elif choice == "302: PDOS Orbital Summation":
            pdos_proc.run_302_pdos_summation()
        elif choice == "303: Optical Constants (n, k, R, alpha)":
            optical_proc.run_303_optical_processing()
        elif choice == "<- Return to Main Menu":
            break

def main():
    while True:
        show_banner()
        category = questionary.select(
            "RESEARCH WORKFLOW DASHBOARD",
            style=custom_style,
            choices=[
                Separator("  --- 🧪 ZONE 1: DISCOVERY & SETUP ---  "),
                "1. Pre-Processing Suite",
                Separator("  --- ⚙️ ZONE 2: SIMULATION ENGINE ---  "),
                "2. Input File Generation",
                Separator("  --- 📊 ZONE 3: ANALYSIS & PHYSICS --- "),
                "3. Post-Processing Suite",
                Separator("  --- 🛠️ SYSTEM --- "),
                "4. Settings & Paths",
                "Exit Application"
            ]
        ).ask()

        if category == "1. Pre-Processing Suite":
            pre_processing_menu()
        elif category == "2. Input File Generation":
            input_generation_menu()
        elif category == "3. Post-Processing Suite":
            post_processing_menu()
        elif category == "4. Settings & Paths":
            current_path = config_manager.get_pseudo_dir()
            print(f"\n[Current Pseudo Directory]: {current_path}")
            
            new_path = questionary.text("Enter new absolute path for PseudoPotentials:").ask()
            if new_path and os.path.isdir(new_path):
                config_manager.save_pseudo_dir(new_path)
                print("[+] Path updated successfully.")
            else:
                print("[!] Invalid path or cancelled. No changes made.")
        elif category == "Exit Application":
            print("\nExiting QE-Kit. Happy Computing!")
            sys.exit()

if __name__ == "__main__":
    main()

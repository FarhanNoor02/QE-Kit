import os
import questionary  
import sys
from questionary import Style, Separator

from modules import (
    structure2in, pseudo_select, kpath_gen, scf_gen, vcrelax_gen, 
    nscf_gen, pdos_gen, bands_gen, optical_gen, phonon_gen, 
    optimized, phonon_proc, optical_proc, pdos_proc, epc_processor,
    thermopw_gen, thermopw_proc, fs_gen, build_pseudo_lib, bands_proc
)

from utils import (
    config_manager, help_manager, bz_plotter, xrd_plotter,
    check_strain, symmetric_mat, citations
) 

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
    # Clears the terminal screen every time the main menu is loaded
    os.system('clear' if os.name == 'posix' else 'cls')

    banner = r"""
    ╔════════════════════════════════════════════════════════════════╗
    ║   ____    _____     _        _____   _______                   ║
    ║  / __ \  |____|     | |/ /  |__  __| |__  __|                  ║
    ║ | |  | | | |__      | ' /     | |      | |                     ║
    ║ | |  | | |  __| --  | <       | |      | |                     ║
    ║ | |__| | | |___     | . \    _| |_     | |                     ║
    ║  \ __\\  |_____|    |_|\_\  |_____|    |_|                     ║
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
                "206: Phonon Setup", 
                "207: Fermi Surface Setup (fs.x)",
                "208: Thermo_PW Setup (elastic/thermo)",
                Separator(""), "<- Return to Main Menu"
            ]
        ).ask()
        
        if "200" in choice: scf_gen.run_scf_gen()
        elif "201" in choice: vcrelax_gen.run_vcrelax_gen()
        elif "202" in choice: nscf_gen.run_nscf_gen()
        elif "203" in choice: pdos_gen.run_pdos_gen()
        elif "204" in choice: bands_gen.run_bands_gen()
        elif "205" in choice: optical_gen.run_optical_gen()
        elif "206" in choice: phonon_gen.run_phonon_gen()
        elif "207" in choice: fs_gen.run_207_fs_gen()
        elif "208" in choice: thermopw_gen.run_208_thermopw_gen()
        elif "<- Return to Main Menu" in choice: break

def post_processing_menu():
    while True:
        choice = questionary.select(
            "--- ZONE 3 | POST-PROCESSING ---",
            style=custom_style,
            choices=[
                "300: Refine Symmetry (cell2ibrav)", 
                "301: Band Structure Post-Processing", # <-- Inserted here
                "302: Phonon Dispersion Setup",
                "303: PDOS Summation", 
                "304: Optical Constant Derivation", 
                "305: e-ph analysis",
                "306: Thermo_PW Analysis",
                Separator(""), "<- Return to Main Menu"
            ]
        ).ask()
        
        # UI Number maps to the existing Function Name
        if "300" in choice: optimized.run_300_structure_refinement()
        elif "301" in choice: bands_proc.run_306_bands_processing() 
        elif "302" in choice: phonon_proc.run_301_phonon_processing()
        elif "303" in choice: pdos_proc.run_302_pdos_processing()
        elif "304" in choice: optical_proc.run_303_optical_processing()
        elif "305" in choice: epc_processor.run_304_epc_processor()
        elif "306" in choice: thermopw_proc.run_305_thermopw_proc()
        elif "<- Return to Main Menu" in choice: break

def automation_visualization_menu():
    utils_dir = os.path.join(os.path.dirname(__file__), "utils")
    while True:
        choice = questionary.select(
            "--- ZONE 4 | DATA VISUALIZATION ---",
            style=custom_style,
            choices=[
                "401: Run HT-Phonon Pipeline (Bash Driver)",
                "402: Overlaid Plotter (Gnuplot)",
                "403: Subplot/Spectral Plotter (Gnuplot)",
                "404: Interactive Brillouin Zone Visualizer (Matplotlib)",
                "405: Calculated XRD Pattern Plotter (Pymatgen)",
                 Separator(""),
                "<- Return to Main Menu"
            ]
        ).ask()

        if "401" in choice: os.system(f"bash {os.path.join(utils_dir, 'ht_phonon.sh')}")
        elif "402" in choice: os.system(f"bash {os.path.join(utils_dir, 'plot-overlay.sh')}")
        elif "403" in choice: os.system(f"bash {os.path.join(utils_dir, 'plot-subplots.sh')}")
        elif "404" in choice: bz_plotter.run_404_bz_plotter() 
        elif "405" in choice: xrd_plotter.run_405_xrd_plotter()
        elif "<- Return to Main Menu" in choice: break

def main():
    # --- Headless CLI & Help Handling ---
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        
        if arg in ["--help", "-h"]:
            help_manager.display_help()
        elif arg == "--206": 
            phonon_gen.run_phonon_gen(automated=True)
        elif arg == "--207": 
            fs_gen.run_207_fs_gen()
        elif arg == "--208": 
            thermopw_gen.run_208_thermopw_gen()
        elif arg == "--300": 
            optimized.run_300_structure_refinement()
        elif arg == "--301": 
            bands_proc.run_306_bands_processing(automated=True)
        elif arg == "--302": 
            phonon_proc.run_301_phonon_processing()
        # You will need to add these to main() if they aren't there already
        elif arg == "--303": 
            pdos_proc.run_302_pdos_summation()
        elif arg == "--304": 
            optical_proc.run_303_optical_processing()
        elif arg == "--305": 
            epc_processor.run_304_epc_processor()
        elif arg == "--306": 
            thermopw_proc.run_305_thermopw_proc()
        elif arg == "--404": 
            bz_plotter.run_404_bz_plotter()
        elif arg == "--405":
            xrd_plotter.run_405_xrd_plotter()
        else:
            print(f"[!] Unknown argument: {arg}")
            help_manager.display_help()
            
        sys.exit(0)

    # --- Interactive Menu ---
    while True:
        show_banner()
        category = questionary.select(
            "Main Menu » Select Workstream",
            style=custom_style,
            choices=[
                Separator("========================================="),
                " ▶ ZONE 1 | Structural Discovery",
                " ▶ ZONE 2 | Input File Architect",
                " ▶ ZONE 3 | Post-Processor",
                " ▶ ZONE 4 | Visualization",
                Separator("========================================="),
                " ⚙ Settings",
                " 📚 References & Citations",
                " ✖ Exit Application"
            ]
        ).ask()

        if "ZONE 1" in category: pre_processing_menu()
        elif "ZONE 2" in category: input_generation_menu()
        elif "ZONE 3" in category: post_processing_menu()
        elif "ZONE 4" in category: automation_visualization_menu()
        elif "Settings" in category:
            while True:
                current_path = config_manager.get_pseudo_dir() or "Not Set"
                print(f"\n[Current Pseudo Directory]: {current_path}")
                
                setting_choice = questionary.select(
                    "Settings Menu",
                    style=custom_style,
                    choices=[
                        "1: Manually Link Existing Library Path",
                        "2: Build/Download New Global Library (Internet Required)",
                        Separator(""),
                        "<- Back to Main Menu"
                    ]
                ).ask()
                
                if "1:" in setting_choice:
                    new_path = questionary.text("Enter absolute path to your existing library:").ask()
                    if new_path and os.path.isdir(new_path):
                        config_manager.save_pseudo_dir(new_path)
                        print("[+] Path updated successfully.")
                    else:
                        print("[!] Invalid path or directory does not exist.")
                        
                elif "2:" in setting_choice:
                    build_pseudo_lib.build_library()
                    
                elif "<-" in setting_choice:
                    # Clear screen to keep the UI clean when returning to Mainmenu
                    os.system('clear' if os.name == 'posix' else 'cls')
                    break
        elif "References" in category:         
            citations.display_citations()
        elif "Exit" in category:
            # Clears the screen one last time before exiting cleanly
            os.system('clear' if os.name == 'posix' else 'cls')
            print("\nExiting QE-Kit. Happy Computing!\n")
            sys.exit()

if __name__ == "__main__":
    main()

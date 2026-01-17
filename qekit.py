import os
import questionary  # <--- THIS LINE IS MISSING OR MISPLACED
import sys
from modules import structure2in, pseudo_select, kpath_gen, scf_gen,vcrelax_gen, nscf_gen
from modules import pdos_gen, bands_gen, optical_gen, phonon_gen, optimized, phonon_proc
from utils import config_manager

def show_banner():
    print("---------######----########----------##----###-----------------")
    print("-------##-----##---##----------------##--###-----#-----##------")
    print("-------##-----##---#######-----------##-##------##--######-----")
    print("-------##-----##---##--------######--##--##-----##----##-------")
    print("-------###----##---##----------------##---###---##----##-------")
    print("---------######----########----------##-----##--##----###------")
    print("------------####-----------------------------------------------")
    print("=====================By: Farhan Noor===========================")
    print("---A Pre-/Post-processing tool for use with Quantum ESPRESSO---")
    print("---------------------------------------------------------------")

def pre_processing_menu():
    while True:
        choice = questionary.select(
            "--- 1. Pre-Processing ---",
            choices=[
                "100: Structure to Input",
                "101: Pseudopotential Selection",
                "102: K-Point Selection",
                "Return to Main Menu"
            ]
        ).ask()

        if choice == "100: Structure to Input":
            structure2in.run_structure2in()
        elif choice == "101: Pseudopotential Selection":
            pseudo_select.run_pseudo_select()
        elif choice == "102: K-Point Selection":
            kpath_gen.run_kpath_gen()
        elif choice == "Return to Main Menu":
            break # Exit the while loop to go back to main()

# Updated sub-menu in qekit.py

def input_generation_menu():
    while True:
        choice = questionary.select(
            "--- 2. Input File Generation ---",
            qmark="-->",
            choices=[
                "200: SCF",
                "201: VC-relax",
                "202: NSCF",
                "203: Density of States",
                "204: Bands",
                "205: Optical",
                "206: Phonon",
                "Return to Main Menu"
            ]
        ).ask()
        
        if choice == "200: SCF":
            scf_gen.run_scf_gen()
        elif choice == "201: VC-relax":
            vcrelax_gen.run_vcrelax_gen()
        elif choice == "202: NSCF":
            nscf_gen.run_nscf_gen()
        elif choice == "203: Density of States":
            pdos_gen.run_pdos_gen()
        elif choice == "204: Bands":
            bands_gen.run_bands_gen() # To be written
        elif choice == "205: Optical":
            optical_gen.run_optical_gen() # To be written
        elif choice == "206: Phonon":
            phonon_gen.run_phonon_gen()
        elif choice == "Return to Main Menu":
            break
        else:
            print(f"\n[!] Module {choice.split(':')[0]} is currently in development.")
            
def post_processing_menu():
    while True:
        choice = questionary.select(
            "--- 3. Post-Processing ---",
            choices=[
                "300: Optimized structure (Update scf.in)", # Update this text
                "301: Phonon Dispersions (q2r & matdyn setup)",
                "Return to Main Menu"
            ]
        ).ask()
        
        if choice == "300: Optimized structure (Update scf.in)":
            optimized.run_300_structure_refinement() # Call the function
        elif choice == "301: Phonon Dispersions (q2r & matdyn setup)":
            phonon_proc.run_301_phonon_processing()
        elif choice == "Return to Main Menu":
            break

def main():
    while True:
        show_banner()
        category = questionary.select(
            "Select a Category:",
            choices=[
                "1. Pre-Processing",
                "2. Input File Generation",
                "3. Post-Processing",
                "4. Settings",  # Simplest choice string
                "Exit"
            ]
        ).ask()

        if category == "1. Pre-Processing":
            pre_processing_menu()
        elif category == "2. Input File Generation":
            input_generation_menu()
        elif category == "3. Post-Processing":
            post_processing_menu()
        elif category == "4. Settings":  # Now matches the choice exactly
            current_path = config_manager.get_pseudo_dir()
            print(f"\nCurrent Parent Pseudo Path: {current_path}")
            
            new_path = questionary.text("Enter new absolute path for PseudoPotentials:").ask()
            if new_path and os.path.isdir(new_path):
                config_manager.save_pseudo_dir(new_path)
                print("[+] Path updated successfully.")
            else:
                print("[!] Invalid path. No changes made.")
        elif category == "Exit":
            print("\nExiting QE-Kit. Goodbye!")
            sys.exit()
if __name__ == "__main__":
    main()
           

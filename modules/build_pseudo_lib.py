# modules/build_pseudo_lib.py
import os
import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import questionary
from utils import config_manager

# --- CONFIGURATION ---
PS_BASE_URL = "https://pseudopotentials.quantum-espresso.org/legacy_tables/ps-library"
FHI_BASE_URL = "https://pseudopotentials.quantum-espresso.org/legacy_tables/fhi-pp-from-abinit-web-site"
ROOT_DOMAIN = "https://pseudopotentials.quantum-espresso.org"

FUNCTIONAL_MAP = {"pbe": "PBE", "pbesol": "PBE", "lda": "LDA", "pz": "LDA"}

def setup_session():
    """Sets up a requests session with robust retries."""
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=1, status_forcelist=[502, 503, 504])
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session

def get_element_links(session, base_url, filter_str):
    """Fetches element sub-pages from a given base URL."""
    try:
        response = session.get(base_url, timeout=15)
        soup = BeautifulSoup(response.text, 'html.parser')
        elements = []
        for link in soup.find_all('a', href=True):
            href = link['href']
            if filter_str in href:
                parts = href.rstrip('/').split('/')
                if len(parts[-1]) <= 3: # Elements like H, Li, Si
                    elements.append(urljoin(ROOT_DOMAIN, href))
        return sorted(list(set(elements)))
    except Exception as e:
        print(f"[!] Error fetching main list: {e}")
        return []

def download_file(session, url, save_path):
    """Downloads a file if it doesn't already exist."""
    if os.path.exists(save_path):
        return False
    try:
        r = session.get(url, stream=True, timeout=20)
        r.raise_for_status()
        with open(save_path, 'wb') as f:
            for chunk in r.iter_content(8192):
                f.write(chunk)
        return True
    except:
        return False

def build_library():
    print("\n" + "="*60)
    print(" 🌐 QE-POTCAR Builder: Automated Library Construction")
    print("="*60)
    print("[!] WARNING: This process downloads hundreds of files.")
    print("    It may take 15-30 minutes depending on your connection.")
    print("    If running on a remote HPC, it is highly recommended to")
    print("    run this inside a 'screen' or 'tmux' session to prevent")
    print("    interruption if your SSH connection drops.\n")

    confirm = questionary.confirm("Proceed with download?").ask()
    if not confirm:
        return

    # Set default path to ~/QE-POTCAR
    default_path = os.path.expanduser("~/QE-POTCAR")
    target_dir = questionary.text(
        "Enter target installation path:", 
        default=default_path
    ).ask()

    if not target_dir:
        print("[!] No directory provided. Aborting.")
        return

    target_dir = os.path.abspath(target_dir)
    os.makedirs(target_dir, exist_ok=True)
    session = setup_session()

    # ---------------------------------------------------------
    # PHASE 1: PS-Library (PAW / USPP, Scalar / Full Relativistic)
    # ---------------------------------------------------------
    print("\n--- [Phase 1/2] Fetching PS-Library (PBE/LDA, PAW/USPP) ---")
    ps_links = get_element_links(session, PS_BASE_URL, "/legacy_tables/ps-library/")
    total_ps = len(ps_links)
    
    for i, el_link in enumerate(ps_links, 1):
        el_name = el_link.rstrip('/').split('/')[-1].upper()
        print(f"[*] [{i}/{total_ps}] Processing {el_name}...")
        
        try:
            response = session.get(el_link, timeout=15)
            soup = BeautifulSoup(response.text, 'html.parser')
            downloaded_for_this_el = set()
            
            all_links = [urljoin(ROOT_DOMAIN, l['href']) for l in soup.find_all('a', href=True) if l['href'].lower().endswith('.upf')]
            
            for full_url in all_links:
                fname = full_url.split('/')[-1].lower()
                
                # Identify parameters
                func_key = next((k for k in FUNCTIONAL_MAP.keys() if k in fname), None)
                if not func_key: continue
                functional = FUNCTIONAL_MAP[func_key]
                
                ptype = "PAW" if "kjpaw" in fname else "USPP" if "rrkjus" in fname else None
                if not ptype: continue
                
                rel = "full" if "rel-" in fname else "scalar"
                
                # Version preference (psl.1.0.0)
                combo_key = f"{functional}_{rel}_{ptype}"
                is_latest = "psl.1.0.0" in fname
                has_latest = any(func_key in l and "psl.1.0.0" in l and ptype.lower() in l for l in all_links)
                
                if has_latest and not is_latest:
                    continue

                if combo_key not in downloaded_for_this_el:
                    # Direct routing to final directory architecture
                    save_dir = os.path.join(target_dir, functional, rel, ptype)
                    os.makedirs(save_dir, exist_ok=True)
                    save_path = os.path.join(save_dir, fname)
                    
                    if download_file(session, full_url, save_path):
                        print(f"    + {functional} | {rel} | {ptype} -> {fname}")
                        downloaded_for_this_el.add(combo_key)
                        
        except Exception as e:
            print(f"    [!] Error on {el_name}: {e}")

    # ---------------------------------------------------------
    # PHASE 2: FHI Library (NC-mt Norm-Conserving)
    # ---------------------------------------------------------
    print("\n--- [Phase 2/2] Fetching FHI Library (NC-mt Norm-Conserving) ---")
    fhi_links = get_element_links(session, FHI_BASE_URL, "/legacy_tables/fhi-pp-from-abinit-web-site/")
    total_fhi = len(fhi_links)
    
    # pseudo_select.py looks for NC-mt at the top level
    nc_dir = os.path.join(target_dir, "NC-mt")
    os.makedirs(nc_dir, exist_ok=True)

    for i, el_link in enumerate(fhi_links, 1):
        el_name = el_link.rstrip('/').split('/')[-1].upper()
        print(f"[*] [{i}/{total_fhi}] Processing {el_name}...")
        
        try:
            response = session.get(el_link, timeout=15)
            soup = BeautifulSoup(response.text, 'html.parser')
            
            links = soup.find_all('a', href=True)
            for link in links:
                href = link['href']
                if not href.lower().endswith('.upf'): continue
                
                full_url = urljoin(ROOT_DOMAIN, href)
                fname = href.split('/')[-1].lower()

                if 'pbe' in fname:
                    save_path = os.path.join(nc_dir, fname)
                    if download_file(session, full_url, save_path):
                        print(f"    + NC-mt -> {fname}")
                        
        except Exception as e:
            print(f"    [!] Error on {el_name}: {e}")

    # ---------------------------------------------------------
    # PHASE 3: Update System Configuration
    # ---------------------------------------------------------
    print("\n" + "="*60)
    print(f"✅ Download Complete! Library established at: {target_dir}")
    print("="*60)
    
    config_manager.save_pseudo_dir(target_dir)
    print(f"[*] QE-Kit Settings automatically updated to use this path.")

if __name__ == "__main__":
    build_library()

#!/bin/bash

# setup_env.sh - QE-Kit Installer and Configurator
echo "==============================================================="
echo "                  QE-Kit Environment Setup                     "
echo "==============================================================="

# 1. Create and Activate Virtual Environment
if [ ! -d ".venv" ]; then
    echo "[+] Creating virtual environment (.venv)..."
    python3 -m venv .venv
fi

echo "[+] Activating virtual environment..."
source .venv/bin/activate

# 2. Install dependencies and qekit in editable mode
echo "[+] Installing dependencies and qekit..."
pip install --upgrade pip
pip install -e .

# 3. Prompt for the Global Pseudo Directory
echo ""
echo "---------------------------------------------------------------"
echo "CONFIGURATION: Pseudopotential Parent Directory"
echo "Example: /home/user/PseudoPot"
echo "---------------------------------------------------------------"
read -p "Enter the absolute path to your Parent Pseudo folder: " PSEUDO_PATH

# 4. Save the path using the Python utility
if [ -d "$PSEUDO_PATH" ]; then
    python3 -c "from utils import config_manager; config_manager.save_pseudo_dir('$PSEUDO_PATH')"
    echo "[v] Path verified and saved."
else
    echo "[!] Warning: The path provided does not exist."
    echo "    You can update this later within the qekit Settings menu."
    # We still save it so the file is created
    python3 -c "from utils import config_manager; config_manager.save_pseudo_dir('$PSEUDO_PATH')"
fi

echo "==============================================================="
echo " INSTALLATION COMPLETE!"
echo " To start using the tool, run:"
echo "   source .venv/bin/activate"
echo "   qekit"
echo "==============================================================="

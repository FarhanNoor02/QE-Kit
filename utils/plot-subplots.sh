#!/bin/bash

# QE-Kit: Specialized Subplot Master
# Designed for non-homogeneous data (e.g., Optical Constants, Dielectric Functions)

echo "========================================================="
echo "          QE-KIT: STACKED SUBPLOT VISUALIZER            "
echo "========================================================="

# 1. Detect .dat files
FILES=(*.dat)
if [ ! -e "${FILES[0]}" ]; then
    echo "[!] Error: No .dat files found in this directory."
    exit 1
fi

echo "[+] Available data files:"
for i in "${!FILES[@]}"; do echo "  $i) ${FILES[$i]}"; done
read -p "[?] Select file index: " FILE_IDX
SELECTED_FILE=${FILES[$FILE_IDX]}

# 2. Setup Columns & Layout
echo -e "\n[+] File Preview (First 2 lines):"
head -n 2 "$SELECTED_FILE"
echo "---------------------------------------------------------"

read -p "[?] Enter X-column index (usually 1 for Energy): " X_COL
read -p "[?] Enter Y-columns as space-separated list (e.g., 4 6 7): " -a Y_ARRAY
read -p "[?] Enter Overall Figure Title: " MAIN_TITLE

NUM_PLOTS=${#Y_ARRAY[@]}

# 3. Execute Gnuplot with Multiplot Logic
gnuplot -persist <<-EOF
    # Initialize Multiplot
    set multiplot layout $NUM_PLOTS, 1 title "$MAIN_TITLE\nFile: $SELECTED_FILE" font ",12"
    
    # Global Margins and Styles
    set lmargin 12
    set rmargin 8
    set tmargin 2
    set bmargin 1
    set grid
    
    # Professional Spectral Color Palette
    array c[5] = ["#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51"]

    # Loop through each selected column to create a separate pane
    do for [i=1:$NUM_PLOTS] {
        col = word("${Y_ARRAY[*]}", i)
        
        # Clean UI: Only the bottom plot gets the X-axis label
        if (i < $NUM_PLOTS) { 
            unset xlabel 
            set format x ""  # Remove numbers from intermediate X-axes
        } else { 
            set xlabel "Energy (eV) / Column $X_COL" 
            set format x "%g"
            set bmargin 3    # Add room for bottom label
        }
        
        set ylabel "Value (Col ".col.")"
        
        # Plotting command
        plot "$SELECTED_FILE" u $X_COL:col with lines lc rgb c[(i-1)%5+1] lw 2 title "Parameter ".col
    }
    unset multiplot
EOF

echo "[✔] Multi-pane plot generated successfully."

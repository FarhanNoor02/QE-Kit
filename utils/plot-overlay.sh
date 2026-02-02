#!/bin/bash

# QE-Kit: Overlaid Plotter
# Best for: PDOS, Band Structure, and Comparing Multiple Data Series

echo "========================================================="
echo "          QE-KIT: OVERLAID DATA VISUALIZER              "
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

# 2. Setup Columns
echo -e "\n[+] File Preview (First 2 lines):"
head -n 2 "$SELECTED_FILE"
echo "---------------------------------------------------------"

read -p "[?] Enter X-column index (usually 1): " X_COL
read -p "[?] Enter Y-columns as space-separated list (e.g., 2 3 4): " -a Y_ARRAY
read -p "[?] Enter Plot Title: " PLOT_TITLE

NUM_LINES=${#Y_ARRAY[@]}

# 3. Execute Gnuplot with Overlay Logic
gnuplot -persist <<-EOF
    set title "$PLOT_TITLE\nFile: $SELECTED_FILE" font ",12"
    
    # Global Styles
    set grid
    set xlabel "Column $X_COL"
    set ylabel "Intensity / Value"
    set key box opaque top right  # Professional legend box
    
    # Professional color cycle
    array colors[6] = ["#E63946", "#457B9D", "#1D3557", "#06D6A0", "#F4A261", "#A8DADC"]

    # Generate the plot command by looping through the array
    # This handles the simplest case (1 column) or complex multi-column overlays
    plot for [i=1:$NUM_LINES] "$SELECTED_FILE" \
         u $X_COL:(column(word("${Y_ARRAY[*]}", i))) \
         with lines \
         lc rgb colors[(i-1)%6+1] \
         lw 2 \
         title "Column ".word("${Y_ARRAY[*]}", i)
EOF

echo "[✔] Overlaid plot generated successfully."

#!/bin/bash

# Define the working directory and the tutorial directory
PROJECT_DIR="/data"
TUTORIAL_DIR="tutorial3"
ZIP_FILE="toturial3.zip"

# Navigate to the working directory
cd "$PROJECT_DIR" || { echo "Failed to change directory to $PROJECT_DIR. Exiting."; exit 1; }
ech $ZIP_FILE
# Unzip the inputs
if [ -f "$ZIP_FILE" ]; then
    echo "Unzipping input files..."
    unzip -o "$ZIP_FILE" -d "$TUTORIAL_DIR"
else
    echo "$ZIP_FILE does not exist in the current directory. Exiting."
    exit 1
fi

# Navigate to the tutorial directory
cd "$TUTORIAL_DIR" || { echo "Failed to change directory to $TUTORIAL_DIR. Exiting."; exit 1; }

# Equilibrate the solvated complex
echo "Step 1: Energy Minimization"
sander -O -i min.in -o min.out -p ras-raf_solvated.prmtop -c ras-raf_solvated.inpcrd -r min.rst -ref ras-raf_solvated.inpcrd

echo "Step 2: Heating the System"
sander -O -i heat.in -o heat.out -p ras-raf_solvated.prmtop -c min.rst -r heat.rst -x heat.mdcrd -ref min.rst
gzip -9 heat.mdcrd

echo "Step 3: Density Equilibration"
sander -O -i density.in -o density.out -p ras-raf_solvated.prmtop -c heat.rst -r density.rst -x density.mdcrd -ref heat.rst
gzip -9 density.mdcrd

echo "Step 4: Production Equilibration"
sander -O -i equil.in -o equil.out -p ras-raf_solvated.prmtop -c density.rst -r equil.rst -x equil.mdcrd
gzip -9 equil.mdcrd

echo "Processing MD output files"
./process_mdout.pl heat.out density.out equil.out

echo "AMBER MD Simulation steps completed successfully."

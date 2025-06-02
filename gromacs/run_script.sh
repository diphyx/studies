#!/bin/bash

# Set the name of the folder that will be created upon unzipping the archive
case_study_folder="automated-lysozyme-water-forcefields-GROMACS-master"

# Navigate to the /volume/data directory where the ZIP archive is located
cd /volume/data || { echo "Failed to change directory to /volume/data. Exiting."; exit 1; }

# Check if the case study folder exists after extraction
if [ ! -d "$case_study_folder" ]; then
    echo "Expected directory '$case_study_folder' not found after extraction. Exiting."
    exit 1
fi

# Change to the case study directory
cd "$case_study_folder" || { echo "Failed to change directory to $case_study_folder. Exiting."; exit 1; }

# Make sure the current working directory is correctly set
echo "Current working directory: $PWD"

# Ensure all shell scripts in the directory are executable
chmod +x *.sh

# Execute the dynamics script with the specified PDB file and output directory
./script-dynamics.sh files/1aki.pdb files/ || { echo "Dynamics script execution failed. Exiting."; exit 1; }

echo "GROMACS simulation setup completed successfully."
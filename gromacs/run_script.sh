#!/bin/bash

# Set the name of the folder that will be created upon unzipping the archive
case_study_folder="automated-lysozyme-water-forcefields-GROMACS-master"

# Define the path to the ZIP archive
zip_archive="/data/automated-lysozyme-water-forcefields-GROMACS-master.zip"

# Navigate to the /data directory where the ZIP archive is located
cd /data || { echo "Failed to change directory to /data. Exiting."; exit 1; }

# Check if the ZIP archive exists and can be read
if [ ! -r "$zip_archive" ]; then
    echo "ZIP archive $zip_archive not found or is not readable. Exiting."
    exit 1
fi

# Unzip the archive, and check if the unzip operation was successful
unzip -o "$zip_archive" || { echo "Failed to unzip $zip_archive. Exiting."; exit 1; }

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

#!/bin/bash

# Ensure the script stops on error
set -e

# Define the project directory and zip archive name
project_directory="/data"
zip_archive_name="motorBikeSteady.zip"
simulation_folder="motorBikeSteady"

# Navigate to the project directory
cd "$project_directory" || { echo "Failed to change directory to $project_directory. Exiting."; exit 1; }

# Check if the ZIP archive exists
if [ ! -f "$zip_archive_name" ]; then
    echo "ZIP archive $zip_archive_name not found in $project_directory. Exiting."
    exit 1
fi

# Unzip the archive
echo "Unzipping $zip_archive_name..."
unzip "$zip_archive_name" -d "$simulation_folder"

# Check if the simulation folder exists after extraction
if [ ! -d "$simulation_folder" ]; then
    echo "Expected simulation directory '$simulation_folder' not found after extraction. Exiting."
    exit 1
fi

# Change to the simulation directory
cd "$simulation_folder" || { echo "Failed to change directory to $simulation_folder. Exiting."; exit 1; }

# Execute OpenFOAM pre-processing and simulation commands
echo "Starting pre-processing with blockMesh..."
blockMesh

echo "============================ Extracting surface features ============================ "
surfaceFeatures

echo "============================ Generating mesh with snappyHexMesh ============================ "
snappyHexMesh -overwrite

echo "============================ Running the simulation with simpleFoam ============================ "
mpirun --allow-run-as-root simpleFoam

echo "Simulation completed successfully."
touch case.foam
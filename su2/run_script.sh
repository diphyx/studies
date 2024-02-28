#!/bin/bash

# Ensure the script stops on error
set -e

# Number of processors
NP=1 # if greater than 1 then it will use parallel_computation.py

# Define the project directory
project_directory="/data"

# Navigate to the project directory
cd "$project_directory"

# Check the number of processors and run SU2 accordingly
if [[ $NP -gt 1 ]]; then
    echo "Running SU2 in parallel with $NP number of processors"
    # Ensure to replace parallel_computation.py with the correct path if necessary
    parallel_computation.py -f lam_buoyancy_cavity.cfg -n $NP
elif [[ $NP -eq 1 ]]; then
    echo "Running SU2 on a single processor"
    # Ensure to replace SU2_CFD with the correct path if necessary
    SU2_CFD lam_buoyancy_cavity.cfg
else
    echo "Invalid number of processors: $NP"
    exit 1
fi

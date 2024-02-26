#!/bin/bash

# Script to run FastQC on a specified .fastq file and output the results to a specified directory

# Define variables
result_dir="results" # Name of the result directory
input_file="/data/test.fastq" # Path to the input .fastq file

# Welcome message
echo "Starting FastQC analysis..."

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file $input_file does not exist."
    exit 1
fi

# Create the result directory if it doesn't exist
mkdir -p $result_dir
echo "Output will be stored in $result_dir"

# Run FastQC
echo "Running FastQC on $input_file..."
fastqc $input_file -o $result_dir

# Check if FastQC ran successfully
if [ $? -eq 0 ]; then
    echo "FastQC analysis completed successfully."
    echo "Check the $result_dir directory for the results."
else
    echo "FastQC analysis failed."
    exit 1
fi

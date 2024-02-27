#!/bin/bash

# Ensure the script stops on error
set -e

# Check if the input files exist
if [ ! -f "ex1.fa" ]; then
    echo "Reference file ex1.fa not found!"
    exit 1
fi

if [ ! -f "ex1.sam.gz" ]; then
    echo "Input SAM file ex1.sam.gz not found!"
    exit 1
fi

# Step 1: Index the reference FASTA file.
echo "Indexing the reference FASTA..."
samtools faidx ex1.fa

# Step 2: Convert SAM to BAM, sort, and index the BAM file.
# Note: Assuming ex1.sam.gz is gzipped. If not, adjust accordingly.
echo "Converting SAM to BAM, sorting, and indexing..."
samtools sort -o ex1.bam ex1_unsorted.bam
samtools index ex1.bam

# Clean up intermediate files
rm ex1_unsorted.bam

echo "Script completed successfully."
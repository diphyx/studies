#!/bin/bash

# Ensure the script stops on error
set -e
cd /data
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

# Step 2: Convert gzipped SAM to BAM using the reference FASTA for sequence dictionary, sort it, and index.
echo "Converting SAM to BAM, sorting, and indexing..."
# Decompress, convert SAM to BAM with reference for sequence dictionary, then sort and output as ex1.bam
samtools view -bT ex1.fa ex1.sam.gz | samtools sort -o ex1.bam -

# Step 3: Index the sorted BAM file.
samtools index ex1.bam

echo "Script completed successfully."

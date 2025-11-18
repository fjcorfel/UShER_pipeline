#!/bin/bash

# Input multifasta file
INPUT_FILE="$1"  # Replace with your multifasta file name
OUTPUT_DIR="$2"  # Directory to save individual FASTA files

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Split the multifasta file
awk '/^>/ { 
    if (seq) { 
        print seq > filename 
    } 
    seq=""; 
    filename=sprintf("%s/%s.fasta", "'$OUTPUT_DIR'", substr($0, 2)); 
    next 
} 
{ seq = seq $0 "\n" } 
END { if (seq) print seq > filename }' "$INPUT_FILE"

echo "Sequences have been saved in the directory: $OUTPUT_DIR"

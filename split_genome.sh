#!/bin/bash

# Check that an input file has been provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fasta_file> <output_directory> <number_of_parts>"
    exit 1
fi

fasta_file="$1"
output_dir="$2"
num_parts="$3"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Split the FASTA file into parts
seqkit split -p "$num_parts" -O "$output_dir" "$fasta_file"


# The files will be in the format TEST.part_001.fasta, TEST.part_002.fasta, etc.
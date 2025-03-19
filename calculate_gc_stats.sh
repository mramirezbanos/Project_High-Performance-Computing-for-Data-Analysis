#!/bin/bash
# Verify arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <tsv_file> <output_stats>"
    exit 1
fi

# Define input and output files
infile="$1"
statsfile="$2"

# Initialize variables
total_gc=0
total_bases=0
first_line=1

# Read the file line by line
while read -r line; do
    # Skip the header
    if [[ $first_line -eq 1 ]]; then
        first_line=0
        continue
    fi
    # Extract GC content and length from the file
    gc=$(echo $line | awk '{print $3}')  # Third column is GC
    length=$(echo $line | awk '{print $2}')  # Second column is length

    # Calculate GC bases for this window
    gc_bases=$(echo "$gc * $length / 100" | bc)

    # Add to total GC and length variables
    total_gc=$(echo "$total_gc + $gc_bases" | bc)
    total_bases=$(echo "$total_bases + $length" | bc)
done < "$infile"

# Calculate global GC content
global_gc=$(echo "scale=4; $total_gc / $total_bases * 100" | bc)

# Calculate GC median with `awk`
median_gc=$(awk 'BEGIN{n=0} {a[n++]=$3} END{if (n%2==1) print a[int(n/2)]; else print (a[n/2-1]+a[n/2])/2}' "$infile")

# Save results
echo -e "Species\tGlobal_GC\tMedian_GC" > "$statsfile"
echo -e "$(basename "$infile" _gc.tsv)\t$global_gc\t$median_gc" >> "$statsfile"

echo "Global GC calculated: $global_gc"
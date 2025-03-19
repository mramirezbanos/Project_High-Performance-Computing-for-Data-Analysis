import os
import time
from gwf import Workflow

#######################################
#
# GWF workflow to analyze genomes
#
# How to run:
# conda config --append channels gfworg
# conda create -n gwf -y gwf
# conda activate gwf
# gwf -f 04.workflow_gwf.py status
#
#######################################

# Create workflow and indicate resources
gwf = Workflow(defaults={'nodes': 1, 'queue': "normal", 'account': "BiRC_HPC_course"})

# Define output and input directories
input_dir = "/faststorage/project/BiRC_HPC_course/people/marinarb/week06/5genomes"
output_dir = "/faststorage/project/BiRC_HPC_course/people/marinarb/week06/5genomes/seq_results"
stats_dir = "/faststorage/project/BiRC_HPC_course/people/marinarb/week06/5genomes/seq_results"
log_dir = "/faststorage/project/BiRC_HPC_course/people/marinarb/week06/5genomes/seq_results"

# Obtain input files (FASTA)
fasta_files = sorted([f for f in os.listdir(input_dir) if f.endswith(".fasta")])

# Variable to track the previous task (for sequential execution)
prev_task = None

# Log the start time of the entire workflow
start_time = time.time()

# List to store all stats files generated in second task
stats_files = []

# Loop through each FASTA file
for fasta in fasta_files:
    # Extract file name without the file extension
    species_name = os.path.splitext(fasta)[0]

    # Define variables using input (FASTA file) and output (GC content and stats) file paths
    infile = f"{input_dir}/{fasta}"  # Input file path
    outfile = f"{output_dir}/{species_name}_gc.tsv"  # Output file path for GC content
    statsfile = f"{stats_dir}/{species_name}_gc_stats.tsv"  # Output file path for statistics
    combined_stats_file = f"{stats_dir}/combined_gc_stats.tsv"
    log_file = f"{log_dir}/{species_name}_time_log.txt"  # Log file to store time information

    # Task 1: Calculate GC content for each window using seqkit
    task = gwf.target(
        name=f"gc_{species_name}",  # Task name 
        inputs=[infile],  # Input file (FASTA)
        outputs=[outfile],  # Output file (GC content per window)
        walltime="02:00:00",  
        memory="16g", 
        cores=1  
    ) << f"""
    echo "Starting GC calculation for {species_name} at $(date)" >> {log_file}  # Log start time
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate hpc_gwf # Activate conda environment
    mkdir -p {output_dir}  # Create the output directory 
    # Calculate GC content for 25,000 base windows with a 5,000 base step
    time (cat {infile} | seqkit sliding -W 25000 -s 5000 | seqkit fx2tab -n -H -l -C N -g -G > {outfile}) &> {log_file}
    echo "Finished GC calculation for {species_name} at $(date)" >> {log_file}  # Log end time
    """

    # Task 2: Calculate GC median using bash script calculate_gc_stats.sh
    stats_task = gwf.target(
        name=f"gc_stats_{species_name}", 
        inputs=[outfile],  # Input GC content file
        outputs=[statsfile],  # Output file for GC statistics
        walltime="03:00:00",  
        memory="4g", 
        cores=1  
    ) << f"""
    echo "Starting GC statistics calculation for {species_name} at $(date)" >> {log_file}
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate hpc_gwf # Activate conda environment
    # Execute bash script
    bash calculate_gc_stats.sh {outfile} {statsfile}
    echo "Finished GC statistics calculation for {species_name} at $(date)" >> {log_file}
    """

    # Ensure that the statistics task only runs after the GC task
    if prev_task:
        stats_task.inputs = [prev_task.outputs[0]]  # stats_task depends on the previous task's output

    # Update the previous task to be the current stats task
    prev_task = stats_task

    # Collect stats files for final merge
    stats_files.append(statsfile)

# Task 3: Combine all GC stats into one file
combine_task = gwf.target(
    name="combine_gc_stats",  # Give the combine task a unique name
    inputs=stats_files,  # Inputs will be the outputs from all stats tasks
    outputs=[combined_stats_file],
    walltime="00:20:00",
    memory="1g",
    cores=1
) << f"""
echo "Starting the combination of GC stats at $(date)" >> {log_file}
echo -e "Species\tGlobal_GC\tAverage_GC\tMedian_GC" > {combined_stats_file}
# Loop through all stats files and append their contents to the combined stats file
for file in {' '.join(stats_files)}; do
    species=$(basename $file _gc_stats.tsv)  # Extract species name from the filename
    tail -n +2 $file | sed "s/^/$species\t/" >> {combined_stats_file}  # Add species name as a new column
done
echo "Finished combining GC stats at $(date)" >> {log_file}
"""

# Ensure the combination task runs only after all stats tasks
if prev_task:
    combine_task.inputs = [prev_task.outputs[0]]  # Combine task depends on the last stats task

# Log the end time of the entire workflow
end_time = time.time()

# Calculate the total wall time taken by the entire workflow
total_wall_time = end_time - start_time

# Output the total wall time to a file
with open("workflow_time.txt", "w") as f:
    f.write(f"Total wall time for the entire workflow: {total_wall_time} seconds\n")
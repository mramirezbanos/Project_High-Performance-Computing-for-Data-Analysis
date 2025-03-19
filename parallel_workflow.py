import os
from gwf import Workflow

# Create the workflow
gwf = Workflow(defaults={'nodes': 1, 'queue': "normal", 'account': "BiRC_HPC_course"})

# Directories
input_dir = "/faststorage/project/BiRC_HPC_course/people/marinarb/week06/5genomes"
output_dir = "/faststorage/project/BiRC_HPC_course/people/marinarb/week06/5genomes/par_results"
log_dir = "/faststorage/project/BiRC_HPC_course/people/marinarb/week06/5genomes/par_results"
N_PARTS = 4  # Number of parts

# List FASTA files in the input directory
fasta_files = sorted([f for f in os.listdir(input_dir) if f.endswith(".fasta")])

# Create the necessary directories if they don't exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

# Log
log_file = f"{log_dir}/split_log.txt"

# Task 0: Log the start of the workflow
start_time_file = f"{log_dir}/workflow_start_time.txt"

start_task = gwf.target(
    name="workflow_start_time",
    inputs=[],
    outputs=[start_time_file],
    walltime="00:01:00",
    memory="1g",
    cores=1
) << f"""
echo "Starting workflow at $(date)" > {start_time_file}
"""

# Iterate over the FASTA files
for fasta_file in fasta_files:
    species_name = os.path.splitext(fasta_file)[0]  # Use the file name as the species name
    split_dir = f"{output_dir}/{species_name}_split"
    combined_gc_file = f"{output_dir}/{species_name}_gc_combined.tsv"
    stats_file = f"{output_dir}/{species_name}_gc_stats.tsv"

    # Create subdirectory for splitting files
    os.makedirs(split_dir, exist_ok=True)

    # Task 1: Split the genome into N parts
    split_task = gwf.target(
        name=f"split_{species_name}",
        inputs=[f"{input_dir}/{fasta_file}"],
        outputs=[f"{split_dir}/{species_name}.part_{str(i).zfill(3)}.fasta" for i in range(1, N_PARTS + 1)],
        walltime="01:00:00",
        memory="16g",
        cores=7
    ) << f"""
    echo "Starting split for {species_name} at $(date)" >> {log_file}
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate hpc_gwf
    bash split_genome.sh {input_dir}/{fasta_file} {split_dir} {N_PARTS}
    echo "Split finished for {species_name} at $(date)" >> {log_file}
    """

    # Tasks 2: Calculate GC for each part in parallel
    part_outputs = []
    for i in range(1, N_PARTS + 1):
        part_file = f"{split_dir}/{species_name}.part_{str(i).zfill(3)}.fasta"
        part_gc_file = f"{split_dir}/{species_name}.part_{str(i).zfill(3)}_gc.tsv"
        part_outputs.append(part_gc_file)
        combined_gc_stats_file = f"{output_dir}/{species_name}_gc_stats_combined.tsv"

        # Task to calculate GC content for each part
        part_gc_task = gwf.target(
            name=f"gc_{species_name}_part{i}",
            inputs=[f"{split_dir}/{species_name}.part_{str(i).zfill(3)}.fasta"],  # Depends on the split parts
            outputs=[part_gc_file],
            walltime="12:00:00",
            memory="16g",
            cores=7
        ) << f"""
        echo "Starting GC calculation for {species_name} part {i} at $(date)" >> {log_file}
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate hpc_gwf
        time (cat {part_file} | seqkit sliding -W 25000 -s 5000 | seqkit fx2tab -n -H -l -C N -g -G > {part_gc_file}) &>> {log_file}
        echo "GC calculation finished for {species_name} part {i} at $(date)" >> {log_file}
        """

    # Task 3: Combine the GC analysis from each part
    combine_task = gwf.target(
        name=f"combine_gc_{species_name}",
        inputs=part_outputs,  # These are the output files from task 2 (GC parts)
        outputs=[combined_gc_file],  # The combined file
        walltime="05:00:00",
        memory="16g",
        cores=7
    ) << f"""
    echo "Combining GC results for {species_name} at $(date)" >> {log_file}
    cat {' '.join(part_outputs)} > {combined_gc_file}  # Combines the GC files
    echo "Combining GC results finished for {species_name} at $(date)" >> {log_file}
    """

    # Task 4: Calculate GC statistics on the combined file
    gc_stats_task = gwf.target(
        name=f"gc_stats_{species_name}",
        inputs=[combined_gc_file],  # The combined GC file
        outputs=[stats_file],  # The GC statistics file
        walltime="03:00:00",
        memory="16g",
        cores=7
    ) << f"""
    echo "Starting GC statistics calculation for {species_name} at $(date)" >> {log_file}
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate hpc_gwf
    bash calculate_gc_stats.sh {combined_gc_file} {stats_file}
    echo "GC statistics calculation finished for {species_name} at $(date)" >> {log_file}
    """

    # Task 5: Combine GC statistics
    stats_files = [stats_file for _ in range(1, N_PARTS + 1)]  # Statistics files generated in task 4
    combined_stats_file = f"{output_dir}/{species_name}_combined_gc_stats.tsv"  # Final combined stats file

    combine_stats_task = gwf.target(
        name=f"combine_stats_{species_name}",
        inputs=stats_files,  # These are the generated statistics files
        outputs=[combined_stats_file],  # The final combined file
        walltime="05:00:00",
        memory="16g",
        cores=7
    ) << f"""
    echo "Combining GC statistics at $(date)" >> {log_file}
    echo -e "Species\tMedian_GC" > {combined_stats_file}
    for file in {' '.join(stats_files)}; do
        species=$(basename $file _gc_stats.tsv)
        tail -n +2 $file | sed "s/^/$species\t/" >> {combined_stats_file}
    done
    echo "Combining GC statistics finished at $(date)" >> {log_file}
    """

# Task 6: Combine all files into one
all_combined_gc_stats_file = f"{output_dir}/all_combined_gc_stats.tsv"

combine_all_task = gwf.target(
    name="combine_all_gc_stats",
    inputs=[f"{output_dir}/{species_name}_combined_gc_stats.tsv" for species_name in [os.path.splitext(f)[0] for f in fasta_files]],  # All the combined files by species
    outputs=[all_combined_gc_stats_file],  # The final combined file
    walltime="05:00:00",
    memory="16g",
    cores=7
) << f"""
echo "Combining all GC files into one at $(date)" >> {log_file}
echo -e "Species\tMedian_GC" > {all_combined_gc_stats_file}
for file in {' '.join([f"{output_dir}/{os.path.splitext(f)[0]}_combined_gc_stats.tsv" for f in fasta_files])}; do
    species=$(basename $file _combined_gc_stats.tsv)
    tail -n +2 $file | sed "s/^/$species\t/" >> {all_combined_gc_stats_file}
done
echo "Combining all GC files finished at $(date)" >> {log_file}
"""

# Task 7: Log the end time of the workflow
end_time_file = f"{log_dir}/workflow_end_time.txt"

end_task = gwf.target(
    name="workflow_end_time",
    inputs=[all_combined_gc_stats_file],  # The combined stats file should be the last one
    outputs=[end_time_file],
    walltime="00:01:00",
    memory="1g",
    cores=1
) << f"""
echo "Ending workflow at $(date)" > {end_time_file}
"""

# Task 8: Calculate total execution time
total_walltime_file = f"{log_dir}/total_walltime.txt"

total_time_task = gwf.target(
    name="total_walltime",
    inputs=[start_time_file, end_time_file],
    outputs=[total_walltime_file],
    walltime="00:01:00",
    memory="1g",
    cores=1
) << f"""
start_time=$(cat {start_time_file})
end_time=$(cat {end_time_file})
total_time=$(echo "$end_time - $start_time" | bc)
echo "Total execution time: $total_time seconds" > {total_walltime_file}

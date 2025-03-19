# High-Performance-Computing-for-Data-Analysis
Final project of the course "High-Performance Computing for Data Analysis", part of the curriculum of the Master's program in Bioinformatics at Aarhus University.

## Project description
The project aims to analyze the given genomes (in .fasta format) of 5 species (Homo sapiens, Bos taurus, Sus scrofa domesticus, Canis lupus familiaris, and Felis catus) using the tool **Seqkit** and implementing workflows working on **GenomeDK cluster**.

The main goal is to provide answers to the following questions:

Question 1: Which species has the highest median GC?

Question 2: If we assume that only GC% > 2/3 (66.6%) can be potentially damaged by a chemical, which species is most at risk?

To accomplish the goal, we implemented a **serial and a parallel workflow using GWF (Genome Workflow Framework)** and examined how these workflows work. We worked under a Conda environment named hpc_gwf, which contains GWF and Seqkit tools.

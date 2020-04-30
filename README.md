# 5-UTR-Mutation-Analysis

A toolkit for analyzing mutations in 5' UTRs for functional impacts on translation or transcription. 

helper.py contains functions to determine whether a mutation effects various types of elements. 

motif_lists used in this analysis are provided in PWM format. 

run_cluster_slurm_array.pl contains a perl script to submit jobs to a server with srun. 

Motif_Probability_Analysis.py loads in a UTR file, loads a dictionary of trinucleotide contexts for the mutations, generates one permutation, and counts the number of each type of element that gets mutated. These output files can be glued together to generate the distribution for each with get_data.py

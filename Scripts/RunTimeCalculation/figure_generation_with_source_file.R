#!/usr/bin/Rscript

source("./function_of_plot_generation_final.R")

file_path <- "/home/veits/devel/Bioinformatics/ELIXIR_EDAM/WOMBAT-P/WOMBAT-P_Processed/PXD020394/dev/"
file_name <- "execution_trace_2023-03-22_10-10-31.txt"
mapping_file_path <- "./"
mapping_file <- "step_name_task_name_mapping_v2023-03-24.txt"
running_times_path <- file_path
running_times_file_name <- "runningTimes.txt"
output_path <- file_path

create_exe_trace_output(file_path=file_path,
                        file_name=file_name,
                        mapping_file_path = mapping_file_path,
                        mapping_file = mapping_file,
                        running_times_path = running_times_path,
                        running_times_file_name = running_times_file_name,
                        output_path = output_path)






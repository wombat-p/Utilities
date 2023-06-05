#
source("D:/dev/wombatp/wombat_pipelines/final_scripts/function_of_plot_generation_final.R")


  
file_path <- "D:/dev/wombatp/wombat_pipelines/pxd009815/"
file_name <- "execution_trace_2023-03-21_10-11-21.txt"
mapping_file_path <-"D:/dev/wombatp/wombat_pipelines/"
mapping_file <- "step_name_task_name_mapping_v2023-03-24.txt"
running_times_path <- "D:/dev/wombatp/wombat_pipelines/pxd009815/"
running_times_file_name <- "runningTimes.txt"
output_path <- "D:/dev/wombatp/wombat_pipelines/pxd009815/outputs/"

create_exe_trace_output(file_path=file_path,
                        file_name=file_name,
                        mapping_file_path = mapping_file_path,
                        mapping_file = mapping_file,
                        running_times_path = running_times_path,
                        running_times_file_name = running_times_file_name,
                        output_path = output_path)






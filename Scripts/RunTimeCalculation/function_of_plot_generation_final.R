
library(dplyr)
library(neonUtilities)
library(stringr)
library(lubridate)
library(tidyr)
library(viridis)
library(ggplot2)

create_exe_trace_output <- function(file_path,
                                    file_name,
                                    mapping_file_path,
                                    mapping_file,
                                    running_times_path,
                                    running_times_file_name,
                                    output_path) {
  
  # This is the script that parse information about the nextflow run
  # and give a final summary output running processes as bar plots.
  # Input: trace_file.txt   #Output: ggplot() the summary data
  #        merging.txt
  #        maxquant_running_times.txt
  
  # From trace_file only some columns are used for the summarization which are;
  
  ## ## https://www.nextflow.io/docs/latest/tracing.html#execution-report ## ## 
  
  # duration and realtime: ## Time elapsed to complete since the submission.
  # cpu: ## Percentage of CPU used by the process
  # peak_rss: ## Peak of real memory.
  # rchar: ## Number of bytes the process read.
  # wchar: ## Number of bytes the process wrote.
  
  all_final <- NULL
  contain_time <-c("duration","realtime")
  time_type <- c("h","m","s","ms")
  contain_bytes <-c("peak_rss","peak_vmem","rchar","wchar")
  byte_type <- c("TB","GB","MB","KB","B")
  
  # Reading of trace file
  raw_file_elixir <- read.delim(paste0(file_path,file_name),
                                header=TRUE, sep = "\t")
  
  submit_conv <- as.data.frame(strptime(raw_file_elixir$submit,
                                        format = "%Y-%m-%d %H:%M:%S"))
  colnames(submit_conv) <-"submit"
  
  # Time Conversion: everything is converted to minute.
  
  for (d in 1:length(contain_time)){
    temp1<-NULL
    
    for (t in 1:length(time_type)){
      assign(paste0("temp_", time_type[t]),str_extract_all(raw_file_elixir[,d+7], 
                                                           time_type[t]))
    }
    for (e in 1:dim(raw_file_elixir)[1]){
      # for hours, minutes and seconds
      if (is.element("h",unlist(temp_h[[e]]))){
        
        query_hms <- str_match(raw_file_elixir[,d+7][e],'(\\d+)h\\s(\\d+)m\\s(\\d+)s')
        temp1[e] <- as.numeric(query_hms[2]) *60 + as.numeric(query_hms[3])
        + as.numeric(query_hms[4])/60}
      # for microseconds only
      else if (is.element("ms",unlist(temp_ms[[e]])) && is.element("s", unlist(temp_s[[e]])) && 
               is.element("m", unlist(temp_m))){
        
        query_msec <- str_extract(raw_file_elixir[,d+7][e],'[[:digit:]]+\\.*[[:digit:]]*')
        temp1[e] <- as.numeric(query_msec) *0.00001666666
      }
      # for minutes only
      else if ((is.element("m",unlist(temp_m[[e]]))) && !is.element("ms",unlist(temp_ms[[e]])) && 
               !is.element("s", unlist(temp_s[[e]]))){
        query_m <- str_match(raw_file_elixir[,d+7][e],'[[:digit:]]+\\.*[[:digit:]]*')
        temp1[e] <- as.numeric(query_m)
        
      }
      # for seconds only
      else if ((is.element("s",unlist(temp_s[[e]]))) && 
               (!is.element("m",unlist(temp_m[[e]]))) && !is.element("ms",unlist(temp_ms[[e]]))){
        query_s <- str_extract(raw_file_elixir[,d+7][e],'[[:digit:]]+\\.*[[:digit:]]*')
        temp1[e] <- as.numeric(query_s)/60
      }
      # for both minutes and seconds
      else if ((is.element("m",unlist(temp_m[[e]]))) && !is.element("ms",unlist(temp_ms[[e]]))){
        query_ms <- str_match(raw_file_elixir[,d+7][e],'(\\d+)m\\s(\\d+)s')
        temp1[e] <- as.numeric(query_ms[2]) + as.numeric(query_ms[3])/60}
      
      assign(paste0(contain_time[d],"_conv"),as.data.frame(temp1))
    }
    rm(temp1)
  }
  
  # Byte Conversion: everything is converted to Gigabyte.
  
  for (j in 1:length(contain_bytes)){
    
    # To prevent NA in the index list of each contain_bytes column
    raw_file_elixir[,j+10][raw_file_elixir[,j+10] == 0] <- "0 B"
    
    df <- assign(paste0("sep_",contain_bytes[j]),
                 as.data.frame(t(as.data.frame(sapply(strsplit(raw_file_elixir[,j+10], " "),
                                                      "[", 1:2)))))
    colnames(df) <- c("num","MB_name")
    temp2 <- NULL
    ## Source for the conversion: https://www.checkyourmath.com/convert/digital_storage/kilobytes_bytes.php
    for (d in 1:length(byte_type)){
      
      indx_list <- str_extract_all(df$MB_name, paste0("^",byte_type[d]))
      for (c in 1:(dim(raw_file_elixir)[1])){
        
        if (is.element("GB",unlist(indx_list[[c]]))){
          
          temp2[c] <- as.numeric(df$num[c])
          
        }else if (is.element("MB",unlist(indx_list[[c]]))){
          
          temp2[c] <- as.numeric(df$num[c])*0.0009765625
          
        }else if(is.element("KB",unlist(indx_list[[c]]))){
          
          temp2[c] <- as.numeric(df$num[c])*0.0000009537
          
        }else if(is.element("B",unlist(indx_list[[c]]))){ 
          
          temp2[c] <- as.numeric(df$num[c])*0.0000000009
          
        }else if(is.element("TB",unlist(indx_list[[c]]))){
          temp2[c] <- as.numeric(df$num[c])*1024
        }
      }
      
    }
    assign(paste0(contain_bytes[j],"_gb"),temp2)
    rm(temp2)
  }
  
  # For percentage of cpu, only "%" symbol is removed.
  xcpu <- gsub("%","",raw_file_elixir$X.cpu)
  xcpu_tmp <- as.numeric(unlist(xcpu))
  
  rmv_num_name <- as.data.frame(sapply(strsplit(raw_file_elixir[,"name"], " "), "[", 1))
  
  # Converted version of actual raw file
  refined_raw_file_elixir <- cbind(rmv_num_name[,1],duration_conv$temp1,
                                   xcpu_tmp,peak_rss_gb,rchar_gb, 
                                   wchar_gb,as.data.frame.POSIXct(submit_conv$submit))
  
  
  colnames(refined_raw_file_elixir) <-  c(colnames(raw_file_elixir[c(4,8)]), 
                                          "cpu",colnames(raw_file_elixir[,c(11,13,14)]),
                                          "submit")
  col_name <- colnames(refined_raw_file_elixir)
  
  ## The part where calculation is done.
  # Summation = wchar and rchar
  # Mean = peak_rss and cpu
  # Duration = max(submit) - min(submit) + max(duration)
  
  final1 <- NULL
  final2 <- NULL
  final3 <- NULL
  
  
  for (i in 3:length(col_name)){
    
    if (col_name[i] == "submit"){
      
      tmpd <- assign(paste0("total_duration",col_name[i]), 
                     as.data.frame(refined_raw_file_elixir) %>% 
                       group_by(name) %>% summarise(as.numeric(difftime(max(submit),
                                                                        min(submit),units = "mins"))+
                                                      max(duration)))
      
      tmpd <- cbind(col_name[i],tmpd)
      final1 <- bind_rows(final1,tmpd)
      
    }else if (col_name[i] == "peak_rss" || col_name[i] == "cpu"){
      tmpps <- assign(paste0("stas_calc_",col_name[i]), 
                      as.data.frame(refined_raw_file_elixir) %>%
                        group_by(name) %>% summarise(mean=mean(get(col_name[i]))))
      tmpps <- cbind(col_name[i],tmpps)
      final2 <- bind_rows(final2,tmpps)
      
    }else{
      tmprest <- assign(paste0("stas_calc_",col_name[i]), 
                        as.data.frame(refined_raw_file_elixir) %>%
                          group_by(name) %>% summarise(sum=sum(get(col_name[i]))))
      tmprest <- cbind(col_name[i],tmprest)
      final3 <- bind_rows(final3,tmprest) 
    }
    
  }
  rm(list = ls()[!(ls() %in% c('all_final','list_wombat','final1','final2','final3',
                               'contain_time','contain_bytes','time_type',
                               'col_name','byte_type','raw_file_elixir',
                               'refined_raw_file_elixir', 'mapping_file_path',
                               'mapping_file','running_times_path',
                               'running_times_file_name','output_path'))])
  all_final <- bind_rows(final1,final2,final3)
  
  # Formatting of Wombat steps
  all_final_sft_name <- all_final %>% separate(name,sep=":", c('Nextflow',"Wombat",
                                                               'Software_name',"Step_name"))
  
  all_final_sft_name$Step_name[is.na(all_final_sft_name$Step_name)] <- "OTHER"
  all_final_sft_name <- all_final_sft_name[,-c(2,3)]
  
  # Merge step_name and software name for further merging below:
  
  all_final_sft_name[,"merged_name"] <- paste(all_final_sft_name$Software_name,
                                              all_final_sft_name$Step_name, sep =  "_")
  
  #colnames(all_final_sft_name) <- c("trace_report","task_name","real_duration","software_name","mean","sum")
  
  colnames(all_final_sft_name)[c(1,4)] <- c("task_name","duration")
  #c("Software_name","Step_name","real_duration_min","mean_gb","sum_gb")
  
  all_final_sft_name$mean[is.na(all_final_sft_name$mean)] <- ""
  all_final_sft_name$sum[is.na(all_final_sft_name$sum)] <- ""
  
  all_final_sft_name$duration <- as.numeric(all_final_sft_name$duration)
  all_final_sft_name$duration[is.na(all_final_sft_name$duration)] <- ""
  
  ### BE CAREFUL of FILE VERSION ###
  
  df_for_merging <- read.delim(paste0(mapping_file_path,
                                      mapping_file),
                               header=TRUE,
                               sep = "\t")
  
  df_for_merging[,"merged_name"] <- paste(df_for_merging$Software_name,
                                          df_for_merging$Step_name,sep = "_")
  
  # Software_name column in the df_for_merging object creates a redundancy 
  # When the task_names are shared with more than one software.
  # We always have to believe "Software_name.x" which comes from 
  # all_final_sft_name object.) But always check these columns to make sure 
  # that everything is merged correctly.
  
  merge_with_conv_min_gb <- merge(all_final_sft_name,df_for_merging,
                                  by="merged_name")
  merge_with_conv_min_gb <- merge_with_conv_min_gb[,c(1:3,5:7,10)]
  
  merge_with_conv_min_gb[,"combine_values"] <- paste0(merge_with_conv_min_gb$duration,
                                                      merge_with_conv_min_gb$mean,
                                                      merge_with_conv_min_gb$sum)
  
  
  ### TO CREATE SEPARATE OBJECT of EACH TASK NAME
  
  submit_summation <- merge_with_conv_min_gb %>% 
    filter(task_name %in% "submit") %>%
    filter(!grepl("CALCBENC|PREP|SDRF",Software_name.x)) %>%
    group_by(optimal_task_name_final,Software_name.x) %>% 
    summarise(sum(as.numeric(combine_values)))
  
  
  merge_with_conv <- merge_with_conv_min_gb %>% 
    filter(!grepl("submit", task_name)) %>%
    filter(!grepl("CALCBENC|PREP|SDRF",Software_name.x)) %>%
    select(task_name,Software_name.x,combine_values) %>%
    group_by(Software_name.x, task_name) %>%
    summarise(values=sum(as.numeric(combine_values)))
  
  ###### DURATION CALCULATION BY INVOLVING MAXQUANT ADDITIONAL FILE ######
  # Remove summation info of Maxquant
  submit_summation_rmv_mq <- submit_summation %>% filter(!grepl("MAXQUANT",Software_name.x))
  
  ###### DURATION CALCULATION BY INVOLVING MAXQUANT ADDITIONAL FILE ######
  
  # Read additional file for merging
  mq_run_time <- read.delim(paste0(running_times_path,running_times_file_name),
                            header=TRUE,
                            sep = "\t")
  
  ### BE CAREFUL WITH WHITE SPACE!! THEY ARE IMPORTANT FOR MATCHING ###
  ### ORIJINAL RUNNING_FILE.TXT FROM MAXQUANT HAS DOUBLE WHITE SPACES, FILTER() FUNCTION IS SENSITIVE THESE SPACES.
  peaklist_creation <- c("MS/MS preparation ",
                         "Combining apl files for first search ",
                         "Read search results for recalibration ",
                         "Mass recalibration ",
                         "Calculating masses ",
                         "Preparing searches "
  )
  MSMS_search <- c("Testing fasta files ",
                   "Combining apl files for main search ",
                   "Combining second peptide files ",
                   "Correcting errors ",
                   "Finish search engine results ",
                   "Finish search engine results (SP) ",
                   "MS/MS first search ",
                   "MS/MS main search ",
                   "MS/MS preparation for main search ",
                   "Preparing combined folder  ",
                   "Preparing reverse hits ",
                   "Reading search engine results ",
                   "Reading search engine results (SP) ",
                   "Second peptide search ",
                   "Assembling second peptide MS/MS "
  )
  
  other <- c("Configuring ","Assemble run info ","Finish run info ")
  
  quantification <- c("Testing raw files ",
                      "Feature detection ",
                      "Deisotoping ",
                      "Calculating peak properties ",
                      "Re-quantification ",
                      "Reporter quantification ",
                      "Retention time alignment ",
                      "Matching between runs 1 ",
                      "Matching between runs 2 ",
                      "Matching between runs 3 ",
                      "Matching between runs 4 ",
                      "Prepare protein assembly ",
                      "Assembling proteins ",
                      "Assembling unidentified peptides ",
                      "Finish protein assembly ",
                      "Updating identifications ",
                      "iBAQ ",
                      "Label-free preparation ",
                      "Label-free normalization ",
                      "Label-free quantification ",
                      "Label-free collect ",
                      "Estimating complexity ")
  
  MSMS_search_validation <- c("Filter identifications (MS/MS) ",
                              "Calculating PEP ",
                              "Copying identifications ",
                              "Applying FDR ",
                              "Filtering identifications (SP) ",
                              "Applying FDR (SP) ")
  
  writing_tables <- c("Testing raw files ",
                      "Writing tables ",
                      "Finish writing tables ")
  
  mapping_category <- c("peaklist_creation", 
                        "MSMS_search",
                        "other",
                        "quantification",
                        "MSMS_search_validation",
                        "writing_tables")
  
  mq_run_time_mapped <- NULL
  
  # Creation for mapping data frame using maxquant input
  for (i in 1:length(mapping_category)){
    temp <- filter(mq_run_time, Job %in% get(mapping_category[i])) %>% 
      mutate(task_category=paste0(mapping_category[i]))
    
    mq_run_time_mapped <- rbind(mq_run_time_mapped,temp)
    
  }
  
  mq_run_time_mapped <- cbind(mq_run_time_mapped,"MAXQUANT")
  
  # Summation of running times based on task_category
  mq_run_time_classed <- mq_run_time_mapped %>% 
    group_by(task_category,`"MAXQUANT"`) %>% 
    summarise(sum(Running.time..min.))
  
  colnames(mq_run_time_classed) <- colnames(submit_summation_rmv_mq) 
  
  # Final output of MaxQuant that contains submit information 
  submit_summation_add_mq <- rbind(submit_summation_rmv_mq,mq_run_time_classed)
  
  colnames(submit_summation_add_mq)[3] <- "values"
  
  order_x_axis_label_duration <- c("workflow_configuration",
                                   "wombatp_prep",
                                   "peaklist_creation",
                                   "MSMS_search",
                                   "MSMS_search_validation",
                                   "quantification",
                                   "statistical_testing",
                                   "writing_tables",
                                   "other",
                                   "benchmark")
  
  ### FUNCTION FOR BAR PLOT GENERATION ###
  final_barplot <- function(data,
                            x,
                            y,
                            order_x_axis_label,
                            fill,
                            title_name,
                            inc,
                            ylabel,
                            ymax) {
    
    ggplot(data ={{data}}, aes(y={{y}}, x=factor({{x}}, levels = order_x_axis_label),fill={{fill}})) + 
geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding=0.01))  +  
      theme_light() + coord_flip() +
      theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            plot.title = element_text(size=15),
            legend.title=element_text(size=15),
            axis.text.x =element_text(size=15,angle = 90),
            axis.text.y = element_text(size = 15, angle = 0, vjust=0.5),
            axis.title=element_text(size=15)
      )  + # scale_y_continuous(breaks = seq(from= 0,to= ymax,by=inc)) +
      ggtitle(title_name) +
      labs(fill="Wombat pipeline steps",x= "", y= ylabel) +
      scale_fill_brewer(type="seq", palette = "Set1") +
      scale_x_discrete(limits = levels(order_x_axis_label),
                       labels = function(order_x_axis_label) str_wrap(order_x_axis_label, width = 10))
    
  }
  
  #### GGPLOT FIGURES ####
  
  # THIS IS FOR DURATION: THIS IS THE ONLY ONE THAT USE DIFFERENT X AXIS LABELS
  assign(paste0("submit_p"),final_barplot(data = submit_summation_add_mq %>% 
                                            filter(!grepl("CALCBENC|PREP|SDRF",Software_name.x)) %>%
                                            filter(!grepl("workflow_configuration|other", optimal_task_name_final)),
                                          x=optimal_task_name_final,
                                          y=values,
                                          order_x_axis_label= order_x_axis_label_duration,
                                          fill=Software_name.x,
                                      title_name = "Running time of different steps",
                                          inc = 20,
                                          ylabel = "Running Time (min)",
                                          ymax = max(submit_summation_add_mq$values)
  ))
  
  
  
  assign(paste0("rest_p"),final_barplot(data = merge_with_conv,
                                        x=merge_with_conv$task_name,
                                        y=merge_with_conv$values,
                                        order_x_axis_label= unique(merge_with_conv$task_name),
                                        fill=merge_with_conv$Software_name.x,
                                    title_name = "Summation of all task names except running time across all MS/MS softwares",
                                        inc = 100,
                                        ylabel = " ",
                                        ymax = max(merge_with_conv$values)
  ))
  
  
  
  ### REST IS CREATED INSIDE FOR LOOP ### 
  
  # ### WITHOUT DURATION ###
  ylabels <- c("GB",
               "%",
               "GB",
               "GB"
  )
  
  ### WITHOUT DURATION ###
  title_names <- c( "real memory peaks",
                    "Used CPUs (%)",
                    "GBs read",
                    "GBs write"
  )
final_title_names <- as.data.frame(sapply(title_names,function(x) paste("Total of",x)))
  
  report_calc <- c("submit","peak_rss","cpu","rchar","wchar")
  report_calc_wo_submit <- report_calc[-1]
  
  for (i in 1:(length(report_calc_wo_submit))){
    
    if (is.element("peak_rss",report_calc_wo_submit[i]) || is.element("wchar",report_calc_wo_submit[i])){
      
      assign(paste0(report_calc_wo_submit[i],"_p"),
             final_barplot(merge_with_conv %>%
                             filter(grepl(report_calc_wo_submit[i],task_name)),
                           x= task_name,
                           y=values,
                           order_x_axis_label = report_calc_wo_submit[i],
                           fill = Software_name.x,
                           title_name = final_title_names[i,],
                           ylabel=ylabels[i],
                           inc = 10,
                           ymax = max(merge_with_conv$values)))
      #print(report_calc_wo_submit[i])
      
    }else{
      assign(paste0(report_calc_wo_submit[i],"_p"),
             final_barplot(merge_with_conv %>%
                             filter(grepl(report_calc_wo_submit[i],task_name)),
                           x= task_name,
                           y=values,
                           order_x_axis_label = report_calc_wo_submit[i],
                           fill = Software_name.x,
                           title_name = final_title_names[i,],
                           ylabel=ylabels[i],
                           inc = 100,
                           ymax = max(merge_with_conv$values)))
      
    }
  }
  
  sapply(report_calc,function(x) ggsave(filename = paste0(x,"_p",".pdf"),
                                        #                                        width = 50, height = 40, 
                                        path = output_path,
                                        units = "cm",
                                        get(paste0(x,"_p")),
                                        device = "pdf", #".svg"
  ))
  
  # MERGE ALL FIGURES INTO ONE
  library(ggpubr)
  gg_all <- ggarrange(ggarrange(get(paste0(report_calc[1],"_p")), heights=2,
                                common.legend=T),
                      ggarrange(get(paste0(report_calc[2],"_p")),
                                get(paste0(report_calc[3],"_p")),
                                get(paste0(report_calc[4],"_p")),
                                get(paste0(report_calc[5],"_p")),
                              common.legend=T),
                      common.legend = TRUE)
  
  ggsave(filename = paste("combined",".pdf",sep = "_"),
     width = 30, height = 20, 
         path =  output_path,
         units = "cm",
         gg_all,
         device = "pdf", #".svg"
  )
  
}


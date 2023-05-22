library("dplyr")
library("neonUtilities")
library("stringr")
#library("kableExtra")
library("reshape2")
library("lubridate")
#library("patchwork")
library("tidyr")

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
raw_file_elixir <- read.delim("execution_trace_file.txt",
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
                             'refined_raw_file_elixir'))])
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

df_for_merging <- read.delim("step_name_task_name_mapping.txt",
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


report_calc <- c("submit","peak_rss","cpu","rchar","wchar")

#test <- merge_with_conv_min_gb %>% filter(task_name %in% "rchar") %>% group_by(optimal_task_name_final) %>% summarise(sum(as.numeric(combine_values)))

for (k in 1:length(report_calc)){
  assign(paste0(report_calc[k],"_summation"), merge_with_conv_min_gb %>% 
           filter(task_name %in% report_calc[k]) %>%
    group_by(optimal_task_name_final,Software_name.x) %>% 
    summarise(sum(as.numeric(combine_values))))
  
}
 
submit_summation_rmv_mq <- submit_summation %>% filter(!grepl("MAXQUANT",Software_name.x))

###### DURATION CALCULATION BY INVOLVING MAXQUANT ADDITIONAL FILE ######

mq_run_time <- read.delim("D:/dev/wombatp/wombat_pipelines/veit_trace_report_maxquant_res/data2_veit_PXD009815_work_6a_92f/#runningTimes.txt",
                          header=TRUE,
                          sep = "\t")
### BE CAREFUL WITH WHITE SPACE!! THEY ARE IMPORTANT FOR MATCHING ###
### ORIJINAL RUNNING_FILE.TXT FROM MAXQUANT HAS DOUBLE WHITE SPACES, FILTER() FUNCTION IS SENSITIVE THESE SPACES.
MSMS_search <- c("Testing fasta files ",
                 "Calculating masses ",
                 "Combining apl files for first search ",
                 "Combining apl files for main search ",
                 "Combining second peptide files ",
                 "Correcting errors ",
                 "Finish search engine results ",
                 "Finish search engine results (SP) ",
                 "Mass recalibration ",
                 "MS/MS first search ",
                 "MS/MS main search ",
                 "MS/MS preparation for main search ",
                 "Preparing combined folder  ",
                 "Preparing reverse hits ",
                 "Read search results for recalibration ",
                 "Reading search engine results ",
                 "Reading search engine results (SP) ",
                 "Second peptide search ",
                 "Assembling second peptide MS/MS ",
                 "MS/MS preparation ")

other <- c("Configuring ","Assemble run info ","Finish run info ")

quantification <- c("Feature detection", 
                    "Deisotoping ",
                    "Calculating peak properties ",
                    "Preparing searches ",
                    "Reporter quantification ",
                    "Prepare protein assembly ",
                    "Assembling proteins ",
                    "Assembling unidentified peptides ",
                    "Finish protein assembly ",
                    "Updating identifications ",
                    "Label-free preparation ",
                    "Label-free normalization ",
                    "Label-free quantification ",
                    "Label-free collect ",
                    "Estimating complexity ",
                    "Prepare writing tables  ",  
                    "Re-quantification ")

statistical_testing <- c("Calculating PEP " ,
               "Filter identifications (MS/MS) ",
               "Filtering identifications (SP) ",
               "Copying identifications ", 
               "Applying FDR ",
               "Applying FDR (SP) ")

writing_tables <- c("Testing raw files ",
             "Writing tables ",
             "Finish writing tables  ")

mapping_category <- c("MSMS_search","other","quantification","statistical_testing","writing_tables")

mq_run_time_mapped <- NULL

# Creation for mapping data frame using maxquant input
for (i in 1:length(mapping_category)){
  temp <- filter(mq_run_time, Job %in% get(mapping_category[i])) %>% mutate(task_category=paste0(mapping_category[i]))
  mq_run_time_mapped <- rbind(mq_run_time_mapped,temp)

}

mq_run_time_mapped <- cbind(mq_run_time_mapped,"MAXQUANT")

# Summation of running times based on task_category
mq_run_time_classed <- mq_run_time_mapped %>% 
  group_by(task_category,`"MAXQUANT"`) %>% summarise(sum(Running.time..min.))

colnames(mq_run_time_classed) <- colnames(submit_summation_rmv_mq) 

submit_summation_add_mq <- rbind(submit_summation_rmv_mq,mq_run_time_classed)


library(ggplot2)
library(hrbrthemes)  

### FUNCTION FOR BOX PLOT GENERATION ###
final_barplot <- function(data, x, y, fill, title_name, ylabel,ymax,order_x_axis_label,inc) {
  ggplot({{data}}, aes(y={{y}}, x=factor({{x}}, levels = order_x_axis_label),fill={{fill}})) + 
    geom_bar(width= .5,stat = "identity", position = position_dodge2(preserve = "single",width = 0.7))  +  
    theme_light() + coord_flip() +
    theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size=20),
          legend.title=element_text(size=15),
          axis.text.x =element_text(size=15,angle = 90),
          axis.text.y = element_text(size = 15),
          axis.title=element_text(size=15)
    ) + scale_y_continuous(breaks = seq(0,ymax,inc),limits = c(0,ymax)) +
    ggtitle(title_name) +
    labs(fill="Wombat pipeline steps",x= "Task Name", y= ylabel) +
    scale_fill_manual(breaks=c("COMPOMICS","PROLINE","TPP","MAXQUANT",
                               "CALCBENCHMARKS_COMPOMICS","CALCBENCHMARKS_PROLINE",
                               "CALCBENCHMARKS_TPP","CALCBENCHMARKS_MQ",
                               "PREPARE_FILES","SDRFMERGE"),
                      values = c("#E69F00","#56B4E9","#009E73","#CC79A7",
                                 "#999999","#999999","#999999","#999999",
                                 "#000000","#000000","#000000","#000000"
                      )) + scale_x_discrete(limits = levels(order_x_axis_label),
                                            labels = function(order_x_axis_label) str_wrap(order_x_axis_label, width = 10))
  
}

order_x_axis_label_duration <- c("wombatp_prep",
                                 "peaklist_creation",
                                 "MSMS_search",
                                 "quantification",
                                 "statistical_testing",
                                 "writing_tables",
                                 "other",
                                 "benchmark")
# For rest of them, MQ has one additional label which is combination of all steps.
order_x_axis_label_rest <- c("wombatp_prep",
                             "peaklist_creation",
                             "MSMS_search + quantification + statistical_testing",
                             "MSMS_search",
                             "quantification",
                             "statistical_testing",
                             "writing_tables",
                             "other",
                             "benchmark")

### WITHOUT DURATION ###
ylabels <- c("peak rss (GB)",
             "percentage of cpu",
             "rchar usage (GB)",
             "wchar usage (GB)"
             )

### WITHOUT DURATION ###
title_names <- c( "number of real memory peaks ",
                  " percentage of used CPUs ",
                  "number of GBs the process read",
                  "number of GBs the process write"
                  )
final_title_names <- as.data.frame(sapply(title_names,function(x) paste("Total",x,"across all MS/MS softwares")))

                                          #### GGPLOT FIGURES ####

# THIS IS FOR DURATION: THIS IS THE ONLY ONE THAT USE DIFFERENT X AXIS LABELS
assign(paste0(report_calc[1],"_p"),
       final_barplot(submit_summation_add_mq %>%
                       filter(!grepl("CALCBENC|PREP|SDRF",Software_name.x)),
                     optimal_task_name_final,
                     `sum(as.numeric(combine_values))`,
                     Software_name.x,
                     "Total of running time across all MS/MS softwares",
                     "Running Time (min)",
                     max(submit_summation_add_mq$`sum(as.numeric(combine_values))`),
                     order_x_axis_label_duration,
                     inc=50)
              )
### REST IS CREATED INSIDE FOR LOOP ### 

report_calc_wo_submit <- report_calc[-1]
#rm(report_calc)

for (i in 1:(length(report_calc_wo_submit))){
  if (is.element("peak_rss",report_calc_wo_submit[i]) || is.element("wchar",report_calc_wo_submit[i])){
    inc_a=10
    assign(paste0(report_calc_wo_submit[i],"_p"),
           final_barplot(get(paste0(report_calc_wo_submit[i],"_summation")) %>%
                           filter(!grepl("CALCBENC|PREP|SDRF",Software_name.x)),
                         optimal_task_name_final,
                         `sum(as.numeric(combine_values))`,
                         Software_name.x,
                         final_title_names[i,],
                         ylabels[i],
                         max(get(paste0(report_calc_wo_submit[i],"_summation"))$`sum(as.numeric(combine_values))`),
                         order_x_axis_label_rest,
                         inc = inc_a))
    #print(report_calc_wo_submit[i])
    
  }else{
    inc_b=100
    assign(paste0(report_calc_wo_submit[i],"_p"),
           final_barplot(get(paste0(report_calc_wo_submit[i],"_summation")) %>%
                           filter(!grepl("CALCBENC|PREP|SDRF",Software_name.x)),
                         optimal_task_name_final,
                         `sum(as.numeric(combine_values))`,
                         Software_name.x,
                         final_title_names[i,],
                         ylabels[i],
                         max(get(paste0(report_calc_wo_submit[i],"_summation"))$`sum(as.numeric(combine_values))`),
                         order_x_axis_label_rest,
                         inc = inc_b))
    
  }
}

sapply(report_calc,function(x) ggsave(filename = paste0(x,"_p",".tiff"),
                                                width = 50, height = 40, 
                                                path = "/The/Path/where/figures/will/be/saved",
                                                units = "cm",
                                                get(paste0(x,"_p")),
                                                device = "tiff", #".svg"
                                                ))

# MERGE ALL FIGURES INTO ONE
gg_all <- ggpubr::ggarrange(get(paste0(report_calc[1],"_p")),
                  get(paste0(report_calc[2],"_p")),
                  get(paste0(report_calc[3],"_p")),
                  get(paste0(report_calc[4],"_p")),
                  get(paste0(report_calc[5],"_p")),
                  ncol = 1,
                  common.legend = TRUE)

ggsave(filename = paste("combined",".tiff",sep = "_"),
       width = 100, height = 80, 
       path =  "/The/Path/where/figures/will/be/saved",
       units = "cm",
       gg_all,
       device = "tiff", #".svg"
)

#### FIGURES ARE CREATED SEPARATELY BY ADDING Y-LABEL and TITLE INDIVIDUVALLY #### 
# 
# 
# assign(paste0(report_calc[2],"_p"),
#        final_barplot(get(paste0(report_calc[2],"_summation")) %>% 
#                        filter(!grepl("CALCBENC|PREP|SDRF",Software_name.x)),
#                                                    optimal_task_name_final,
#                                                    `sum(as.numeric(combine_values))`,
#                                                    Software_name.x,
#                                                    "Total number of real memory peaks across all MS/MS softwares",
#                                                    "peak rss (GB)",
#                                                    max(get(paste0(report_calc[2],"_summation"))$`sum(as.numeric(combine_values))`),
#                                                    order_x_axis_label_rest,
#                                                    inc = 10))
# 
# ##!!!!## Rather than summation of CPU, it was taken mean of % CPU. ##!!!!## 
# assign(paste0(report_calc[3],"_p"),
#        final_barplot(get(paste0(report_calc[3],"_summation")) %>% 
#                        filter(!grepl("CALCBENC|PREP|SDRF",Software_name.x)),
#                      optimal_task_name_final,
#                      `sum(as.numeric(combine_values))`,
#                      Software_name.x,
#                      "Total percentage of used CPUs across all MS/MS softwares",
#                      "percentage of cpu",
#                      max(get(paste0(report_calc[3],"_summation"))$`sum(as.numeric(combine_values))`),
#                      order_x_axis_label_rest,
#                      inc=100))
# 
# assign(paste0(report_calc[4],"_p"),
#        final_barplot(get(paste0(report_calc[4],"_summation")) %>% 
#                        filter(!grepl("CALCBENC|PREP|SDRF",Software_name.x)),
#                      optimal_task_name_final,
#                      `sum(as.numeric(combine_values))`,
#                      Software_name.x,
#                      "Total number of GBs the process read across all MS/MS softwares",
#                      "rchar usage (GB)",
#                      max(get(paste0(report_calc[4],"_summation"))$`sum(as.numeric(combine_values))`),
#                      order_x_axis_label_rest,
#                      inc=100))
# 
# assign(paste0(report_calc[5],"_p"),
#        final_barplot(get(paste0(report_calc[5],"_summation")) %>% 
#                        filter(!grepl("CALCBENC|PREP|SDRF",Software_name.x)),
#                      optimal_task_name_final,
#                      `sum(as.numeric(combine_values))`,
#                      Software_name.x,
#                      "Total number of GBs the process write across all MS/MS softwares",
#                      "wchar usage (GB)",
#                      max(get(paste0(report_calc[5],"_summation"))$`sum(as.numeric(combine_values))`),
#                      order_x_axis_label_rest,
#                      inc=10))


#### EXAMPLE OF BOXPLOT FUNCTION

# #p <- submit_summation_add_mq %>% filter(!grepl("CALCBENC",Software_name.x)) %>% 
# ggplot(., aes(y=`sum(as.numeric(combine_values))`,
#               x=optimal_task_name_final,
#               fill=Software_name.x)) + 
#   geom_bar(width= .5, stat = "identity", position = position_dodge2(preserve = "single",width = 0.7))  +  
#   theme_light() + coord_flip() +
#   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
#         axis.title.x = element_text(size = 15),
#         axis.title.y = element_text(size = 15),
#         plot.title = element_text(size=20),
#         legend.title=element_text(size=15),
#         axis.text.x =element_text(size=15,angle = 90),
#         axis.title=element_text(size=15)
#   ) + scale_y_continuous(breaks = seq(0,700,50),limits = c(0,700)) + #scale_fill_brewer(palette = "Set2") +
#   ggtitle("Total of running time across all MS/MS softwares") +
#   labs(fill="Wombat pipeline steps",x= "Task Name", y= "Running Time (min)") +
#   scale_fill_manual(breaks=c("COMPOMICS","PROLINE","TPP","MAXQUANT",
#                              "CALCBENCHMARKS_COMPOMICS","CALCBENCHMARKS_PROLINE",
#                              "CALCBENCHMARKS_TPP","CALCBENCHMARKS_MQ",
#                              "PREPARE_FILES","SDRFMERGE"),
#                     values = c("#E69F00","#56B4E9","#009E73","#CC79A7",
#                                "#999999","#999999","#999999","#999999",
#                                "#000000","#000000","#000000","#000000"
#                     )) + scale_x_discrete(limits = c("wombatp_prep",
#                                                      "peaklist_creation",
#                                                      "MSMS_search",
#                                                      "quantification",
#                                                      "statistical_testing",
#                                                      "writing_tables",
#                                                      "other",
#                                                      "benchmark"))
# 



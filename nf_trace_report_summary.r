library("dplyr")
library("neonUtilities")
library("stringr")
library("kableExtra")
library("reshape2")
library("lubridate")
raw_file_elixir <- read.csv("compomics_dev-20210111-VS_DK_VS_DK.txt", sep="\t")

wombat_dir <- setwd("D:/dev/R/WombatP-Utilities/data/")
wombat_dir <- setwd("D:/dev/R/WombatP-Utilities/data/")

list_wombat <- list.files(wombat_dir)[-3]
all_final <- NULL
contain_time <-c("duration","realtime")
time_type <- c("h","m","s")

for (a in 1:length(list_wombat)){
  assign(paste0("raw_file_elixir",a),read.delim(paste0(list_wombat[a]),header=TRUE, sep = "\t"))
  #
  assign(paste0("raw_file_elixir"), get(paste0("raw_file_elixir",a)))
  submit_conv <- as.data.frame(strptime(raw_file_elixir$submit,format = "%Y-%M-%d  %H:%M:%OS"))
  colnames(submit_conv) <-"submit"
  #
  for (d in 1:length(contain_time)){
    temp1<-NULL
    temp_h <- str_extract_all(raw_file_elixir[,d+7], time_type[1])
    temp_m <- str_extract_all(raw_file_elixir[,d+7], time_type[2])
    temp_s <- str_extract_all(raw_file_elixir[,d+7], time_type[3])
    for (e in 1:dim(raw_file_elixir)[1]){
      if (is.element("h",unlist(temp_h[[e]]))){
        
        query_hms <- str_match(raw_file_elixir[,d+7][e],'(\\d+)h\\s(\\d+)m\\s(\\d+)s')
        temp1[e] <- as.numeric(query_hms[2]) *60*60 + as.numeric(query_hms[3])*60 + as.numeric(query_hms[4])}
      
      else if ((is.element("m",unlist(temp_m[[e]])))){
        query_ms <- str_match(raw_file_elixir[,d+7][e],'(\\d+)m\\s(\\d+)s')
        temp1[e] <- as.numeric(query_ms[2]) *60 + as.numeric(query_ms[3])}
      else if ((is.element("s",unlist(temp_s[[e]]))) || 
               (!is.element("m",unlist(temp_m[[e]])))){
        query_s <- str_extract(raw_file_elixir[,d+7][e],'[[:digit:]]+\\.*[[:digit:]]*')
        temp1[e] <- as.numeric(query_s)
      }
      
      assign(paste0(contain_time[d],"_conv"),as.data.frame(temp1))
      
    }
    rm(temp1)
  }
  
    
  contain_bytes <-c("peak_rss","peak_vmem","rchar","wchar")
  byte_type <- c("GB","MB","KB","B")
  
  for (j in 1:length(contain_bytes)){
    
    #to prevent NA in the index list of each contain_bytes column
    raw_file_elixir[,j+10][raw_file_elixir[,j+10] == 0] <- "0 B"
    
    df <- assign(paste0("sep_",contain_bytes[j]),
                 as.data.frame(t(as.data.frame(sapply(strsplit(raw_file_elixir[,j+10], " "),
                                                      "[", 1:2)))))
    colnames(df) <- c("num","MB_name")
    temp2 <- NULL
    for (d in 1:length(byte_type)){
      
      indx_list <- str_extract_all(df$MB_name, paste0("^",byte_type[d]))
      for (c in 1:(dim(raw_file_elixir)[1])){
        
        if (is.element("GB",unlist(indx_list[[c]]))){
          
          temp2[c] <- as.numeric(df$num[c])*1000000000
          
        }else if (is.element("MB",unlist(indx_list[[c]]))){
          
          temp2[c] <- as.numeric(df$num[c])*1000000
          
        }else if(is.element("KB",unlist(indx_list[[c]]))){
          
          temp2[c] <- as.numeric(df$num[c])*1000
          
        }else if(is.element("B",unlist(indx_list[[c]]))){ 
          
          temp2[c] <- as.numeric(df$num[c])
        }
        
      }
      
    }
    assign(paste0(contain_bytes[j],"_gb"),temp2)
    rm(temp2)
  }  
  
  
  xcpu <- gsub("%","",raw_file_elixir$X.cpu)
  xcpu_tmp <- as.numeric(unlist(xcpu))
  
  rmv_num_name <- as.data.frame(sapply(strsplit(raw_file_elixir[,"name"], " "), "[", 1))
  
  refined_raw_file_elixir <- cbind(rmv_num_name[,1],duration_conv$temp1,
                                   xcpu_tmp,peak_rss_gb,rchar_gb, wchar_gb,as.data.frame.POSIXct(submit_conv$submit))
  refined_raw_file_elixir[is.na(refined_raw_file_elixir)] <-0 
  
  colnames(refined_raw_file_elixir) <-  c(colnames(raw_file_elixir[c(4,8,10,11,13,14)]),"submit")
  col_name <- colnames(refined_raw_file_elixir)
  
  final1 <- NULL
  final2 <- NULL
  final3 <- NULL
  
  for (i in 3:length(col_name)){

    if (col_name[i] == "submit"){
      
      tmpd <- assign(paste0("total_duration",col_name[i]), as.data.frame(refined_raw_file_elixir) %>% 
                       group_by(name) %>% summarise(duration=max(submit)-min(submit) + max(duration(duration,"seconds"))))
      tmpd <- cbind(col_name[i],tmpd,paste0(list_wombat[a]))
      final1 <- bind_rows(final1,tmpd)
      
    }else if (col_name[i] == "peak_rss"){
      tmpps <- assign(paste0("stas_calc_",col_name[i]), as.data.frame(refined_raw_file_elixir) %>%
                        group_by(name) %>% summarise(mean=mean(get(col_name[i]))))
      tmpps <- cbind(col_name[i],tmpps,paste0(list_wombat[a]))
      final2 <- bind_rows(final2,tmpps)
      
    }else{
      tmprest <- assign(paste0("stas_calc_",col_name[i]), as.data.frame(refined_raw_file_elixir) %>%
                          group_by(name) %>% summarise(sum=sum(get(col_name[i]))))
      tmprest <- cbind(col_name[i],tmprest,paste0(list_wombat[a]))
      final3 <- bind_rows(final3,tmprest) 
      
      # These two lines create a separate object that includes the sum,avg, and std of each column.
      #colnames(tmpd)[1] <- paste0("name_",col_name[i])
      #assign(paste0("stas_calc_",col_name[i]),d)
      
      #tmpd <- cbind(col_name[i],tmpd,paste0(list_wombat[a]))
      #final <- cbind2(final1,final2,final3)
      #colnames(tmpd)[6] <- "software_name"
    }
    
  }
  rm(list = ls()[!(ls() %in% c('all_final','list_wombat','final1','final2','final3',
                               'contain_time','contain_byte','time_type','col_name'))])
  all_final <- bind_rows(all_final,final1,final2,final3)
  
}

colnames(all_final) <- c("trace_report","task_name","real_duration", "software_name","mean","sum")
write.table(final, file = "all_list_duration_real_time_cpu_avg_std_sum_new.txt", sep = "\t")

test <- all_final %>% group_by(trace_report,software_name) %>% filter(trace_report == "submit")



# final %>%
#   kbl(caption = "Summary of raw file from Elixir") %>%
#   kable_classic(full_width = T, html_font = "Cambria")
all_final %>% kbl(caption = "Summary of raw file from Elixir") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))



melt <- melt(all_final)
colnames(melt)[c(1,3)] <- c("col_name","software_id")
melt1 <- cbind(melt[,c(2,3,5)],as.data.frame(str_c(melt$col_name,"_",melt$variable)))
colnames(melt1)[4] <- c("col_name_calc_type") 
#melt_soft_id <- melt %>% group_by(software_id) %>% summarise()
library(ggplot2)


library("ggplot2")
ggplot(melt) +
  geom_bar( aes(x=variable, y=value,fill=col_name[i]), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=name, ymin=value-sd(value), ymax=value+sd(value)), width=0.4, colour="orange", alpha=0.9, size=1.3)

final_melt <- melt %>% group_by(col_name,variable) %>% filter(name == "convert_raw_mzml")

ggplot(test, aes(fill=software_name, y=real_duration, x=task_name)) + 
  geom_bar(position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


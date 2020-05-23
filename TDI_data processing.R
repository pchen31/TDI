# Load packages 
pacman::p_load(pacman, dplyr, tidyr, stringr, readxl, writexl, outliers)

# Rawdata file name
Name_rawdata_file <- list.files(pattern = "*.xlsx", full.names = T)

# Load raw data
Rawdata_3A4 <- read_excel(Name_rawdata_file, sheet = "3A4_Midazolam_rawdata")                       
Rawdata_2C8 <- read_excel(Name_rawdata_file, sheet = "2C8_Paclitaxel_rawdata")
Rawdata_2C9 <- read_excel(Name_rawdata_file, sheet = "2C9_Diclofenac_rawdata")
Rawdata_2D6 <- read_excel(Name_rawdata_file, sheet = "2D6_Dextromethorphan_rawdata")
Df_Cpd_List <- read_excel(Name_rawdata_file, sheet = "Compound List")
#View(Df_Cpd_List)

# Clean up compound list
Df_Cpd_List <- Df_Cpd_List[6:37, 2:5]
colnames(Df_Cpd_List) <- c("ID", "Project", "Batch", "Concentration (uM)")

################ Define variables ########################
Val_replicate <- 3                                                                                  # 3 for triplicates
Val_pct_CV_cutoff <- 20                                                                             # Cut-off value of %CV  
Val_inc_time = 30

Df_Cpd_List <- mutate(Df_Cpd_List,
                      "Structure" = "",
                      "Pre-inc time (min)" = Val_inc_time)
Df_Cpd_List <- Df_Cpd_List[c(32, 31, 1:30), c(1,3,2, 4:6)]

# Remove MeOH wash samples
Rawdata_3A4 <- Rawdata_3A4[!grepl("MeOH", Rawdata_3A4$`Sample Name`), ]                     
Rawdata_2C8 <- Rawdata_2C8[!grepl("MeOH", Rawdata_2C8$`Sample Name`), ]                     
Rawdata_2C9 <- Rawdata_2C9[!grepl("MeOH", Rawdata_2C9$`Sample Name`), ]     
Rawdata_2D6 <- Rawdata_2D6[!grepl("MeOH", Rawdata_2D6$`Sample Name`), ]       

# Split strings in sample name
Rawdata_3A4 <- Rawdata_3A4 %>% separate(`Sample Name`, c("ID", "NADPH", "Replicate"), "_")
Rawdata_2C8 <- Rawdata_2C8 %>% separate(`Sample Name`, c("ID", "NADPH", "Replicate"), "_")
Rawdata_2C9 <- Rawdata_2C9 %>% separate(`Sample Name`, c("ID", "NADPH", "Replicate"), "_")
Rawdata_2D6 <- Rawdata_2D6 %>% separate(`Sample Name`, c("ID", "NADPH", "Replicate"), "_")

# Calculate PAR (Analyte to ISTD peak area ratio)
Rawdata_3A4 <- mutate(Rawdata_3A4, PAR = Rawdata_3A4$`Analyte Peak Area (counts)`/Rawdata_3A4$`IS Peak Area (counts)`)
Rawdata_2C8 <- mutate(Rawdata_2C8, PAR = Rawdata_2C8$`Analyte Peak Area (counts)`/Rawdata_2C8$`IS Peak Area (counts)`)
Rawdata_2C9 <- mutate(Rawdata_2C9, PAR = Rawdata_2C9$`Analyte Peak Area (counts)`/Rawdata_2C9$`IS Peak Area (counts)`)
Rawdata_2D6 <- mutate(Rawdata_2D6, PAR = Rawdata_2D6$`Analyte Peak Area (counts)`/Rawdata_2D6$`IS Peak Area (counts)`)

# Split (+) and (-) NADPH samples
Rawdata_3A4_plus <- filter(Rawdata_3A4, Rawdata_3A4$NADPH == "(+)NADPH")
Rawdata_3A4_minus <- filter(Rawdata_3A4, Rawdata_3A4$NADPH == "(-)NADPH")
Rawdata_2C8_plus <- filter(Rawdata_2C8, Rawdata_2C8$NADPH == "(+)NADPH")
Rawdata_2C8_minus <- filter(Rawdata_2C8, Rawdata_2C8$NADPH == "(-)NADPH")
Rawdata_2C9_plus <- filter(Rawdata_2C9, Rawdata_2C9$NADPH == "(+)NADPH")
Rawdata_2C9_minus <- filter(Rawdata_2C9, Rawdata_2C9$NADPH == "(-)NADPH")
Rawdata_2D6_plus <- filter(Rawdata_2D6, Rawdata_2D6$NADPH == "(+)NADPH")
Rawdata_2D6_minus <- filter(Rawdata_2D6, Rawdata_2D6$NADPH == "(-)NADPH")

# Define variables
# Num_sample_3A4_plus <- nrow(Rawdata_3A4_plus)                                                     # Number of samples per CYP isoform w or w/o NADPH
# List_num_by_replicate <- seq(1, Num_sample_3A4_plus, Val_replicate)                    
# List_cpd <- Rawdata_3A4_plus$ID[List_num_by_replicate]                                            # Get compound list to use in final results
# Val_cpd_list <- Num_sample_3A4_plus/Val_replicate   

# Create function to calculate avg, sd and %CV without the value most differing from mean
mean.rm.o <- function(Data) { zval_1 <- rm.outlier(Data)
                              zval_2 <- mean(zval_1, na.rm = TRUE)
                              return(zval_2)} 

pct_CV <- function(Data) { zval_1 <- mean(Data, na.rm = TRUE)
                           zval_2 <- sd(Data, na.rm = TRUE)
                           zval_3 <- zval_2/zval_1*100
                           zval_4 <- round(zval_3, 1)
                           return(zval_4)} 

pct_CV.rm.o <- function(Data) { zval_1 <- rm.outlier(Data)
                                zval_2 <- mean(zval_1, na.rm = TRUE)
                                zval_3 <- sd(zval_1, na.rm = TRUE)
                                zval_4 <- zval_3/zval_2*100
                                zval_5 <- round(zval_4, 1)
                                return(zval_5)} 


######################### Calculate % activity remaining for 3A4 (+) NADPH samples #############################
Df_CYP3A4_plus <- Rawdata_3A4_plus %>%                                                                   
                  group_by(ID) %>%
                  summarize(PAR_avg = mean(PAR, na.rm = TRUE),
                            Pct_CV = pct_CV(PAR),
                            PAR_avg_rm_o = mean.rm.o(PAR),
                            Pct_CV_rm_o = pct_CV.rm.o(PAR))

Val_row_num_DMSO <- str_which(Df_CYP3A4_plus$ID, "DMSO")                                         # find the row number containing "DMSO" in sample name

Df_CYP3A4_plus <- Df_CYP3A4_plus %>%                                                                           
                      mutate("PAR_avg (+)" = ifelse(Df_CYP3A4_plus$Pct_CV < Val_pct_CV_cutoff, Df_CYP3A4_plus$PAR_avg, 
                                                ifelse(Df_CYP3A4_plus$Pct_CV_rm_o < Val_pct_CV_cutoff, Df_CYP3A4_plus$PAR_avg_rm_o, NA))
                             )
Df_CYP3A4_plus <- Df_CYP3A4_plus %>% mutate("PAR_avg_DMSO (+)" = Df_CYP3A4_plus$`PAR_avg (+)`[Val_row_num_DMSO])
Df_CYP3A4_plus <- Df_CYP3A4_plus %>% mutate("3A4 % Activity remaining (+)" = round(Df_CYP3A4_plus$`PAR_avg (+)`/Df_CYP3A4_plus$`PAR_avg_DMSO (+)`*100, 1))


######################### Calculate % activity remaining for 3A4 (-) NADPH samples #############################
Df_CYP3A4_minus <- Rawdata_3A4_minus %>%                                                                   
  group_by(ID) %>%
  summarize(PAR_avg = mean(PAR, na.rm = TRUE),
            Pct_CV = pct_CV(PAR),
            PAR_avg_rm_o = mean.rm.o(PAR),
            Pct_CV_rm_o = pct_CV.rm.o(PAR))

Val_row_num_DMSO <- str_which(Df_CYP3A4_minus$ID, "DMSO")                                         # find the row number containing "DMSO" in sample name

Df_CYP3A4_minus <- Df_CYP3A4_minus %>%                                                                           
  mutate("PAR_avg (-)" = ifelse(Df_CYP3A4_minus$Pct_CV < Val_pct_CV_cutoff, Df_CYP3A4_minus$PAR_avg, 
                                ifelse(Df_CYP3A4_minus$Pct_CV_rm_o < Val_pct_CV_cutoff, Df_CYP3A4_minus$PAR_avg_rm_o, NA))
  )
Df_CYP3A4_minus <- Df_CYP3A4_minus %>% mutate("PAR_avg_DMSO (-)" = Df_CYP3A4_minus$`PAR_avg (-)`[Val_row_num_DMSO])
Df_CYP3A4_minus <- Df_CYP3A4_minus %>% mutate("3A4 % Activity remaining (-)" = round(Df_CYP3A4_minus$`PAR_avg (-)`/Df_CYP3A4_minus$`PAR_avg_DMSO (-)`*100, 1))


######################### Calculate % activity remaining for 2C8 (+) NADPH samples #############################
Df_CYP2C8_plus <- Rawdata_2C8_plus %>%                                                                   
  group_by(ID) %>%
  summarize(PAR_avg = mean(PAR, na.rm = TRUE),
            Pct_CV = pct_CV(PAR),
            PAR_avg_rm_o = mean.rm.o(PAR),
            Pct_CV_rm_o = pct_CV.rm.o(PAR))

Val_row_num_DMSO <- str_which(Df_CYP2C8_plus$ID, "DMSO")                                          # find the row number containing "DMSO" in sample name

Df_CYP2C8_plus <- Df_CYP2C8_plus %>%                                                                           
  mutate("PAR_avg (+)" = ifelse(Df_CYP2C8_plus$Pct_CV < Val_pct_CV_cutoff, Df_CYP2C8_plus$PAR_avg, 
                                ifelse(Df_CYP2C8_plus$Pct_CV_rm_o < Val_pct_CV_cutoff, Df_CYP2C8_plus$PAR_avg_rm_o, NA))
  )
Df_CYP2C8_plus <- Df_CYP2C8_plus %>% mutate("PAR_avg_DMSO (+)" = Df_CYP2C8_plus$`PAR_avg (+)`[Val_row_num_DMSO])
Df_CYP2C8_plus <- Df_CYP2C8_plus %>% mutate("2C8 % Activity remaining (+)" = round(Df_CYP2C8_plus$`PAR_avg (+)`/Df_CYP2C8_plus$`PAR_avg_DMSO (+)`*100, 1))


######################### Calculate % activity remaining for 2C8 (-) NADPH samples #############################
Df_CYP2C8_minus <- Rawdata_2C8_minus %>%                                                                   
  group_by(ID) %>%
  summarize(PAR_avg = mean(PAR, na.rm = TRUE),
            Pct_CV = pct_CV(PAR),
            PAR_avg_rm_o = mean.rm.o(PAR),
            Pct_CV_rm_o = pct_CV.rm.o(PAR))

Val_row_num_DMSO <- str_which(Df_CYP2C8_minus$ID, "DMSO")                                         # find the row number containing "DMSO" in sample name

Df_CYP2C8_minus <- Df_CYP2C8_minus %>%                                                                           
  mutate("PAR_avg (-)" = ifelse(Df_CYP2C8_minus$Pct_CV < Val_pct_CV_cutoff, Df_CYP2C8_minus$PAR_avg, 
                                ifelse(Df_CYP2C8_minus$Pct_CV_rm_o < Val_pct_CV_cutoff, Df_CYP2C8_minus$PAR_avg_rm_o, NA))
  )
Df_CYP2C8_minus <- Df_CYP2C8_minus %>% mutate("PAR_avg_DMSO (-)" = Df_CYP2C8_minus$`PAR_avg (-)`[Val_row_num_DMSO])
Df_CYP2C8_minus <- Df_CYP2C8_minus %>% mutate("2C8 % Activity remaining (-)" = round(Df_CYP2C8_minus$`PAR_avg (-)`/Df_CYP2C8_minus$`PAR_avg_DMSO (-)`*100, 1))


######################### Calculate % activity remaining for 2C9 (+) NADPH samples #############################
Df_CYP2C9_plus <- Rawdata_2C9_plus %>%                                                                   
  group_by(ID) %>%
  summarize(PAR_avg = mean(PAR, na.rm = TRUE),
            Pct_CV = pct_CV(PAR),
            PAR_avg_rm_o = mean.rm.o(PAR),
            Pct_CV_rm_o = pct_CV.rm.o(PAR))

Val_row_num_DMSO <- str_which(Df_CYP2C9_plus$ID, "DMSO")                                          # find the row number containing "DMSO" in sample name

Df_CYP2C9_plus <- Df_CYP2C9_plus %>%                                                                           
  mutate("PAR_avg (+)" = ifelse(Df_CYP2C9_plus$Pct_CV < Val_pct_CV_cutoff, Df_CYP2C9_plus$PAR_avg, 
                                ifelse(Df_CYP2C9_plus$Pct_CV_rm_o < Val_pct_CV_cutoff, Df_CYP2C9_plus$PAR_avg_rm_o, NA))
  )
Df_CYP2C9_plus <- Df_CYP2C9_plus %>% mutate("PAR_avg_DMSO (+)" = Df_CYP2C9_plus$`PAR_avg (+)`[Val_row_num_DMSO])
Df_CYP2C9_plus <- Df_CYP2C9_plus %>% mutate("2C9 % Activity remaining (+)" = round(Df_CYP2C9_plus$`PAR_avg (+)`/Df_CYP2C9_plus$`PAR_avg_DMSO (+)`*100, 1))


######################### Calculate % activity remaining for 2C9 (-) NADPH samples #############################
Df_CYP2C9_minus <- Rawdata_2C9_minus %>%                                                                   
  group_by(ID) %>%
  summarize(PAR_avg = mean(PAR, na.rm = TRUE),
            Pct_CV = pct_CV(PAR),
            PAR_avg_rm_o = mean.rm.o(PAR),
            Pct_CV_rm_o = pct_CV.rm.o(PAR))

Val_row_num_DMSO <- str_which(Df_CYP2C9_minus$ID, "DMSO")                                        # find the row number containing "DMSO" in sample name

Df_CYP2C9_minus <- Df_CYP2C9_minus %>%                                                                           
  mutate("PAR_avg (-)" = ifelse(Df_CYP2C9_minus$Pct_CV < Val_pct_CV_cutoff, Df_CYP2C9_minus$PAR_avg, 
                                ifelse(Df_CYP2C9_minus$Pct_CV_rm_o < Val_pct_CV_cutoff, Df_CYP2C9_minus$PAR_avg_rm_o, NA))
  )
Df_CYP2C9_minus <- Df_CYP2C9_minus %>% mutate("PAR_avg_DMSO (-)" = Df_CYP2C9_minus$`PAR_avg (-)`[Val_row_num_DMSO])
Df_CYP2C9_minus <- Df_CYP2C9_minus %>% mutate("2C9 % Activity remaining (-)" = round(Df_CYP2C9_minus$`PAR_avg (-)`/Df_CYP2C9_minus$`PAR_avg_DMSO (-)`*100, 1))


######################### Calculate % activity remaining for 2D6 (+) NADPH samples #############################
Df_CYP2D6_plus <- Rawdata_2D6_plus %>%                                                                   
  group_by(ID) %>%
  summarize(PAR_avg = mean(PAR, na.rm = TRUE),
            Pct_CV = pct_CV(PAR),
            PAR_avg_rm_o = mean.rm.o(PAR),
            Pct_CV_rm_o = pct_CV.rm.o(PAR))

Val_row_num_DMSO <- str_which(Df_CYP2D6_plus$ID, "DMSO")                                         # find the row number containing "DMSO" in sample name

Df_CYP2D6_plus <- Df_CYP2D6_plus %>%                                                                           
  mutate("PAR_avg (+)" = ifelse(Df_CYP2D6_plus$Pct_CV < Val_pct_CV_cutoff, Df_CYP2D6_plus$PAR_avg, 
                                ifelse(Df_CYP2D6_plus$Pct_CV_rm_o < Val_pct_CV_cutoff, Df_CYP2D6_plus$PAR_avg_rm_o, NA))
  )
Df_CYP2D6_plus <- Df_CYP2D6_plus %>% mutate("PAR_avg_DMSO (+)" = Df_CYP2D6_plus$`PAR_avg (+)`[Val_row_num_DMSO])
Df_CYP2D6_plus <- Df_CYP2D6_plus %>% mutate("2D6 % Activity remaining (+)" = round(Df_CYP2D6_plus$`PAR_avg (+)`/Df_CYP2D6_plus$`PAR_avg_DMSO (+)`*100, 1))


######################### Calculate % activity remaining for 2D6 (-) NADPH samples #############################
Df_CYP2D6_minus <- Rawdata_2D6_minus %>%                                                                   
  group_by(ID) %>%
  summarize(PAR_avg = mean(PAR, na.rm = TRUE),
            Pct_CV = pct_CV(PAR),
            PAR_avg_rm_o = mean.rm.o(PAR),
            Pct_CV_rm_o = pct_CV.rm.o(PAR))

Val_row_num_DMSO <- str_which(Df_CYP2D6_minus$ID, "DMSO")                                         # find the row number containing "DMSO" in sample name

Df_CYP2D6_minus <- Df_CYP2D6_minus %>%                                                                           
  mutate("PAR_avg (-)" = ifelse(Df_CYP2D6_minus$Pct_CV < Val_pct_CV_cutoff, Df_CYP2D6_minus$PAR_avg, 
                                ifelse(Df_CYP2D6_minus$Pct_CV_rm_o < Val_pct_CV_cutoff, Df_CYP2D6_minus$PAR_avg_rm_o, NA))
  )
Df_CYP2D6_minus <- Df_CYP2D6_minus %>% mutate("PAR_avg_DMSO (-)" = Df_CYP2D6_minus$`PAR_avg (-)`[Val_row_num_DMSO])
Df_CYP2D6_minus <- Df_CYP2D6_minus %>% mutate("2D6 % Activity remaining (-)" = round(Df_CYP2D6_minus$`PAR_avg (-)`/Df_CYP2D6_minus$`PAR_avg_DMSO (-)`*100, 1))


######################### Final summary file #########################

Processed_data <- do.call("cbind", list(Df_CYP3A4_plus[,c(1,8)], Df_CYP3A4_minus[8], Df_CYP2C8_plus[8], Df_CYP2C8_minus[8], 
                                       Df_CYP2C9_plus[8], Df_CYP2C9_minus[8], Df_CYP2D6_plus[8], Df_CYP2D6_minus[8]))

Processed_data <- Processed_data %>% mutate("3A4 % Activity Loss" = `3A4 % Activity remaining (-)` - `3A4 % Activity remaining (+)`,
                                          "2C8 % Activity Loss" = `2C8 % Activity remaining (-)` - `2C8 % Activity remaining (+)`,
                                          "2C9 % Activity Loss" = `2C9 % Activity remaining (-)` - `2C9 % Activity remaining (+)`,
                                          "2D6 % Activity Loss" = `2D6 % Activity remaining (-)` - `2D6 % Activity remaining (+)`)

Processed_data$`3A4 % Activity Loss` <- pmax(Processed_data$`3A4 % Activity Loss`,0)
Processed_data$`2C8 % Activity Loss` <- pmax(Processed_data$`2C8 % Activity Loss`,0)
Processed_data$`2C9 % Activity Loss` <- pmax(Processed_data$`2C9 % Activity Loss`,0)
Processed_data$`2D6 % Activity Loss` <- pmax(Processed_data$`2D6 % Activity Loss`,0)

Processed_data <- Processed_data[, c(1:3, 10, 4:5, 11, 6:7, 12, 8:9, 13)]

# Prepare dataframe for summary file
Data_Summary <- left_join(Df_Cpd_List, Processed_data, "ID") 

# Prepare CDD upload file
Num_cpd <- nrow(Df_CYP3A4_plus)
CDD_upload_3A4 <- left_join(Df_Cpd_List[3:Num_cpd, c(1:2,4,6)], Processed_data[, c(1:4)], "ID")
CDD_upload_2C8 <- left_join(Df_Cpd_List[3:Num_cpd, c(1:2,4,6)], Processed_data[, c(1,5:7)], "ID")
CDD_upload_2C9 <- left_join(Df_Cpd_List[3:Num_cpd, c(1:2,4,6)], Processed_data[, c(1,8:10)], "ID")
CDD_upload_2D6 <- left_join(Df_Cpd_List[3:Num_cpd, c(1:2,4,6)], Processed_data[, c(1,11:13)], "ID")

colnames(CDD_upload_3A4) <- c("Molecule_Name", "Batch", "Concentration (µM)", "Pre_Incubation_Time (min)", "Pct_Plus_NADPH", "Pct_Minus_NADPH", "Pct_Activity_Loss")
colnames(CDD_upload_2C8) <- c("Molecule_Name", "Batch", "Concentration (µM)", "Pre_Incubation_Time (min)", "Pct_Plus_NADPH", "Pct_Minus_NADPH", "Pct_Activity_Loss")
colnames(CDD_upload_2C9) <- c("Molecule_Name", "Batch", "Concentration (µM)", "Pre_Incubation_Time (min)", "Pct_Plus_NADPH", "Pct_Minus_NADPH", "Pct_Activity_Loss")
colnames(CDD_upload_2D6) <- c("Molecule_Name", "Batch", "Concentration (µM)", "Pre_Incubation_Time (min)", "Pct_Plus_NADPH", "Pct_Minus_NADPH", "Pct_Activity_Loss")

CDD_upload_all <- bind_rows(CDD_upload_3A4, CDD_upload_2C8)
CDD_upload_all <- bind_rows(CDD_upload_all, CDD_upload_2C9)
CDD_upload_all <- bind_rows(CDD_upload_all, CDD_upload_2D6)
CDD_upload_all <- mutate(CDD_upload_all, 
                         "Enzyme (eg 3A4, 2D6, 1A2)" = rep(c("3A4", "2C8", "2C9", "2D6"), each = Num_cpd - 2),
                         "Comment" = "")
CDD_upload_all <- CDD_upload_all[ ,c(1:2, 8, 3:7, 9)]

# Export final results to excel files
Val_current_date <- Sys.Date()                                                                                # Get the current date to attach to the file name

List_Summary <- list("Summary" = Data_Summary)
write_xlsx(List_Summary,                                                                                      # Export data for summary to an excel file
           paste(Val_current_date, " CYP TDI 30min assay - 3A4 2C8 2C9 2D6 - Summary", ".xlsx", sep = ""))  

List_CDD_upload <- list("DatabaseResults" = CDD_upload_all)
write_xlsx(List_CDD_upload,                                                                                   # Export data for summary to an excel file
           paste(Val_current_date, " CYP TDI 30min assay - 3A4 2C8 2C9 2D6 - CDD upload", ".xlsx", sep = ""))  

##### CLEAN UP ######

rm(list = ls())   # Clear environment
p_unload(all)     # Remove all add-ons
cat("\014")       # ctrl+L # Clear console
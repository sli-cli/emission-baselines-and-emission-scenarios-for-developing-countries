####################################################################################################
rm(list = ls())
cat("\14")

#====================================================================================================
wdir='F:/research/1668008131197-AR6_Scenarios_Database_ISO3_v1.1.csv' # set working directory
setwd(paste(wdir))
library(reshape2)
library(openxlsx)
library(writexl)
library(dplyr)
library(abind)
#====================================================================================================      
#====================================================================================================      
#---Read raw metadata file---meta_Ch3vetted_with_climate
temp1=openxlsx::read.xlsx('./AR6_Scenarios_Database_metadata_indicators_v1.1.xlsx',sheet=2,colNames=TRUE)
#---Subset C1-C8 scenarios
ms = temp1[which(temp1$Category %in% c("C1",'C3','C5','C6','C7','C8')), c(1:4,8,12)]; dim(ms)
ms$ms = paste(ms$Model, ms$Scenario, sep='--') 
ms %>% dplyr::count(Category)
#====================================================================================================   
#---Add new "Model_family" categorization
ms$Model_family = "Others"
ms[which(ms$Model %in% c("AIM/CGE 2.0","AIM/CGE 2.1","AIM/CGE 2.2", "AIM/Hub-Global 2.0")), ncol(ms)] = "AIM"
ms[which(ms$Model %in% c("C-ROADS-5.005")), ncol(ms)] = "C-ROADS"
ms[which(ms$Model %in% c("COFFEE 1.1")), ncol(ms)] = "COFFEE"
ms[which(ms$Model %in% c("EPPA 6")), ncol(ms)] = "EPPA"
ms[which(ms$Model %in% c("GCAM-PR 5.3", "GCAM 4.2","GCAM 5.3","GCAM 5.2")), ncol(ms)] = "GCAM"
ms[which(ms$Model %in% c("GEM-E3_V2021")), ncol(ms)] = "GEM-E3"
ms[which(ms$Model %in% c("IMAGE 3.0","IMAGE 3.0.1","IMAGE 3.0.2","IMAGE 3.2")), ncol(ms)] = "IMAGE"
ms[which(ms$Model %in% c("MESSAGE-GLOBIOM 1.0","MESSAGEix-GLOBIOM 1.0","MESSAGEix-GLOBIOM_1.1",
                         "MESSAGEix-GLOBIOM_1.2","MESSAGEix-GLOBIOM_GEI 1.0")), ncol(ms)] = "MESSAGE"
ms[which(ms$Model %in% c("POLES ADVANCE","POLES EMF30","POLES CD-LINKS","POLES EMF33","POLES ENGAGE","POLES GECO2019")), ncol(ms)] = "POLES"
ms[which(ms$Model %in% c("REMIND 1.6","REMIND 1.7","REMIND 2.1","REMIND-Buildings 2.0","REMIND-MAgPIE 1.5","REMIND-MAgPIE 1.7-3.0",
                         "REMIND-MAgPIE 2.0-4.1","REMIND-MAgPIE 2.1-4.2","REMIND-MAgPIE 2.1-4.3","REMIND-Transport 2.1")), ncol(ms)] = "REMIND"
ms[which(ms$Model %in% c("TIAM-ECN 1.1")), ncol(ms)] = "TIAM"
ms[which(ms$Model %in% c("WITCH 4.6","WITCH 5.0","WITCH-GLOBIOM 3.1","WITCH-GLOBIOM 4.2","WITCH-GLOBIOM 4.4")), ncol(ms)] = "WITCH"
#which(ms$Model_family=="Others")
#ms <- ms[!(ms$Model_family %in% c('COFFEE', 'GEM-E3', 'TIAM', 'WITCH')), ]
length(unique(ms$Model_family)) # 12 model families

ms %>% dplyr::count(Category,Model_family)
ms %>% dplyr::count(Category, Project_study)
wdir='F:/research/code_and_data' # set working directory
setwd(paste(wdir))
save(ms, file=paste(wdir,'./data_generation/AR6-v1.1_C1-C8_metadata.Rdata',sep=''))

#====================================================================================================  
####################################################################################################
rm(list = ls())
cat("\14")

#====================================================================================================
wdir='F:/research/1668008131197-AR6_Scenarios_Database_ISO3_v1.1.csv' 
setwd(paste(wdir))
library(reshape2)
library(openxlsx)
library(writexl)
library(plyr)
library(abind)
library(zoo)
#====================================================================================================      
#====================================================================================================  
#--- Load previously compiled metadata——ISO
load('F://research/code_and_data/data_generation/AR6-v1.1_C1-C8_metadata.Rdata') 
ind = which(colnames(ms) %in% c('ms','Category','Model_family','Project_study','Scenario_project','IMP_marker'))
ms2 = ms[,ind]
temp = read.csv('./AR6_Scenarios_Database_ISO3_v1.1.csv', header=TRUE) 
dim(temp)
colnames(temp)
unique(temp$Region) 
vars=unique(temp$Variable)
temp$ms = paste(temp$Model, temp$Scenario, sep='--')

data=temp[which(temp$ms %in% ms$ms),]
data = merge(ms2, data, by='ms')
colnames(data) <- gsub("^X", "", colnames(data))
wdir='F:/research/code_and_data' # set working directory
setwd(paste(wdir))
save(data, file='./data_generation/AR6-v1.1_C1-C8_allvars.Rdata')

#====================================================================================================      
#--- Select Emissions|CO2|Energy——ISO
load('./data_generation/AR6-v1.1_C1-C8_allvars.Rdata') # data
sv=sort(c("Emissions|CO2|Energy"))

data=data[which(data$Variable %in% sv),]
years=c(2010:2100)
itv=seq(2010,2100,by=5)
ind=which(colnames(data) %in% itv)
ind2=which(colnames(data)=="Unit")
data=data[,c(1:ind2,ind)]


write.csv(data,'./data_generation/ISO_ar6_mainmodel.csv',row.names = F)
#==================================================================================================== 
#==================================================================================================== 
rm(list = ls())
cat("\14")

#====================================================================================================
wdir='F:/research/1668008228539-AR6_Scenarios_Database_R5_regions_v1.1.csv' 
setwd(paste(wdir))
library(reshape2)
library(openxlsx)
library(writexl)
library(plyr)
library(abind)
library(zoo)
#====================================================================================================      
#====================================================================================================      
#--- Load previously compiled metadata——R5
load('F://research/code_and_data/data_generation/AR6-v1.1_C1-C8_metadata.Rdata') 
ind = which(colnames(ms) %in% c('ms','Category','Model_family','Project_study','Scenario_project','IMP_marker'))
ms2 = ms[,ind]

temp = read.csv('./AR6_Scenarios_Database_R5_regions_v1.1.csv', header=TRUE) 
dim(temp)
colnames(temp)
unique(temp$Region) 
vars=unique(temp$Variable) 
temp$ms = paste(temp$Model, temp$Scenario, sep='--')

data=temp[which(temp$ms %in% ms$ms),]
data = merge(ms2, data, by='ms')
colnames(data) <- gsub("^X", "", colnames(data))
wdir='F:/research/code_and_data' # set working directory
setwd(paste(wdir))
save(data, file='./data_generation/R5/AR6-v1.1_C1-C8_allvars.Rdata')

#====================================================================================================      
#--- Select Emissions|CO2|Energy——R5
load('./data_generation/R5/AR6-v1.1_C1-C8_allvars.Rdata') 
sv=sort(c("Emissions|CO2|Energy"))

data=data[which(data$Variable %in% sv),]

years=c(2010:2100)
itv=seq(2010,2100,by=5)
ind=which(colnames(data) %in% itv)
ind2=which(colnames(data)=="Unit")
data=data[,c(1:ind2,ind)]

write.csv(data,'./data_generation/R5_ar6_mainmodel.csv',row.names = F)
# Max Min Median Figure 1------------------------------------------------------
rm(list = ls())
cat("\14")
library(reshape2)
library(openxlsx)
library(writexl)
library(plyr)
library(abind)
library(zoo)
library(readxl)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)
library(patchwork)

wdir='F:/research/code_and_data/data_generation' 
allIAMPath = "F:/research/code_and_data/data_generation/emissions_total.xlsx"
allIAM = read_excel(path = allIAMPath, sheet = 'origin_models')

grouped_data <- allIAM %>%
  group_by(Region, Model_family)

a <- which(colnames(grouped_data) == '2010')
b <- which(colnames(grouped_data) == '2020')
selected_columns <- c(colnames(grouped_data)[a:b])

summary_stats_by_group <- grouped_data %>%
  summarize(across(all_of(selected_columns), 
                   list(Q1 = ~ quantile(., 0.25, na.rm = TRUE),
                        Min = ~ min(., na.rm = TRUE),
                        Max = ~ max(., na.rm = TRUE),
                        Median = ~ median(., na.rm = TRUE),
                        Q3 = ~ quantile(., 0.75, na.rm = TRUE))))

median <- summary_stats_by_group %>%
  select(1:2, contains("_median"))

median_long <- median %>%
  pivot_longer(cols = starts_with("20"),  
               names_to = "year",           
               values_to = "value")        
median_long <- median_long %>%
  separate(col = year, into = c("year", "median"), sep = "_", remove = TRUE) %>%
  select(-median)  
median_long <- median_long %>%
  pivot_wider(names_from = Model_family,  
              values_from = value)   
median_long <- median_long %>%
  select(Region, year, AIM, GCAM, IMAGE, MESSAGE, POLES,REMIND,TIAM,WITCH,COFFEE,EPPA,'GEM-E3')

Generate_Max <- summary_stats_by_group %>%
  select(1:2, contains("_Max"))

Generate_Max_long <- Generate_Max %>%
  pivot_longer(cols = starts_with("20"),  
               names_to = "year",          
               values_to = "value")        
Generate_Max_long <- Generate_Max_long %>%
  separate(col = year, into = c("year", "Max"), sep = "_", remove = TRUE) %>%
  select(-Max)  
Generate_Max_long <- Generate_Max_long %>%
  pivot_wider(names_from = Model_family,  
              values_from = value)        
Generate_Max_long <- Generate_Max_long %>%
  select(Region, year, AIM, GCAM, IMAGE, MESSAGE, POLES,REMIND,TIAM,WITCH,COFFEE,EPPA,'GEM-E3')

Generate_Min <- summary_stats_by_group %>%
  select(1:2, contains("_Min"))

Generate_Min_long <- Generate_Min %>%
  pivot_longer(cols = starts_with("20"),  
               names_to = "year",           
               values_to = "value")        
Generate_Min_long <- Generate_Min_long %>%
  separate(col = year, into = c("year", "Min"), sep = "_", remove = TRUE) %>%
  select(-Min)  
Generate_Min_long <- Generate_Min_long %>%
  pivot_wider(names_from = Model_family,  
              values_from = value)        
Generate_Min_long <- Generate_Min_long %>%
  select(Region, year, AIM, GCAM, IMAGE, MESSAGE, POLES,REMIND,TIAM,WITCH,COFFEE,EPPA,'GEM-E3')

Generate_Q1 <- summary_stats_by_group %>%
  select(1:2, contains("_Q1"))

Generate_Q1_long <- Generate_Q1 %>%
  pivot_longer(cols = starts_with("20"),  
               names_to = "year",           
               values_to = "value")        
Generate_Q1_long <- Generate_Q1_long %>%
  separate(col = year, into = c("year", "Q1"), sep = "_", remove = TRUE) %>%
  select(-Q1)  
Generate_Q1_long <- Generate_Q1_long %>%
  pivot_wider(names_from = Model_family,  
              values_from = value)        
Generate_Q1_long <- Generate_Q1_long %>%
  select(Region, year, AIM, GCAM, IMAGE, MESSAGE, POLES,REMIND,TIAM,WITCH,COFFEE,EPPA,'GEM-E3')

Generate_Q3 <- summary_stats_by_group %>%
  select(1:2, contains("_Q3"))

Generate_Q3_long <- Generate_Q3 %>%
  pivot_longer(cols = starts_with("20"),  
               names_to = "year",           
               values_to = "value")        
Generate_Q3_long <- Generate_Q3_long %>%
  separate(col = year, into = c("year", "Q3"), sep = "_", remove = TRUE) %>%
  select(-Q3)  
Generate_Q3_long <- Generate_Q3_long %>%
  pivot_wider(names_from = Model_family,  
              values_from = value)        
Generate_Q3_long <- Generate_Q3_long %>%
  select(Region, year, AIM, GCAM, IMAGE, MESSAGE, POLES,REMIND,TIAM,WITCH,COFFEE,EPPA,'GEM-E3')
write.csv(Generate_Min_long,'F:/research/code_and_data/data_generation/figure1_min.csv',row.names = F)
write.csv(Generate_Max_long,'F:/research/code_and_data/data_generation/figure1_max.csv',row.names = F)
write.csv(median_long,'F:/research/code_and_data/data_generation/figure1_median.csv',row.names = F)
write.csv(Generate_Q3_long,'F:/research/code_and_data/data_generation/figure1_Q3.csv',row.names = F)
write.csv(Generate_Q1_long,'F:/research/code_and_data/data_generation/figure1_Q1.csv',row.names = F)

# change Figure2 2010 2015 2020------------------------------------------------------
rm(list = ls())
cat("\14")
library(reshape2)
library(openxlsx)
library(writexl)
library(plyr)
library(abind)
library(zoo)
library(readxl)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)
library(patchwork)

wdir='F:/research/code_and_data/data_generation' 
allIAMPath = "F:/research/code_and_data/data_generation/emissions_total.xlsx"
allIAM = read_excel(path = allIAMPath, sheet = 'origin_models')
allIAM$country = substr(allIAM$Region, nchar(allIAM$Region)-2, nchar(allIAM$Region))


grouped_data <- allIAM %>%
  group_by(country, Model_family)


a <- which(colnames(grouped_data) == '2010')
b <- which(colnames(grouped_data) == '2020')
selected_columns <- c(colnames(grouped_data)[a:b])

summary_stats_by_group <- grouped_data %>%
  summarize(across(all_of(selected_columns), 
                   list(Q1 = ~ quantile(., 0.25, na.rm = TRUE),
                        Min = ~ min(., na.rm = TRUE),
                        Max = ~ max(., na.rm = TRUE),
                        Median = ~ median(., na.rm = TRUE),
                        Q3 = ~ quantile(., 0.75, na.rm = TRUE))))

median <- summary_stats_by_group %>%
  select(1:2, matches("2010_Median|2015_Median|2020_Median"))
write.csv(median,'F:/research/code_and_data/data_generation/figure2_change.csv',row.names = F)
# Max Min Median Figure S9------------------------------------------------------
rm(list = ls())
cat("\14")
library(reshape2)
library(openxlsx)
library(writexl)
library(plyr)
library(abind)
library(zoo)
library(readxl)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)
library(patchwork)


allIAMPath = "F:/research/code_and_data/data_generation/final_result_scenarios.xlsx"
allIAM = read_excel(path = allIAMPath, sheet = 'CO2_match')
allIAM$country = substr(allIAM$Region, nchar(allIAM$Region)-2, nchar(allIAM$Region))


grouped_data <- allIAM %>%
  group_by(country, Category)


a <- which(colnames(grouped_data) == '2010')
b <- which(colnames(grouped_data) == '2100')
selected_columns <- c(colnames(grouped_data)[a:b])


summary_stats_by_group <- grouped_data %>%
  summarize(across(all_of(selected_columns), 
                   list(Q1 = ~ quantile(., 0.25, na.rm = TRUE),
                        Min = ~ min(., na.rm = TRUE),
                        Max = ~ max(., na.rm = TRUE),
                        Median = ~ median(., na.rm = TRUE),
                        Q3 = ~ quantile(., 0.75, na.rm = TRUE))))

median <- summary_stats_by_group %>%
  select(1:2, contains("_median"))

median_long <- median %>%
  pivot_longer(cols = starts_with("2"),   
               names_to = "year",          
               values_to = "value")         
median_long <- median_long %>%
  separate(col = year, into = c("year", "median"), sep = "_", remove = TRUE) %>%
  select(-median) 
median_long <- median_long %>%
  pivot_wider(names_from = Category,  
              values_from = value)   
median_long <- median_long %>%
  filter(year %in% c(2030, 2050, 2080,2100))

Generate_Max <- summary_stats_by_group %>%
  select(1:2, contains("_Max"))

Generate_Max_long <- Generate_Max %>%
  pivot_longer(cols = starts_with("2"),  
               names_to = "year",          
               values_to = "value")        
Generate_Max_long <- Generate_Max_long %>%
  separate(col = year, into = c("year", "Max"), sep = "_", remove = TRUE) %>%
  select(-Max) 
Generate_Max_long <- Generate_Max_long %>%
  pivot_wider(names_from = Category,  
              values_from = value)        
Generate_Max_long <- Generate_Max_long %>%
  filter(year %in% c(2030, 2050, 2080,2100))

Generate_Min <- summary_stats_by_group %>%
  select(1:2, contains("_Min"))

Generate_Min_long <- Generate_Min %>%
  pivot_longer(cols = starts_with("2"),  
               names_to = "year",          
               values_to = "value")         
Generate_Min_long <- Generate_Min_long %>%
  separate(col = year, into = c("year", "Min"), sep = "_", remove = TRUE) %>%
  select(-Min) 
Generate_Min_long <- Generate_Min_long %>%
  pivot_wider(names_from = Category,  
              values_from = value)      
Generate_Min_long <- Generate_Min_long %>%
  filter(year %in% c(2030, 2050, 2080,2100))
write.csv(Generate_Min_long,'F:/research/code_and_data/data_generation/figureS9_min.csv',row.names = F)
write.csv(Generate_Max_long,'F:/research/code_and_data/data_generation/figureS9_max.csv',row.names = F)
write.csv(median_long,'F:/research/code_and_data/data_generation/figureS9_median.csv',row.names = F)



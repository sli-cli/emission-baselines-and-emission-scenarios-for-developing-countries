####################################################################################################
rm(list = ls())
cat("\14")

#====================================================================================================
wdir='F:/ISO3_v1.1' 
setwd(paste(wdir))
library(reshape2)
library(openxlsx)
library(writexl)
library(dplyr)
library(abind)
#====================================================================================================      
#====================================================================================================      

temp1=openxlsx::read.xlsx('./AR6_Scenarios_Database_metadata_indicators_v1.1.xlsx',sheet=2,colNames=TRUE)
ms = temp1[which(temp1$Category %in% c("C1",'C3','C5','C6','C7','C8')), c(1:4,8,12)]; dim(ms)
ms$ms = paste(ms$Model, ms$Scenario, sep='--') 
ms %>% dplyr::count(Category)
#====================================================================================================      
ms$Model_family = "Others"
ms[which(ms$Model %in% c("AIM/CGE 2.0","AIM/CGE 2.1","AIM/CGE 2.2", "AIM/Hub-Global 2.0")), ncol(ms)] = "AIM"
ms[which(ms$Model %in% c("C-ROADS-5.005")), ncol(ms)] = "C-ROADS"
ms[which(ms$Model %in% c("COFFEE 1.1")), ncol(ms)] = "COFFEE"
ms[which(ms$Model %in% c("EPPA 6")), ncol(ms)] = "EPPA"
ms[which(ms$Model %in% c("GCAM-PR 5.3", "GCAM 4.2","GCAM 5.3")), ncol(ms)] = "GCAM"
ms[which(ms$Model %in% c("GEM-E3_V2021")), ncol(ms)] = "GEM-E3"
ms[which(ms$Model %in% c("IMAGE 3.0","IMAGE 3.0.1","IMAGE 3.0.2","IMAGE 3.2")), ncol(ms)] = "IMAGE"
ms[which(ms$Model %in% c("MESSAGE-GLOBIOM 1.0","MESSAGEix-GLOBIOM 1.0","MESSAGEix-GLOBIOM_1.1",
                         "MESSAGEix-GLOBIOM_1.2","MESSAGEix-GLOBIOM_GEI 1.0")), ncol(ms)] = "MESSAGE"
ms[which(ms$Model %in% c("POLES ADVANCE","POLES EMF30","POLES EMF33","POLES ENGAGE","POLES GECO2019")), ncol(ms)] = "POLES"
ms[which(ms$Model %in% c("REMIND 1.6","REMIND 1.7","REMIND 2.1","REMIND-Buildings 2.0","REMIND-MAgPIE 1.5","REMIND-MAgPIE 1.7-3.0",
                         "REMIND-MAgPIE 2.0-4.1","REMIND-MAgPIE 2.1-4.2","REMIND-MAgPIE 2.1-4.3","REMIND-Transport 2.1")), ncol(ms)] = "REMIND"
ms[which(ms$Model %in% c("TIAM-ECN 1.1")), ncol(ms)] = "TIAM"
ms[which(ms$Model %in% c("WITCH 5.0","WITCH-GLOBIOM 3.1","WITCH-GLOBIOM 4.4")), ncol(ms)] = "WITCH"
#which(ms$Model_family=="Others")
ms <- ms[!(ms$Model_family %in% c('COFFEE', 'GEM-E3', 'TIAM', 'WITCH')), ]
length(unique(ms$Model_family)) # 12 model families

ms %>% dplyr::count(Category,Model_family)
ms %>% dplyr::count(Category, Project_study)
save(ms, file=paste(wdir,'./data generation/AR6-v1.1_C1-C8_metadata.Rdata',sep=''))

#====================================================================================================  
####################################################################################################
rm(list = ls())
cat("\14")

#====================================================================================================
wdir='F:/ISO3_v1.1' 
setwd(paste(wdir))
library(reshape2)
library(openxlsx)
library(writexl)
library(plyr)
library(abind)
library(zoo)
#====================================================================================================      
#====================================================================================================      
load('./data generation/AR6-v1.1_C1-C8_metadata.Rdata') 
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
save(data, file='./data generation/AR6-v1.1_C1-C3_allvars.Rdata')

#====================================================================================================      

load('./data generation/AR6-v1.1_C1-C3_allvars.Rdata') # data
sv=sort(c("Emissions|CO2|Energy"))

data=data[which(data$Variable %in% sv),]
years=c(2010:2100)
itv=seq(2010,2100,by=5)
ind=which(colnames(data) %in% itv)
ind2=which(colnames(data)=="Unit")
data=data[,c(1:ind2,ind)]


icol = which(colnames(data)=="2010")
temp = data[,icol:ncol(data)]
temp2 = matrix(NA, nrow=nrow(temp), ncol=ncol(temp))

for(i in 1:nrow(temp)){
  v = as.numeric(temp[i,])
  temp2[i,]=na.approx(v, na.rm=FALSE)
}

colnames(temp2)=seq(2010,2100,5)
count_na_func <- function(x) sum(is.na(x))
apply(data[,12:ncol(data)], 2, count_na_func) 
apply(temp2, 2, count_na_func) 

data = cbind(data[,1:(icol-1)], temp2)

pw = sort(unique(data$Category))
sv = sort(unique(data$Variable))
itv=seq(2010,2100,by=5)
icols = which(colnames(data) %in% itv)
ind2=which(colnames(data)=="Unit")
df = data[which(data$Variable %in% sv),c(1:ind2,icols)]
write.csv(df,'./data generation/ISO_ar6_mainmodel.csv',row.names = F)
#==================================================================================================== 
#==================================================================================================== 
rm(list = ls())
cat("\14")

#====================================================================================================
wdir='F:/R5_regions_v1.1' 
setwd(paste(wdir))
library(reshape2)
library(openxlsx)
library(writexl)
library(plyr)
library(abind)
library(zoo)
#====================================================================================================      
#====================================================================================================      
load('F:/ISO3_v1.1/data generation/AR6-v1.1_C1-C8_metadata.Rdata')
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
save(data, file='./data generation/AR6-v1.1_C1-C3_allvars.Rdata')

#====================================================================================================      
load('./data generation/AR6-v1.1_C1-C3_allvars.Rdata') 
sv=sort(c("Emissions|CO2|Energy"))

data=data[which(data$Variable %in% sv),]

years=c(2010:2100)
itv=seq(2010,2100,by=5)
ind=which(colnames(data) %in% itv)
ind2=which(colnames(data)=="Unit")
data=data[,c(1:ind2,ind)]

icol = which(colnames(data)=="2010")
temp = data[,icol:ncol(data)]
temp2 = matrix(NA, nrow=nrow(temp), ncol=ncol(temp))

for(i in 1:nrow(temp)){
  v = as.numeric(temp[i,])
  temp2[i,]=na.approx(v, na.rm=FALSE)
}

colnames(temp2)=seq(2010,2100,5)
count_na_func <- function(x) sum(is.na(x))
apply(data[,12:ncol(data)], 2, count_na_func) 
apply(temp2, 2, count_na_func) 

data = cbind(data[,1:(icol-1)], temp2)
#==================================================================================================== 
pw = sort(unique(data$Category))
sv = sort(unique(data$Variable))
itv=seq(2010,2100,by=5)
icols = which(colnames(data) %in% itv)
ind2=which(colnames(data)=="Unit")
df = data[which(data$Variable %in% sv),c(1:ind2,icols)]

write.csv(df,'./data generation/R5_ar6_mainmodel.csv',row.names = F)
# Max Min Median------------------------------------------------------
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

wdir='F:/ISO3_v1.1/data generation' 
allIAMPath = "F:/ISO3_v1.1/data generation/emissions_total.xlsx"
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
  select(country, year, AIM, GCAM, IMAGE, MESSAGE, POLES,REMIND)

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
  select(country, year, AIM, GCAM, IMAGE, MESSAGE, POLES,REMIND)

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
  select(country, year, AIM, GCAM, IMAGE, MESSAGE, POLES,REMIND)
write.csv(Generate_Min_long,'F:/ISO3_v1.1/data generation/figure1_min.csv',row.names = F)
write.csv(Generate_Max_long,'F:/ISO3_v1.1/data generation/figure1_max.csv',row.names = F)
write.csv(median_long,'F:/ISO3_v1.1/data generation/figure1_median.csv',row.names = F)

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

wdir='F:/ISO3_v1.1/data generation' 
allIAMPath = "F:/ISO3_v1.1/data generation/emissions_total.xlsx"
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
write.csv(median,'F:/ISO3_v1.1/data generation/figure2_change.csv',row.names = F)
# Max Min Median Figure S3------------------------------------------------------
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

wdir='F:/ISO3_v1.1/data generation' 
allIAMPath = "F:/ISO3_v1.1/data generation/final_result.xlsx"
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
write.csv(Generate_Min_long,'F:/ISO3_v1.1/data generation/figure3_min.csv',row.names = F)
write.csv(Generate_Max_long,'F:/ISO3_v1.1/data generation/figure3_max.csv',row.names = F)
write.csv(median_long,'F:/ISO3_v1.1/data generation/figure3_median.csv',row.names = F)



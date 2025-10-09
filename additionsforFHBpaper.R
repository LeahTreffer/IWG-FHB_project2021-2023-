# Find which Year had the lowest mean values for disease related traits 
# and which Year had highest values for disease related traits. 

library(data.table)

rm(list=ls(all=TRUE))
#set current working directory
#setwd()
setwd('/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS')

#load in phenotype tables
pheno <- read.delim('data/Intermediate_Files/Phenotypic_Format_long_2.txt') #read in file with all cycles and years and locations
#create separate tables for each cycle, year, location
#Cycle 
C7 <- pheno[str_detect(pheno$Cycle, "C7"),] #keep only data from C7 population
C10 <- pheno[str_detect(pheno$Cycle, "C10"),] #keep only data from C10 population
#Year
C7_2018 <- C7[str_detect(C7$phenotype_year, "2018"),]
C7_2019 <- C7[str_detect(C7$phenotype_year, "2019"),]
C7_2020 <- C7[str_detect(C7$phenotype_year, "2020"),]
C10_2021 <- C10[str_detect(C10$phenotype_year, "2021"),]
C10_2022 <- C10[str_detect(C10$phenotype_year, "2022"),]
#Location
C7_2018_S <- C7_2018[str_detect(C7_2018$Site, "SAL"),] #Salina
C7_2019_S <- C7_2019[str_detect(C7_2019$Site, "SAL"),] #Salina
C7_2020_S <- C7_2020[str_detect(C7_2020$Site, "SAL"),] #Salina
C10_2021_O <- C10_2021[str_detect(C10_2021$Site, "OLA"),] #Olathe
C10_2021_S <- C10_2021[str_detect(C10_2021$Site, "SAL"),] #Salina
C10_2022_O <- C10_2022[str_detect(C10_2022$Site, "OLA"),] #Olathe
C10_2022_S <- C10_2022[str_detect(C10_2022$Site, "SAL"),] #Salina



mean_list <- setNames(data.table(matrix(nrow = 0, ncol = 3)), c("trait", "mean", "year"))

#2018
eighteen <- split(C7_2018_S, C7_2018_S$trait_id) # Split data frame to get lists for each chr
unique(C7_2018_S$trait_id) # print the traits that are a part of this table
# 2018: "FHBDON"     "FHBD3G"     "FHBD3G.DON" "FHBDISIND"  "FHBSPKINC" "GLBLSEV"  "BLSSEV"     "HSATIVLSEV" 
traits <- c("FHBDON", "FHBD3G", "FHBD3G.DON", "FHBDISIND", "FHBSPKINC", "GLBLSEV", "BLSSEV", "HSATIVLSEV") #list of disease traits specific to 2018
for (i in traits){
  temp = eighteen[[i]] #table of one chromosome info
  m = mean(temp$phenotype_value, na.rm=TRUE)
  data=data.table(trait=c(i), mean=c(m), year="2018")
  mean_list <- rbindlist(list(mean_list, data))
}
rm(traits,i,temp,m,data)
#2019
nineteen <- split(C7_2019_S, C7_2019_S$trait_id) # Split data frame to get lists for each chr
unique(C7_2019_S$trait_id) # print the traits that are a part of this table
# 2019: "FHBDON"     "FHBD3G"     "FHBD3G.DON" "FHBDISIND"  "FHBSPKINC" "GLBLSEV"  "BLSSEV"     "HSATIVLSEV" 
traits <- c("FHBDON", "FHBD3G", "FHBD3G.DON", "FHBDISIND", "FHBSPKINC", "GLBLSEV", "BLSSEV", "HSATIVLSEV") #list of disease traits 
for (i in traits){
  temp = nineteen[[i]] #table of one chromosome info
  m = mean(temp$phenotype_value, na.rm=TRUE)
  data=data.table(trait=c(i), mean=c(m), year="2019")
  mean_list <- rbindlist(list(mean_list, data))
}
rm(traits,i,temp,m,data)
#2020
twenty <- split(C7_2020_S, C7_2020_S$trait_id) # Split data frame to get lists for each chr
unique(C7_2020_S$trait_id) # print the traits that are a part of this table
# 2019: "FHBDON"     "FHBD3G"     "FHBD3G.DON" "FHBDISIND"  "FHBSPKINC" "GLBLSEV"  "BLSSEV"     "HSATIVLSEV" 
traits <- c("FHBDON", "FHBD3G", "FHBD3G.DON", "FHBDISIND", "FHBSPKINC", "GLBLSEV", "BLSSEV", "HSATIVLSEV") #list of disease traits 
for (i in traits){
  temp = twenty[[i]] #table of one chromosome info
  m = mean(temp$phenotype_value, na.rm=TRUE)
  data=data.table(trait=c(i), mean=c(m), year="2020")
  mean_list <- rbindlist(list(mean_list, data))
}
rm(traits,i,temp,m,data)
#2021S
twentyoneS <- split(C10_2021_S, C10_2021_S$trait_id) # Split data frame to get lists for each chr
unique(C10_2021_S$trait_id) # print the traits that are a part of this table
# 2019: "FHBDON"     "FHBD3G"     "FHBD3G.DON" "FHBDISIND"  "FHBSPKINC" "GLBLSEV"  "BLSSEV"     "HSATIVLSEV" 
traits <- c("FHBDON", "FHBD3G", "FHBD3G.DON", "FHBDISIND", "FHBSPKINC", "GLBLSEV", "BLSSEV", "HSATIVLSEV") #list of disease traits 
for (i in traits){
  temp = twentyoneS[[i]] #table of one chromosome info
  m = mean(temp$phenotype_value, na.rm=TRUE)
  data=data.table(trait=c(i), mean=c(m), year="2021_SAL")
  mean_list <- rbindlist(list(mean_list, data))
}
rm(traits,i,temp,m,data)
#2021O
twentyoneO <- split(C10_2021_O, C10_2021_O$trait_id) # Split data frame to get lists for each chr
unique(C10_2021_O$trait_id) # print the traits that are a part of this table
# 2019: "FHBDON"     "FHBD3G"     "FHBD3G.DON" "FHBDISIND"  "FHBSPKINC" "GLBLSEV"  "BLSSEV"     "HSATIVLSEV" 
traits <- c("FHBDON", "FHBD3G", "FHBD3G.DON", "FHBDISIND", "FHBSPKINC", "GLBLSEV", "BLSSEV", "HSATIVLSEV") #list of disease traits 
for (i in traits){
  temp = twentyoneO[[i]] #table of one chromosome info
  m = mean(temp$phenotype_value, na.rm=TRUE)
  data=data.table(trait=c(i), mean=c(m), year="2021_OLA")
  mean_list <- rbindlist(list(mean_list, data))
}

rm(traits,i,temp,m,data)
#2022S
twentytwoS <- split(C10_2022_S, C10_2022_S$trait_id) # Split data frame to get lists for each chr
unique(C10_2022_S$trait_id) # print the traits that are a part of this table
# 2019: "FHBDON"     "FHBD3G"     "FHBD3G.DON" "FHBDISIND"  "FHBSPKINC" "GLBLSEV"  "BLSSEV"     "HSATIVLSEV" 
traits <- c("FHBDON", "FHBD3G", "FHBD3G.DON", "FHBDISIND", "FHBSPKINC", "GLBLSEV", "BLSSEV", "HSATIVLSEV", "FHBZEA") #list of disease traits 
for (i in traits){
  temp = twentytwoS[[i]] #table of one chromosome info
  m = mean(temp$phenotype_value, na.rm=TRUE)
  data=data.table(trait=c(i), mean=c(m), year="2022_SAL")
  mean_list <- rbindlist(list(mean_list, data))
}
rm(traits,i,temp,m,data)
#2022O
twentytwoO <- split(C10_2022_O, C10_2022_O$trait_id) # Split data frame to get lists for each chr
unique(C10_2022_O$trait_id) # print the traits that are a part of this table
# 2019: "FHBDON"     "FHBD3G"     "FHBD3G.DON" "FHBDISIND"  "FHBSPKINC" "GLBLSEV"  "BLSSEV"     "HSATIVLSEV" 
traits <- c("FHBDON", "FHBD3G", "FHBD3G.DON", "FHBDISIND", "FHBSPKINC", "GLBLSEV", "BLSSEV", "HSATIVLSEV", "FHBZEA") #list of disease traits 
for (i in traits){
  temp = twentytwoO[[i]] #table of one chromosome info
  m = mean(temp$phenotype_value, na.rm=TRUE)
  data=data.table(trait=c(i), mean=c(m), year="2022_OLA")
  mean_list <- rbindlist(list(mean_list, data))
}

mean_list <- mean_list[order(mean_list$trait, -mean_list$mean), ]

write.csv(mean_list, file="/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/DiseaseTraits_means.csv")



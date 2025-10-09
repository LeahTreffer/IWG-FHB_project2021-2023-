rm(list=ls(all=TRUE))


UA <- read.csv(file = '/Volumes/Backup Plus/FHB2022/unfinished_iwg.csv', header = TRUE, stringsAsFactors = FALSE) #unfinished
C10 <- read.csv(file = '~/Documents/2022_FHBC10_Salina_Fieldbookimport.csv', header = TRUE, stringsAsFactors = FALSE) #fieldbook upload from googledrive
C11 <- read.csv(file = '~/Documents/2022_FHBC11_Salina_Fieldbookimport.csv', header = TRUE, stringsAsFactors = FALSE) #fieldbook upload from googledrive

colnames(C11)[4] <- "PlantID" #source_id to 'plant_id'
#All <- merge(UA[,c(1)], C10[,c(3,4,6:8)], C11[,c(3,4,7:9)], by="CPGID", all=TRUE)
#pheno20 <- merge(pheno20[,c(2,4:5,7:14,16,18:20,23)], pheno[,c(1,2)], by="plant_id", all=TRUE) #merge files, using traits with data and removing empty traits

#C10.2 <- as.data.frame(match(UA, C10, nomatch = NA_integer_, incomparables = NULL))
ColC10= semi_join(C10, UA, by = "CPGID") #12
ColC11 = semi_join(C11[,-5], UA, by = "CPGID") #58

n_distinct(UA$CPGID) # there are 70 unique CPGID
#CPGID = rbind(C10[,c(1,3)],C11[,c(1,3)])
#noCollect = anti_join(UA, CPGID, by = "CPGID") # 2 of the CPGID in UA are not found in C10/C11 files, so the 68 CPGID total above adds up

write.csv(ColC10, file='/Volumes/Backup Plus/FHB2022/C10_LaterHarvestList.csv', row.names=FALSE, quote=FALSE)
write.csv(ColC11, file='/Volumes/Backup Plus/FHB2022/C11_LaterHarvestList.csv', row.names=FALSE, quote=FALSE)

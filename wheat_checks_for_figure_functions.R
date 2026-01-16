# wheat checks
# creates a list of tables for each fhb trait with a column for each year and column for the average phenotype rating  
phenow <- read.csv('data/Original_Files/plants.csv')
phenow$Site <- gsub("_Fusarium.*", "", phenow$experiment_id) #Remove characters after first allele
phenow$Site <- gsub(".*_","",phenow$Site) #Remove all before and up to "/"

overley <- phenow[str_detect(phenow$germplasm_id, "Overley"),] #
bacup <- phenow[str_detect(phenow$germplasm_id, "Bacup"),] #
roblin <- phenow[str_detect(phenow$germplasm_id, "Roblin"),] #
everest <- phenow[str_detect(phenow$germplasm_id, "Everest"),] #

checks <- data.frame()
checks <- rbind(overley,bacup,roblin,everest)

fhbdisind <- checks[str_detect(checks$trait_id, "FHBDISIND"),]
fhbspkinc <- checks[str_detect(checks$trait_id, "FHBSPKINC"),]
fhbdon <- checks[str_detect(checks$trait_id, "FHBDON"),]
fhbd3gdon <- checks[str_detect(checks$trait_id, "FHBD3G.DON"),]
#fhbzea <- checks[str_detect(checks$trait_id, "FHBZEA"),]

fhbdisind$phenotype_value <- as.numeric(fhbdisind$phenotype_value)
fhbspkinc$phenotype_value <- as.numeric(fhbspkinc$phenotype_value)
fhbdon$phenotype_value <- as.numeric(fhbdon$phenotype_value)
fhbd3gdon$phenotype_value <- as.numeric(fhbd3gdon$phenotype_value)
#fhbzea$phenotype_value <- as.numeric(fhbzea$phenotype_value)

fhbdisind$year_site <- paste(fhbdisind$phenotype_year,"_",fhbdisind$Site)
fhbspkinc$year_site <- paste(fhbspkinc$phenotype_year,"_",fhbspkinc$Site)
fhbdon$year_site <- paste(fhbdon$phenotype_year,"_",fhbdon$Site)
fhbd3gdon$year_site <- paste(fhbd3gdon$phenotype_year,"_",fhbd3gdon$Site)
#fhbzea$year_site <- paste(fhbzea$phenotype_year,"_",fhbzea$Site)

# 2018 FHBSPKINC in wrong format
fhbspkinc <- fhbspkinc %>%
  mutate(phenotype_value = if_else(phenotype_year == "2018", phenotype_value * 100, phenotype_value)) # multiply 2018 values by 100

STATS_fhbdisind<- fhbdisind %>% group_by(year_site) %>% 
  dplyr::summarise(trait_mean = mean(na.omit(phenotype_value)))
STATS_fhbspkinc<- fhbspkinc %>% group_by(year_site) %>% 
  dplyr::summarise(trait_mean = mean(na.omit(phenotype_value)))
STATS_fhbdon<- fhbdon %>% group_by(year_site) %>% 
  dplyr::summarise(trait_mean = mean(na.omit(phenotype_value)))
STATS_fhbd3gdon<- fhbd3gdon %>% group_by(year_site) %>% 
  dplyr::summarise(trait_mean = mean(na.omit(phenotype_value)))
#STATS_fhbzea<- fhbzea %>% group_by(year_site) %>% 
  #dplyr::summarise(trait_mean = mean(na.omit(phenotype_value)))


STATS_fhbdisind$Plot <- c("C07_SAL_2018","C07_SAL_2019","C07_OLA_2020","C07_SAL_2020","C10_OLA_2021","C10_SAL_2021","C10_OLA_2022","C10_SAL_2022")
STATS_fhbspkinc$Plot <- c("C07_SAL_2018","C07_SAL_2019","C07_OLA_2020","C07_SAL_2020","C10_OLA_2021","C10_SAL_2021","C10_OLA_2022","C10_SAL_2022")
STATS_fhbdon$Plot <- c("C07_SAL_2018","C07_SAL_2019","C10_OLA_2021","C10_SAL_2021")
STATS_fhbd3gdon$Plot <- c("C07_SAL_2018","C07_SAL_2019","C10_OLA_2021","C10_SAL_2021")

STATS_fhbdisind <- STATS_fhbdisind[c(1:2,4:8),c(3,2)]
STATS_fhbspkinc <- STATS_fhbspkinc[c(1:2,4:8),c(3,2)]
STATS_fhbdon <- STATS_fhbdon[,c(3,2)]
STATS_fhbd3gdon <- STATS_fhbd3gdon[,c(3,2)]

List2 <- list(FHBDISIND = STATS_fhbdisind, FHBSPKINC = STATS_fhbspkinc, FHBDON = STATS_fhbdon, FHBD3GDON = STATS_fhbd3gdon)

# phenotype rating changed from numeric continuous to discrete
# I only care about doing this for FHBSPKINC and FHBDISIND for the bar plot figure so only writing code for those
STATS_fhbspkinc$select <- STATS_fhbspkinc$trait_mean
STATS_fhbspkinc$select <- ifelse(STATS_fhbspkinc$trait_mean <= 0, "0", STATS_fhbspkinc$select) 
STATS_fhbspkinc$select <- ifelse(STATS_fhbspkinc$trait_mean > 0 & STATS_fhbspkinc$trait_mean <= 10, "1-10", STATS_fhbspkinc$select) 
STATS_fhbspkinc$select <- ifelse(STATS_fhbspkinc$trait_mean > 10 & STATS_fhbspkinc$trait_mean <= 20, "11-20", STATS_fhbspkinc$select) 
STATS_fhbspkinc$select <- ifelse(STATS_fhbspkinc$trait_mean > 20 & STATS_fhbspkinc$trait_mean <= 30, "21-30", STATS_fhbspkinc$select) 
STATS_fhbspkinc$select <- ifelse(STATS_fhbspkinc$trait_mean > 30 & STATS_fhbspkinc$trait_mean <= 40, "31-40", STATS_fhbspkinc$select) 
STATS_fhbspkinc$select <- ifelse(STATS_fhbspkinc$trait_mean > 40 & STATS_fhbspkinc$trait_mean <= 50, "41-50", STATS_fhbspkinc$select) 
STATS_fhbspkinc$select <- ifelse(STATS_fhbspkinc$trait_mean > 50 & STATS_fhbspkinc$trait_mean <= 60, "51-60", STATS_fhbspkinc$select) 
STATS_fhbspkinc$select <- ifelse(STATS_fhbspkinc$trait_mean > 60 & STATS_fhbspkinc$trait_mean <= 70, "61-70", STATS_fhbspkinc$select) 
STATS_fhbspkinc$select <- ifelse(STATS_fhbspkinc$trait_mean > 70 & STATS_fhbspkinc$trait_mean <= 80, "71-80", STATS_fhbspkinc$select) 
STATS_fhbspkinc$select <- ifelse(STATS_fhbspkinc$trait_mean > 80 & STATS_fhbspkinc$trait_mean <= 90, "81-90", STATS_fhbspkinc$select) 
STATS_fhbspkinc$select <- ifelse(STATS_fhbspkinc$trait_mean > 90 & STATS_fhbspkinc$trait_mean <= 100, "91-100", STATS_fhbspkinc$select) 

STATS_fhbdisind$select <- STATS_fhbdisind$trait_mean
STATS_fhbdisind$select <- ifelse(STATS_fhbdisind$trait_mean <= 0, "0", STATS_fhbdisind$select) 
STATS_fhbdisind$select <- ifelse(STATS_fhbdisind$trait_mean > 0 & STATS_fhbdisind$trait_mean <= 0.5, "0.01-0.5", STATS_fhbdisind$select) 
STATS_fhbdisind$select <- ifelse(STATS_fhbdisind$trait_mean > 0.5 & STATS_fhbdisind$trait_mean <= 1, "0.51-1", STATS_fhbdisind$select) 
STATS_fhbdisind$select <- ifelse(STATS_fhbdisind$trait_mean > 1 & STATS_fhbdisind$trait_mean <= 2, "1.01-2", STATS_fhbdisind$select) 
STATS_fhbdisind$select <- ifelse(STATS_fhbdisind$trait_mean > 2 & STATS_fhbdisind$trait_mean <= 3, "2.01-3", STATS_fhbdisind$select) 
STATS_fhbdisind$select <- ifelse(STATS_fhbdisind$trait_mean > 3 & STATS_fhbdisind$trait_mean <= 4, "3.01-4", STATS_fhbdisind$select) 
STATS_fhbdisind$select <- ifelse(STATS_fhbdisind$trait_mean > 4 & STATS_fhbdisind$trait_mean <= 5, "4.01-5", STATS_fhbdisind$select) 

List3 <- list(FHBDISIND = STATS_fhbdisind, FHBSPKINC = STATS_fhbspkinc, FHBDON = STATS_fhbdon)

save(List2, file = "data/Intermediate_Files/wheat_means.RData")



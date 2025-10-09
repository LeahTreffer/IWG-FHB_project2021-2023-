rm(list=ls(all=TRUE))
#set current working directory
#setwd()
setwd('/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS')

library(tidyverse)
library(ggpubr)

data <- read.delim("data/Intermediate_Files/Phenotypic_Format_wide.txt")

C07 <- data[str_detect(data$Cycle, "C07"),] #keep only data from C7 population
C10 <- data[str_detect(data$Cycle, "C10"),] #keep only data from C7 population

#log transformations
C07$logFHBSPKINC <- log(C07$FHBSPKINC)
C07$logFHBSPKINC[is.na(C07$logFHBSPKINC) | C07$logFHBSPKINC=="-Inf"] = NA
C07$logFHBDISIND <- log(C07$FHBDISIND)
C07$logFHBDISIND[is.na(C07$logFHBDISIND) | C07$logFHBDISIND=="-Inf"] = NA
C07$logGLBLSEV <- log(C07$GLBLSEV)
C07$logGLBLSEV[is.na(C07$logGLBLSEV) | C07$logGLBLSEV=="-Inf"] = NA
C07$logFHBDON <- log(C07$FHBDON)
C07$logFHBDON[is.na(C07$logFHBDON) | C07$logFHBDON=="-Inf"] = NA

C07$logSPKYLD <- log(C07$SPKYLD)
C07$logSPKYLD[is.na(C07$logSPKYLD) | C07$logSPKYLD=="-Inf"] = NA
C07$logPTHT <- log(C07$PTHT)
C07$logPTHT[is.na(C07$logPTHT) | C07$logPTHT=="-Inf"] = NA
C07$logZDK <- log(C07$ZDK)
C07$logZDK[is.na(C07$logZDK) | C07$logZDK=="-Inf"] = NA

hist(C07$logSPKYLD)
hist(C07$SPKYLD)
hist(C07$logPTHT)
hist(C07$PTHT)
hist(C07$logZDK)
hist(C07$ZDK)
hist(C07$logFHBDISIND)
hist(C07$FHBDISIND)
hist(C07$logFHBSPKINC)
hist(C07$FHBSPKINC)

# sqrt transform
C07$sqrtFHBSPKINC <- sqrt(C07$FHBSPKINC)
C07$sqrtFHBSPKINC[is.na(C07$sqrtFHBSPKINC) | C07$sqrtFHBSPKINC=="-Inf"] = NA
C07$sqrtFHBDISIND <- sqrt(C07$FHBDISIND)
C07$sqrtFHBDISIND[is.na(C07$sqrtFHBDISIND) | C07$sqrtFHBDISIND=="-Inf"] = NA
C07$sqrtGLBLSEV <- sqrt(C07$GLBLSEV)
C07$sqrtGLBLSEV[is.na(C07$sqrtGLBLSEV) | C07$sqrtGLBLSEV=="-Inf"] = NA
C07$sqrtFHBDON <- sqrt(C07$FHBDON)
C07$sqrtFHBDON[is.na(C07$sqrtFHBDON) | C07$sqrtFHBDON=="-Inf"] = NA

C07$sqrtSPKYLD <- sqrt(C07$SPKYLD)
C07$sqrtSPKYLD[is.na(C07$sqrtSPKYLD) | C07$sqrtSPKYLD=="-Inf"] = NA
C07$sqrtPTHT <- sqrt(C07$PTHT)
C07$sqrtPTHT[is.na(C07$sqrtPTHT) | C07$sqrtPTHT=="-Inf"] = NA
C07$sqrtZDK <- sqrt(C07$ZDK)
C07$sqrtZDK[is.na(C07$sqrtZDK) | C07$sqrtZDK=="-Inf"] = NA

hist(C07$sqrtSPKYLD, breaks = 100) #
hist(C07$logSPKYLD, breaks = 100)
hist(C07$SPKYLD, breaks = 100) #
hist(C07$sqrtPTHT, breaks = 100) #
hist(C07$logPTHT, breaks = 100) #
hist(C07$PTHT, breaks = 100) #
hist(C07$sqrtZDK, breaks = 100) #
hist(C07$logZDK, breaks = 100) #
hist(C07$ZDK, breaks = 100) #
hist(C07$sqrtFHBDISIND, breaks = 100) #
hist(C07$logFHBDISIND, breaks = 100)
hist(C07$FHBDISIND, breaks = 100)
hist(C07$sqrtFHBSPKINC, breaks = 100) 
hist(C07$logFHBSPKINC, breaks = 100)
hist(C07$FHBSPKINC, breaks = 100) #


C07_2018_SAL <- C07[str_detect(C07$phenotype_year, "2018"),]
C07_2019_SAL <- C07[str_detect(C07$phenotype_year, "2019"),]
C07_2020_SAL <- C07[str_detect(C07$phenotype_year, "2020"),]

C07$phenotype_year <- as.factor(C07$phenotype_year)

a <- ggscatter(
  C07, x = "FHBSPKINC", y = "SPKYLD", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line",
  xlab = FALSE
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 0.25) +
  stat_regline_equation(label.y = 0.23)

b <- ggscatter(
  C07, x = "FHBSPKINC", y = "PTHT", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line",
  xlab = FALSE
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 180) +
  stat_regline_equation(label.y = 175)

c <- ggscatter(
  C07, x = "FHBSPKINC", y = "ZDK", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 85) +
  stat_regline_equation(label.y = 84)

figure <- ggarrange(a, b, c,
                    #labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3, 
                    common.legend = TRUE)
figure

#Export to figures file as regressions_for_Paper.pdf
## Dimensions : 15.69 x 14.68

d <- ggscatter(
  C07, x = "FHBDISIND", y = "SPKYLD", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line",
  xlab = FALSE
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 0.26) +
  stat_regline_equation(label.y = 0.25)

e <- ggscatter(
  C07, x = "FHBDISIND", y = "PTHT", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line",
  xlab = FALSE
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 182) +
  stat_regline_equation(label.y = 178)

f <- ggscatter(
  C07, x = "FHBDISIND", y = "ZDK", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 86) +
  stat_regline_equation(label.y = 84)

figure <- ggarrange(d, e, f,
                    #labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3, 
                    common.legend = TRUE)
figure

#Export to figures file as regressions_for_Paper2.pdf
## Dimensions : 15.69 x 14.68






# try with log values

g <- ggscatter(
  C07, x = "logFHBDISIND", y = "logSPKYLD", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line",
  xlab = FALSE
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 0) +
  stat_regline_equation(label.y = -1)

h <- ggscatter(
  C07, x = "logFHBDISIND", y = "logPTHT", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line",
  xlab = FALSE
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 5.75) +
  stat_regline_equation(label.y = 5.5)

i <- ggscatter(
  C07, x = "logFHBDISIND", y = "logZDK", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 4.6) +
  stat_regline_equation(label.y = 4.53)

figure <- ggarrange(g, h, i,
                    #labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3, 
                    common.legend = TRUE)
figure


# try with sqrt values

j <- ggscatter(
  C07, x = "sqrtFHBDISIND", y = "sqrtSPKYLD", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 0.6) +
  stat_regline_equation(label.y = 0.55)

k <- ggscatter(
  C07, x = "sqrtFHBDISIND", y = "sqrtPTHT", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 16) +
  stat_regline_equation(label.y = 15)

l <- ggscatter(
  C07, x = "sqrtSPKYLD", y = "sqrtPTHT", 
  color = "phenotype_year", palette = "jco",
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 15) +
  stat_regline_equation(label.y = 14)

figure <- ggarrange(j, k, l,
                    #labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3, 
                    common.legend = TRUE)
figure

#Export to figures file as regressions_for_Paper3.pdf
## Dimensions : 15.69 x 14.68



####
#log transformations
C07$logFHBSPKINC <- log(C07$FHBSPKINC)
C07$logFHBSPKINC[is.na(C07$logFHBSPKINC) | C07$logFHBSPKINC=="-Inf"] = NA
C07$logFHBDISIND <- log(C07$FHBDISIND)
C07$logFHBDISIND[is.na(C07$logFHBDISIND) | C07$logFHBDISIND=="-Inf"] = NA
C07$logGLBLSEV <- log(C07$GLBLSEV)
C07$logGLBLSEV[is.na(C07$logGLBLSEV) | C07$logGLBLSEV=="-Inf"] = NA
C07$logFHBDON <- log(C07$FHBDON)
C07$logFHBDON[is.na(C07$logFHBDON) | C07$logFHBDON=="-Inf"] = NA

C07$logSPKYLD <- log(C07$SPKYLD)
C07$logSPKYLD[is.na(C07$logSPKYLD) | C07$logSPKYLD=="-Inf"] = NA
C07$logPTHT <- log(C07$PTHT)
C07$logPTHT[is.na(C07$logPTHT) | C07$logPTHT=="-Inf"] = NA
C07$logZDK <- log(C07$ZDK)
C07$logZDK[is.na(C07$logZDK) | C07$logZDK=="-Inf"] = NA

hist(C07$logSPKYLD)
hist(C07$SPKYLD)
hist(C07$logPTHT)
hist(C07$PTHT)
hist(C07$logZDK)
hist(C07$ZDK)
hist(C07$logFHBDISIND)
hist(C07$FHBDISIND)
hist(C07$logFHBSPKINC)
hist(C07$FHBSPKINC)

# sqrt transform
C07$sqrtFHBSPKINC <- sqrt(C07$FHBSPKINC)
C07$sqrtFHBSPKINC[is.na(C07$sqrtFHBSPKINC) | C07$sqrtFHBSPKINC=="-Inf"] = NA
C07$sqrtFHBDISIND <- sqrt(C07$FHBDISIND)
C07$sqrtFHBDISIND[is.na(C07$sqrtFHBDISIND) | C07$sqrtFHBDISIND=="-Inf"] = NA
C07$sqrtGLBLSEV <- sqrt(C07$GLBLSEV)
C07$sqrtGLBLSEV[is.na(C07$sqrtGLBLSEV) | C07$sqrtGLBLSEV=="-Inf"] = NA
C07$sqrtFHBDON <- sqrt(C07$FHBDON)
C07$sqrtFHBDON[is.na(C07$sqrtFHBDON) | C07$sqrtFHBDON=="-Inf"] = NA

C07$sqrtSPKYLD <- sqrt(C07$SPKYLD)
C07$sqrtSPKYLD[is.na(C07$sqrtSPKYLD) | C07$sqrtSPKYLD=="-Inf"] = NA
C07$sqrtPTHT <- sqrt(C07$PTHT)
C07$sqrtPTHT[is.na(C07$sqrtPTHT) | C07$sqrtPTHT=="-Inf"] = NA
C07$sqrtZDK <- sqrt(C07$ZDK)
C07$sqrtZDK[is.na(C07$sqrtZDK) | C07$sqrtZDK=="-Inf"] = NA

hist(C07$sqrtSPKYLD, breaks = 100) #
hist(C07$logSPKYLD, breaks = 100)
hist(C07$SPKYLD, breaks = 100) #
hist(C07$sqrtPTHT, breaks = 100) #
hist(C07$logPTHT, breaks = 100) #
hist(C07$PTHT, breaks = 100) #
hist(C07$sqrtZDK, breaks = 100) #
hist(C07$logZDK, breaks = 100) #
hist(C07$ZDK, breaks = 100) #
hist(C07$sqrtFHBDISIND, breaks = 100) #
hist(C07$logFHBDISIND, breaks = 100)
hist(C07$FHBDISIND, breaks = 100)
hist(C07$sqrtFHBSPKINC, breaks = 100) 
hist(C07$logFHBSPKINC, breaks = 100)
hist(C07$FHBSPKINC, breaks = 100) #


C07_2018_SAL <- C07[str_detect(C07$phenotype_year, "2018"),]
C07_2019_SAL <- C07[str_detect(C07$phenotype_year, "2019"),]
C07_2020_SAL <- C07[str_detect(C07$phenotype_year, "2020"),]

C07$phenotype_year <- as.factor(C07$phenotype_year)

# C10 
# C10 only has traits: FHBDON, FHBDISIND, FHBSPKINC, PTHT, ZDK, BLSSEV, HSATIVLSEV, 
#log transformations
C10$logFHBSPKINC <- log(C10$FHBSPKINC)
C10$logFHBSPKINC[is.na(C10$logFHBSPKINC) | C10$logFHBSPKINC=="-Inf"] = NA
C10$logFHBDISIND <- log(C10$FHBDISIND)
C10$logFHBDISIND[is.na(C10$logFHBDISIND) | C10$logFHBDISIND=="-Inf"] = NA
C10$logFHBDON <- log(C10$FHBDON)
C10$logFHBDON[is.na(C10$logFHBDON) | C10$logFHBDON=="-Inf"] = NA
C10$logPTHT <- log(C10$PTHT)
C10$logPTHT[is.na(C10$logPTHT) | C10$logPTHT=="-Inf"] = NA
C10$logZDK <- log(C10$ZDK)
C10$logZDK[is.na(C10$logZDK) | C10$logZDK=="-Inf"] = NA

hist(C10$logPTHT)
hist(C10$PTHT)
hist(C10$logZDK)
hist(C10$ZDK)
hist(C10$logFHBDISIND)
hist(C10$FHBDISIND)
hist(C10$logFHBSPKINC)
hist(C10$FHBSPKINC)

# sqrt transform
C10$sqrtFHBSPKINC <- sqrt(C10$FHBSPKINC)
C10$sqrtFHBSPKINC[is.na(C10$sqrtFHBSPKINC) | C10$sqrtFHBSPKINC=="-Inf"] = NA
C10$sqrtFHBDISIND <- sqrt(C10$FHBDISIND)
C10$sqrtFHBDISIND[is.na(C10$sqrtFHBDISIND) | C10$sqrtFHBDISIND=="-Inf"] = NA
C10$sqrtFHBDON <- sqrt(C10$FHBDON)
C10$sqrtFHBDON[is.na(C10$sqrtFHBDON) | C10$sqrtFHBDON=="-Inf"] = NA
C10$sqrtPTHT <- sqrt(C10$PTHT)
C10$sqrtPTHT[is.na(C10$sqrtPTHT) | C10$sqrtPTHT=="-Inf"] = NA
C10$sqrtZDK <- sqrt(C10$ZDK)
C10$sqrtZDK[is.na(C10$sqrtZDK) | C10$sqrtZDK=="-Inf"] = NA

hist(C10$sqrtPTHT, breaks = 100) #
hist(C10$logPTHT, breaks = 100) #
hist(C10$PTHT, breaks = 100) #
hist(C10$sqrtZDK, breaks = 100) #
hist(C10$logZDK, breaks = 100) #
hist(C10$ZDK, breaks = 100) #
hist(C10$sqrtFHBDISIND, breaks = 100) #
hist(C10$logFHBDISIND, breaks = 100)
hist(C10$FHBDISIND, breaks = 100)
hist(C10$sqrtFHBSPKINC, breaks = 100) 
hist(C10$logFHBSPKINC, breaks = 100)
hist(C10$FHBSPKINC, breaks = 100) #


C10_2021_SAL <- C10[str_detect(C10$phenotype_year, "2021"),] 
C10_2021_SAL <- C10_2021_SAL[str_detect(C10_2021_SAL$Site, "SAL"),]
C10_2022_SAL <- C10[str_detect(C10$phenotype_year, "2022"),]
C10_2022_SAL <- C10_2022_SAL[str_detect(C10_2022_SAL$Site, "SAL"),]
C10_2021_OLA <- C10[str_detect(C10$phenotype_year, "2021"),]
C10_2021_OLA <- C10_2021_OLA[str_detect(C10_2021_OLA$Site, "OLA"),]
C10_2022_OLA <- C10[str_detect(C10$phenotype_year, "2022"),]
C10_2022_OLA <- C10_2022_OLA[str_detect(C10_2022_OLA$Site, "OLA"),]

C10$phenotype_year <- as.factor(C10$phenotype_year)

# Adjust the margins: c(bottom, left, top, right)
par(mar=c(4, 4, 2, 2))
# Set up the plotting area to be 3 plots across and 4 plots down
par(mfrow=c(4, 3))
hist(C10$sqrtPTHT, breaks = 100) #
hist(C10$logPTHT, breaks = 100) #
hist(C10$PTHT, breaks = 100) #
hist(C10$sqrtZDK, breaks = 100) #
hist(C10$logZDK, breaks = 100) #
hist(C10$ZDK, breaks = 100) #
hist(C10$sqrtFHBDISIND, breaks = 100) #
hist(C10$logFHBDISIND, breaks = 100)
hist(C10$FHBDISIND, breaks = 100)
hist(C10$sqrtFHBSPKINC, breaks = 100) 
hist(C10$logFHBSPKINC, breaks = 100)
hist(C10$FHBSPKINC, breaks = 100) #
par(mfrow=c(1, 1))


C10$phenotype_year_site <- paste(C10$phenotype_year,C10$Site)




j <- ggscatter(
  C07, x = "sqrtFHBDISIND", y = "sqrtSPKYLD", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 0.6) +
  stat_regline_equation(label.y = 0.55)

k <- ggscatter(
  C07, x = "sqrtFHBDISIND", y = "sqrtPTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 16) +
  stat_regline_equation(label.y = 15) +
  theme(legend.position = "none")  # This line removes the legend

l <- ggscatter(
  C07, x = "sqrtSPKYLD", y = "sqrtPTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 15) +
  stat_regline_equation(label.y = 14)+
  theme(legend.position = "none")  # This line removes the legend


# no spkyld in C10

n <- ggscatter(
  C10, x = "sqrtFHBDISIND", y = "sqrtPTHT", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 16) +
  stat_regline_equation(label.y = 15)

# no spkyld in c10




figure <- ggarrange(j, NULL, k, n, l, NULL,
                    #labels = c("A", "B", "C"),
                    ncol = 2, nrow = 3, 
                    common.legend = FALSE)
figure

#Export to figures file as regressions_for_Paper3.pdf
## Dimensions : 15.69 x 14.68


# try regression by trait and include all possible years
#FHBDON
#FHBD3G
#FHBD3G.DON
#FHBDISIND
#FHBSPKINC

#PTHT
#STMANG
#HDANG
#SPKLNG
#SPKYLD
#SPKDEN




j <- ggscatter(
  C07, x = "FHBDON", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)

k <- ggscatter(
  C07, x = "FHBD3G", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185) +
  theme(legend.position = "none")  # This line removes the legend

l <- ggscatter(
  C07, x = "FHBD3G.DON", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

m <- ggscatter(
  C07, x = "FHBDISIND", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

n <- ggscatter(
  C07, x = "FHBSPKINC", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

o <- ggscatter(
  C10, x = "FHBDON", y = "PTHT", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  labs(x = NULL, y = NULL)  # This line removes the axis labels

p <- ggscatter(
  C10, x = "FHBD3G", y = "PTHT", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none")  # This line removes the legend

q <- ggscatter(
  C10, x = "FHBD3G.DON", y = "PTHT", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none")  # This line removes the legend

r <- ggscatter(
  C10, x = "FHBDISIND", y = "PTHT", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none")  # This line removes the legend

s <- ggscatter(
  C10, x = "FHBSPKINC", y = "PTHT", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none")  # This line removes the legend

t <- ggscatter(
  C07, x = "SPKYLD", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")+  # This line removes the legend
  theme(legend.position = "none")  # This line removes the legend


figure <- ggarrange(j, o, k, p, l, q, m, r, n, s, t, NULL,
                    #labels = c("A", "B", "C"),
                    ncol = 2, nrow = 6, 
                    common.legend = FALSE)
figure




j <- ggscatter(
  C07, x = "FHBDON", y = "STMANG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)

k <- ggscatter(
  C07, x = "FHBD3G", y = "STMANG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185) +
  theme(legend.position = "none")  # This line removes the legend

l <- ggscatter(
  C07, x = "FHBD3G.DON", y = "STMANG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

m <- ggscatter(
  C07, x = "FHBDISIND", y = "STMANG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

n <- ggscatter(
  C07, x = "FHBSPKINC", y = "STMANG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

figure <- ggarrange(j, k, l, m, n,
                    #labels = c("A", "B", "C"),
                    ncol = 1, nrow = 5, 
                    common.legend = FALSE)
figure





j <- ggscatter(
  C07, x = "FHBDON", y = "HDANG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)

k <- ggscatter(
  C07, x = "FHBD3G", y = "HDANG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185) +
  theme(legend.position = "none")  # This line removes the legend

l <- ggscatter(
  C07, x = "FHBD3G.DON", y = "HDANG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

m <- ggscatter(
  C07, x = "FHBDISIND", y = "HDANG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

n <- ggscatter(
  C07, x = "FHBSPKINC", y = "HDANG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

figure <- ggarrange(j, k, l, m, n,
                    #labels = c("A", "B", "C"),
                    ncol = 1, nrow = 5, 
                    common.legend = FALSE)
figure






j <- ggscatter(
  C07, x = "FHBDON", y = "SPKLNG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)

k <- ggscatter(
  C07, x = "FHBD3G", y = "SPKLNG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185) +
  theme(legend.position = "none")  # This line removes the legend

l <- ggscatter(
  C07, x = "FHBD3G.DON", y = "SPKLNG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

m <- ggscatter(
  C07, x = "FHBDISIND", y = "SPKLNG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

n <- ggscatter(
  C07, x = "FHBSPKINC", y = "SPKLNG", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

figure <- ggarrange(j, k, l, m, n,
                    #labels = c("A", "B", "C"),
                    ncol = 1, nrow = 5, 
                    common.legend = FALSE)
figure





j <- ggscatter(
  C07, x = "FHBDON", y = "SPKYLD", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)

k <- ggscatter(
  C07, x = "FHBD3G", y = "SPKYLD", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185) +
  theme(legend.position = "none")  # This line removes the legend

l <- ggscatter(
  C07, x = "FHBD3G.DON", y = "SPKYLD", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

m <- ggscatter(
  C07, x = "FHBDISIND", y = "SPKYLD", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

n <- ggscatter(
  C07, x = "FHBSPKINC", y = "SPKYLD", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

figure <- ggarrange(j, k, l, m, n,
                    #labels = c("A", "B", "C"),
                    ncol = 1, nrow = 5, 
                    common.legend = FALSE)
figure





j <- ggscatter(
  C07, x = "FHBDON", y = "SPKDEN", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)

k <- ggscatter(
  C07, x = "FHBD3G", y = "SPKDEN", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185) +
  theme(legend.position = "none")  # This line removes the legend

l <- ggscatter(
  C07, x = "FHBD3G.DON", y = "SPKDEN", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

m <- ggscatter(
  C07, x = "FHBDISIND", y = "SPKDEN", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

n <- ggscatter(
  C07, x = "FHBSPKINC", y = "SPKDEN", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none")  # This line removes the legend

figure <- ggarrange(j, k, l, m, n,
                    #labels = c("A", "B", "C"),
                    ncol = 1, nrow = 5, 
                    common.legend = FALSE)
figure







j <- ggscatter(
  C07, x = "FHBDISIND", y = "FHBDON", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 30) +
  stat_regline_equation(label.y = 25)

k <- ggscatter(
  C07, x = "FHBSPKINC", y = "FHBDON", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 30) +
  stat_regline_equation(label.y = 25) +
  theme(legend.position = "none")  # This line removes the legend

l <- ggscatter(
  C07, x = "FHBDISIND", y = "FHBSPKINC", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 150) +
  stat_regline_equation(label.y = 145)+
  theme(legend.position = "none")  # This line removes the legend

m <- ggscatter(
  C07, x = "FHBD3G", y = "FHBDON", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 30) +
  stat_regline_equation(label.y = 25)+
  theme(legend.position = "none")  # This line removes the legend

n <- ggscatter(
  C07, x = "FHBD3G.DON", y = "FHBDON", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 30) +
  stat_regline_equation(label.y = 25)+
  theme(legend.position = "none")  # This line removes the legend

o <- ggscatter(
  C10, x = "FHBDISIND", y = "FHBDON", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 150) +
  stat_regline_equation(label.y = 145)+
  labs(x = NULL, y = NULL)  # This line removes the axis labels

p <- ggscatter(
  C10, x = "FHBSPKINC", y = "FHBDON", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 150) +
  stat_regline_equation(label.y = 145)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none")  # This line removes the legend

q <- ggscatter(
  C10, x = "FHBDISIND", y = "FHBSPKINC", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 150) +
  stat_regline_equation(label.y = 145)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none")  # This line removes the legend

r <- ggscatter(
  C10, x = "FHBD3G", y = "FHBDON", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none")  # This line removes the legend

s <- ggscatter(
  C10, x = "FHBD3G.DON", y = "FHBDON", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 190) +
  stat_regline_equation(label.y = 185)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none")  # This line removes the legend


figure <- ggarrange(j, o, k, p, l, q, m, r, n, s,
                    #labels = c("A", "B", "C"),
                    ncol = 2, nrow = 5, 
                    common.legend = FALSE)
figure

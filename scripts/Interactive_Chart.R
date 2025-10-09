# Interactive Chart

# Broman KW (2015) R/qtlcharts: interactive graphics for quantitative trait locus mapping. Genetics 199:359-361.

knitr::opts_chunk$set(fig.width=8, fig.height=6, message=FALSE)

install.packages("qtl")
library(qtl)
install.packages("qtlcharts")
library(qtlcharts)

Sig <- read.csv('data/Intermediate_Files/SigM.csv', head=TRUE)

Sig <- Sig[,c("Chr","Pos","ID","SNP")]
colnames(Sig)[1] <- "group"
colnames(Sig)[2] <- "postion"
colnames(Sig)[3] <- "id"
colnames(Sig)[4] <- "locus"

j01 <- Sig[str_detect(Sig$group, "J01"),]
J01 <- split(j01$postion, j01$locus)
j02 <- Sig[str_detect(Sig$group, "J02"),]
J02 <- split(j02$postion, j02$locus)
j03 <- Sig[str_detect(Sig$group, "J03"),]
J03 <- split(j03$postion, j03$locus)
j04 <- Sig[str_detect(Sig$group, "J04"),]
J04 <- split(j04$postion, j04$locus)
j05 <- Sig[str_detect(Sig$group, "J05"),]
J05 <- split(j05$postion, j05$locus)
j06 <- Sig[str_detect(Sig$group, "J06"),]
J06 <- split(j06$postion, j06$locus)
j07 <- Sig[str_detect(Sig$group, "J07"),]
J07 <- split(j07$postion, j07$locus)


map3 <- list(J01 = J01,J02 = J02, J03 = J03, J04 = J04, J05 = J05, J06 = J06, J07 = J07)

iplotMap(map3)

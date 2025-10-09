rm(list=ls(all=TRUE))
#set current working directory
#setwd()
setwd('/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS')

###
FHB_for_GS_GWAS <- read.table(file = 'data/Original_Files/FHB_for_GS_GWAS.txt', header = TRUE, row.names = 1, stringsAsFactors = FALSE) #load genotype calls
#make a copy of master file to work with
GENO <- FHB_for_GS_GWAS
geno <- GENO + 1 # change geno numeric format from -101 to 012

PHENO <- read.delim('data/Intermediate_Files/Phenotypic_Format_wide_2.txt')

#get list of significant markers for traits
LIST <- read.csv('data/Final_Files/output_COMBINED.csv') #master list
#make sig list of each year, this is helpful since phenotype values are stored in by-year files
LIST18 <- LIST[str_detect(LIST$Plot, "2018"),]
LIST19 <- LIST[str_detect(LIST$Plot, "2019"),]
LIST20 <- LIST[str_detect(LIST$Plot, "2020"),]
LIST21O <- LIST[str_detect(LIST$Plot, "OLA_2021"),]
LIST21S <- LIST[str_detect(LIST$Plot, "SAL_2021"),]
LIST22O <- LIST[str_detect(LIST$Plot, "OLA_2022"),]
LIST22S <- LIST[str_detect(LIST$Plot, "SAL_2022"),]
###

#make list with traits for each year
#2018
LIST_2018 <- split(LIST18, LIST18$Trait) # Split data frame to get lists for each chr
trt <- unique(LIST18$Trait) #list of all traits 
for (i in trt){
  temp = LIST_2018[[i]] #table of infor for sig loci for one trait
  new <- subset(geno[, colnames(geno) %in% temp$SNP, drop=FALSE])
  new <- tibble::rownames_to_column(new, "taxa") #convert row names into first column of data
  assign(paste0(i,'_2018_snps'), new, envir=parent.frame()) #save table to environment 
}
#add phenotypic value scores to the geno matrix 
new_str <- gsub('FarmCPU.','',trt)
PHENO2018 <- PHENO[str_detect(PHENO$phenotype_year, "2018"),]
PHENO2018 <- PHENO2018[,c(2,7:ncol(PHENO2018))]
p2018 <- subset(PHENO2018[, colnames(PHENO2018) %in% new_str]) #traits that are in the significant trait list
p2018 <- cbind(PHENO2018[,c(1)], p2018)
colnames(p2018)[1] <- "taxa"



#SPKDEN
SM_SPKDEN_2018 <- p2018[,c("taxa", "SPKDEN")] #get taxa and trait values
SM_SPKDEN_2018 <- left_join(FarmCPU.SPKDEN_2018_snps,SM_SPKDEN_2018, by="taxa") #combine the phenotypic values to genomic matrix
colnames(SM_SPKDEN_2018)[colnames(SM_SPKDEN_2018)=="SPKDEN"] <- "traitval" #trait name
CSigMar_stk<-reshape2::melt(SM_SPKDEN_2018,id=c('traitval','taxa'))
colnames(CSigMar_stk)[colnames(CSigMar_stk)=="variable"] <- "marker"
colnames(CSigMar_stk)[colnames(CSigMar_stk)=="value"] <- "score"
####

# For each marker
markers <- unique(CSigMar_stk$marker)
list(markers)
#SJ01_10162383  SJ01_400972679 SS04_329197898 SS07_134934580 SV01_7732815   SV02_172479310 SV04_80560865 SJ03_118524570 SJ03_598048237

L <- data.frame("Marker"=markers)
E <- data.frame(Effect=NA)

########

doe <- CSigMar_stk[CSigMar_stk$marker == "SJ01_10162383",]
Zero = doe[doe$score == "0",]
One=doe[doe$score == "1",]
Two=doe[doe$score == "2",]

Z = mean(Zero$traitval, na.rm = T)
O = mean(One$traitval, na.rm = T)
D = mean(Two$traitval, na.rm = T)
if(Z=='NaN'){Z <- 0}
if(O=='NaN'){O <- 0}
if(D=='NaN'){D <- 0}

Effect = ((2*Z)+O) - ((2*D)+O)

E <- rbind(E,Effect)
###
doe <- CSigMar_stk[CSigMar_stk$marker == "SJ01_400972679",]
Zero = doe[doe$score == "0",]
One=doe[doe$score == "1",]
Two=doe[doe$score == "2",]

Z = mean(Zero$traitval, na.rm = T) # mean is zero so modify Effect
O = mean(One$traitval, na.rm = T)
D = mean(Two$traitval, na.rm = T)
if(Z=='NaN'){Z <- 0}
if(O=='NaN'){O <- 0}
if(D=='NaN'){D <- 0}

Effect = ((2*D)+O) - ((2*Z)+O) 

E <- rbind(E,Effect)
####
doe <- CSigMar_stk[CSigMar_stk$marker == "SS04_329197898",]
Zero = doe[doe$score == "0",]
One=doe[doe$score == "1",]
Two=doe[doe$score == "2",]

Z = mean(Zero$traitval, na.rm = T)
O = mean(One$traitval, na.rm = T)
D = mean(Two$traitval, na.rm = T)

Effect = ((2*Z)+O) - ((2*D)+O)

E <- rbind(E,Effect)
####
doe <- CSigMar_stk[CSigMar_stk$marker == "SS07_134934580",]
Zero = doe[doe$score == "0",]
One=doe[doe$score == "1",]
Two=doe[doe$score == "2",]

Z = mean(Zero$traitval, na.rm = T)
O = mean(One$traitval, na.rm = T)
D = mean(Two$traitval, na.rm = T)

Effect = ((2*Z)+O) - ((2*D)+O)

E <- rbind(E,Effect)
####
doe <- CSigMar_stk[CSigMar_stk$marker == "SV01_7732815",]
Zero = doe[doe$score == "0",]
One=doe[doe$score == "1",]
Two=doe[doe$score == "2",]

Z = mean(Zero$traitval, na.rm = T)
O = mean(One$traitval, na.rm = T)
D = mean(Two$traitval, na.rm = T)

Effect = ((2*Z)+O) - ((2*D)+O)

E <- rbind(E,Effect)
####
doe <- CSigMar_stk[CSigMar_stk$marker == "SV02_172479310",]
Zero = doe[doe$score == "0",]
One=doe[doe$score == "1",]
Two=doe[doe$score == "2",]

Z = mean(Zero$traitval, na.rm = T)
O = mean(One$traitval, na.rm = T)
D = mean(Two$traitval, na.rm = T)

Effect = ((2*Z)+O) - ((2*D)+O)

E <- rbind(E,Effect)
####
doe <- CSigMar_stk[CSigMar_stk$marker == "SV04_80560865",]
Zero = doe[doe$score == "0",]
One=doe[doe$score == "1",]
Two=doe[doe$score == "2",]

Z = mean(Zero$traitval, na.rm = T)
O = mean(One$traitval, na.rm = T)
D = mean(Two$traitval, na.rm = T)

Effect = ((2*Z)+O) - ((2*D)+O)

E <- rbind(E,Effect)
####
doe <- CSigMar_stk[CSigMar_stk$marker == "SJ03_118524570",]
Zero = doe[doe$score == "0",]
One=doe[doe$score == "1",]
Two=doe[doe$score == "2",]

Z = mean(Zero$traitval, na.rm = T) # zero so modify Effect
O = mean(One$traitval, na.rm = T)
D = mean(Two$traitval, na.rm = T)

Effect = ((2*0)+O) - ((2*D)+O)

E <- rbind(E,Effect)
####
doe <- CSigMar_stk[CSigMar_stk$marker == "SJ03_598048237",]
Zero = doe[doe$score == "0",]
One=doe[doe$score == "1",]
Two=doe[doe$score == "2",]

Z = mean(Zero$traitval, na.rm = T) # zero so modigy effect
O = mean(One$traitval, na.rm = T)
D = mean(Two$traitval, na.rm = T)

Effect = ((2*0)+O) - ((2*D)+O)

E <- rbind(E,Effect)
########
EE <- as.data.frame(E[-1,])
colnames(EE)[1] <- "Effect"
EE$Effect <- abs(EE$Effect)
LEE <- cbind(L,EE)



# Make panneled box plot with:
## SPKYLD 2018  V06_102610885
## FHBDON 2019  S05_89380677
## FHBDON 2021S J07_491374653 

LIST_2018 <- split(LIST18, LIST18$Trait) # Split data frame to get lists for each chr
trt <- unique(LIST18$Trait) #list of all traits 
for (i in trt){
  temp = LIST_2018[[i]] #table of infor for sig loci for one trait
  new <- subset(geno[, colnames(geno) %in% temp$SNP, drop=FALSE])
  new <- tibble::rownames_to_column(new, "taxa") #convert row names into first column of data
  assign(paste0(i,'_2018_snps'), new, envir=parent.frame()) #save table to environment 
}
#add phenotypic value scores to the geno matrix 
new_str <- gsub('FarmCPU.','',trt)
PHENO2018 <- PHENO[str_detect(PHENO$phenotype_year, "2018"),]
PHENO2018 <- PHENO2018[,c(2,7:ncol(PHENO2018))]
p2018 <- subset(PHENO2018[, colnames(PHENO2018) %in% new_str]) #traits that are in the significant trait list
p2018 <- cbind(PHENO2018[,c(1)], p2018)
colnames(p2018)[1] <- "taxa"

## SPKYLD 2018  V06_102610885

#SPKYLD
SM_SPKYLD_2018 <- p2018[,c("taxa", "SPKYLD")] #get taxa and trait values
SM_SPKYLD_2018 <- left_join(FarmCPU.SPKYLD_2018_snps,SM_SPKYLD_2018, by="taxa") #combine the phenotypic values to genomic matrix
colnames(SM_SPKYLD_2018)[colnames(SM_SPKYLD_2018)=="SPKYLD"] <- "traitval" #trait name
CSigMar_stk<-reshape2::melt(SM_SPKYLD_2018,id=c('traitval','taxa'))
colnames(CSigMar_stk)[colnames(CSigMar_stk)=="variable"] <- "marker"
colnames(CSigMar_stk)[colnames(CSigMar_stk)=="value"] <- "score"

CSigMar_stk_NEW <- CSigMar_stk[CSigMar_stk$marker == "SV06_102610885",] # just the marker I want for figure

a <-ggplot(CSigMar_stk_NEW,aes(x=factor(score), y=traitval)) + 
  geom_boxplot(varwidth = TRUE) +
  labs(title="SPKYLD_2018     SV06_102610885", x="Allele", y="Phenotype Rating") + 
  theme(strip.text.x = element_text(size = 6))
rm(SM_SPKYLD_2018,CSigMar_stk, CSigMar_stk_NEW,FarmCPU.SPKYLD_2018_snps)

#2019
LIST_2019 <- split(LIST19, LIST19$Trait) # Split data frame to get lists for each chr
trt <- unique(LIST19$Trait) #list of all traits 
for (i in trt){
  temp = LIST_2019[[i]] #table of infor for sig loci for one trait
  new <- subset(geno[, colnames(geno) %in% temp$SNP, drop=FALSE])
  new <- tibble::rownames_to_column(new, "taxa") #convert row names into first column of data
  assign(paste0(i,'_2019_snps'), new, envir=parent.frame()) #save table to environment 
}
#add phenotypic value scores to the geno matrix 
new_str <- gsub('FarmCPU.','',trt)
PHENO2019 <- PHENO[str_detect(PHENO$phenotype_year, "2019"),]
PHENO2019 <- PHENO2019[,c(2,7:ncol(PHENO2019))]
p2019 <- subset(PHENO2019[, colnames(PHENO2019) %in% new_str]) #traits that are in the significant trait list
p2019 <- cbind(PHENO2019[,c(1)], p2019)
colnames(p2019)[1] <- "taxa"

## FHBDON 2019  S05_89380677

#FHBDON
SM_FHBDON_2019 <- p2019[,c("taxa", "FHBDON")] #get taxa and trait values
SM_FHBDON_2019 <- left_join(FarmCPU.FHBDON_2019_snps,SM_FHBDON_2019, by="taxa") #combine the phenotypic values to genomic matrix
colnames(SM_FHBDON_2019)[colnames(SM_FHBDON_2019)=="FHBDON"] <- "traitval" #trait name
CSigMar_stk<-reshape2::melt(SM_FHBDON_2019,id=c('traitval','taxa'))
colnames(CSigMar_stk)[colnames(CSigMar_stk)=="variable"] <- "marker"
colnames(CSigMar_stk)[colnames(CSigMar_stk)=="value"] <- "score"

CSigMar_stk_NEW <- CSigMar_stk[CSigMar_stk$marker == "SS05_89380677",]

b <-ggplot(CSigMar_stk_NEW,aes(x=factor(score), y=traitval)) + 
  geom_boxplot(varwidth = TRUE) +
  labs(title="FHBDON_2019     S05_89380677", x="Allele", y="Phenotype Rating") + 
  theme(strip.text.x = element_text(size = 6))
rm(SM_FHBDON_2019,CSigMar_stk,CSigMar_stk_NEW,FarmCPU.FHBDON_2019_snps)

#2021 SAL
LIST_2021S <- split(LIST21S, LIST21S$Trait) # Split data frame to get lists for each chr
trt <- unique(LIST21S$Trait) #list of all traits 
for (i in trt){
  temp = LIST_2021S[[i]] #table of infor for sig loci for one trait
  new <- subset(geno[, colnames(geno) %in% temp$SNP, drop=FALSE])
  new <- tibble::rownames_to_column(new, "taxa") #convert row names into first column of data
  assign(paste0(i,'_2021S_snps'), new, envir=parent.frame()) #save table to environment 
}
#add phenotypic value scores to the geno matrix 
new_str <- gsub('FarmCPU.','',trt)
PHENO2021S <- PHENO[str_detect(PHENO$Site, "SAL"),]
PHENO2021S <- PHENO2021S[str_detect(PHENO2021S$phenotype_year, "2021"),]
PHENO2021S <- PHENO2021S[,c(2,7:ncol(PHENO2021S))]
p2021S <- subset(PHENO2021S[, colnames(PHENO2021S) %in% new_str]) #traits that are in the significant trait list
p2021S <- cbind(PHENO2021S[,c(1)], p2021S)
colnames(p2021S)[1] <- "taxa"

## FHBDON 2021S J07_491374653 

#FHBDON
SM_FHBDON_2021S <- p2021S[,c("taxa", "FHBDON")] #get taxa and trait values
SM_FHBDON_2021S <- left_join(FarmCPU.FHBDON_2021S_snps,SM_FHBDON_2021S, by="taxa") #combine the phenotypic values to genomic matrix
colnames(SM_FHBDON_2021S)[colnames(SM_FHBDON_2021S)=="FHBDON"] <- "traitval" #trait name
CSigMar_stk<-reshape2::melt(SM_FHBDON_2021S,id=c('traitval','taxa'))
colnames(CSigMar_stk)[colnames(CSigMar_stk)=="variable"] <- "marker"
colnames(CSigMar_stk)[colnames(CSigMar_stk)=="value"] <- "score"

CSigMar_stk_NEW <- CSigMar_stk[CSigMar_stk$marker == "SJ07_491374653",]

c <-ggplot(CSigMar_stk_NEW,aes(x=factor(score), y=traitval)) + 
  geom_boxplot(varwidth = TRUE) +
  labs(title="FHBDON_2021_Salina     J07_491374653", x="Allele", y="Phenotype Rating") + 
  theme(strip.text.x = element_text(size = 6))
rm(SM_FHBDON_2021S,CSigMar_stk,CSigMar_stk_NEW,FarmCPU.FHBDON_2021S_snps)

pdf("data/Figures/BoxPlots_Allelic_Mean_for_Paper.pdf", width = 10, height = 7)
a
b
c
dev.off()

figure <- ggarrange(a, b, c,
                    #labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
figure

#Export to figures file as BoxPlots_Allelic_Mean_for_Paper2.pdf
## Dimensions : 15.69 x 4.68

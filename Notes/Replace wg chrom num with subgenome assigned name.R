#Change chrom numeric notation back to chromosome names
#Merge chrom num to chrom name -- AM1032

rm(list=ls(all=TRUE)) #clears variables
setwd("C:\\Users\\kturner\\utils\\output\\hitsblastout")

#geno <- read.table("iwg+wdIWGSC_7094Unk100pi+D.txt", header=TRUE,as.is=TRUE)
geno <- read.table("iwg+wdIWGSC_6216Unk+D100pi.txt", header=TRUE,as.is=TRUE)

geno[c(1:5),c(1:5)]
x <- geno[,c(2)] #to replace values in column 4

x[x == "Chr01"] <- "TI_1J"
x[x == "Chr06"] <- "TI_2J"
x[x == "Chr07"] <- "TI_3J"
x[x == "Chr10"] <- "TI_4J"
x[x == "Chr14"] <- "TI_5J"
x[x == "Chr18"] <- "TI_6J"
x[x == "Chr20"] <- "TI_7J"
x[x == "Chr02"] <- "TI_1S"
x[x == "Chr04"] <- "TI_2S"
x[x == "Chr08"] <- "TI_3S"
x[x == "Chr11"] <- "TI_4S"
x[x == "Chr13"] <- "TI_5S"
x[x == "Chr17"] <- "TI_6S"
x[x == "Chr21"] <- "TI_7S"
x[x == "Chr03"] <- "TI_1Vs"
x[x == "Chr05"] <- "TI_2Vs"
x[x == "Chr09"] <- "TI_3Vs"
x[x == "Chr12"] <- "TI_4Vs"
x[x == "Chr15"] <- "TI_5Vs"
x[x == "Chr16"] <- "TI_6Vs"
x[x == "Chr19"] <- "TI_7Vs"

#add TA to triticum aestivum chrom
y <- gsub("chr", "TA_", x)

genox<-cbind(geno[,c(1)],y,geno[,c(3:ncol(geno))])
genox[c(1:5),c(1:5)]

colnames(genox)[1] <- "contig"
colnames(genox)[2] <- "chrom"

#order the sequences by chromosome and position on chromosome
ord <- genox[order(genox$chrom, genox$subject_start),]

#write.table(genox, "iwg+wdIWGSC_7094Unk100pi+D.txt", sep="\t") 
write.table(ord, "iwg+wdIWGSC_6216Unk+D100pi_named.txt", sep="\t") 


#this was to match the mismapped file, but not able to match without whole field match
geno[c(1:5),c(1:5)]
x <- geno[,c(2)] #to replace values in column 4

x <- mm[,c(1)] #to replace values in column 4
x[x	=="*ti1"]	<-	"*Chr01-1J"
x[x	==	"*ti6"]	<-	"*Chr06_2J"
x[x	==	"*ti7"]	<-	"*Chr07_3J"
x[x	==	"*ti10"]	<-	"*Chr10_4J"
x[x	==	"*ti14"]	<-	"*Chr14_5J"
x[x	==	"*ti18"]	<-	"*Chr18_6J"
x[x	==	"*ti20"]	<-	"*Chr20_7J"
x[x	==	"*ti2"]	<-	"*Chr02_1S"
x[x	==	"*ti4"]	<-	"*Chr04_2S"
x[x	==	"*ti8"]	<-	"*Chr08_3S"
x[x	==	"*ti11"]	<-	"*Chr11_4S"
x[x	==	"*ti13"]	<-	"*Chr13_5S"
x[x	==	"*ti17"]	<-	"*Chr17_6S"
x[x	==	"*ti21"]	<-	"*Chr21_7S"
x[x	==	"*ti3"]	<-	"*Chr03_1VS"
x[x	==	"*ti5"]	<-	"*Chr05_2VS"
x[x	==	"*ti9"]	<-	"*Chr09_3VS"
x[x	==	"*ti12"]	<-	"*Chr12_4VS"
x[x	==	"*ti15"]	<-	"*Chr15_5VS"
x[x	==	"*ti16"]	<-	"*Chr16_6VS"
x[x	==	"*ti19"]	<-	"*Chr19_7VS"

mmx<-cbind(x,mm[,c(2:ncol(mm))])
mmx[c(1:5),c(1:5)]
colnames(genox)[1] <- "Name"

#this might work...
#gsub("cheap", "sheep's", "A wolf in cheap clothing")
mm <- read.table("C:\\Users\\kturner\\utils\\output\\seq\\wg+wd\\bc1\\bc1n_mismapped.txt",header=T, as.is=TRUE)

mm$Name <- gsub("ti1$", "Chr01-1J", mm$Name)
mm$Name <- gsub("ti6$", "Chr06_2J", mm$Name)
mm$Name <- gsub("ti7$", "Chr07_3J", mm$Name)
mm$Name <- gsub("ti10$", "Chr10_4J", mm$Name)
mm$Name <- gsub("ti14$", "Chr14_5J", mm$Name)
mm$Name <- gsub("ti18$", "Chr18_6J", mm$Name)
mm$Name <- gsub("ti20$", "Chr20_7J", mm$Name)
mm$Name <- gsub("ti2$", "Chr02_1S", mm$Name)
mm$Name <- gsub("ti4$", "Chr04_2S", mm$Name)
mm$Name <- gsub("ti8$", "Chr08_3S", mm$Name)
mm$Name <- gsub("ti11$", "Chr11_4S", mm$Name)
mm$Name <- gsub("ti13$", "Chr13_5S", mm$Name)
mm$Name <- gsub("ti17$", "Chr17_6S", mm$Name)
mm$Name <- gsub("ti21$", "Chr21_7S", mm$Name)
mm$Name <- gsub("ti3$", "Chr03_1VS", mm$Name)
mm$Name <- gsub("ti5$", "Chr05_2VS", mm$Name)
mm$Name <- gsub("ti9$", "Chr09_3VS", mm$Name)
mm$Name <- gsub("ti12$", "Chr12_4VS", mm$Name)
mm$Name <- gsub("ti15$", "Chr15_5VS", mm$Name)
mm$Name <- gsub("ti16$", "Chr16_6VS", mm$Name)
mm$Name <- gsub("ti19$", "Chr19_7VS", mm$Name)






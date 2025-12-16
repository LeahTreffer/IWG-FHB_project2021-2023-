# Add information about Major Allele to the Phenotypic Mean table
# add ref and allele columns then also find and include which of those alleles is minor (0) and major (2)
# used in 2.7_Phenotypic_Means_2025.Rmd

# Function
freqphenomean <- function(frequency_table,trait_name,mean_table,table_name){
  #Get same list of SNPs in both tables then merge
  bread <- frequency_table[frequency_table$trait == trait_name, ]
  bread <- semi_join(bread, mean_table, by="SNP") #shorten frequencies table to the same SNPs that are in the phenotypic mean table 
  bread <- merge(mean_table, bread, by="SNP") #merge
  
  MinorAllele <- bread[,c("SNP", "Ref", "alleles", "allelic_frequency.A", "allelic_frequency.T", "allelic_frequency.G", "allelic_frequency.C")] #subset allele frequencies
  MinorAllele$allele1 <- gsub("/.*", "", MinorAllele$alleles) #Remove after first allele
  MinorAllele$allele2 <- gsub(".*/","",MinorAllele$alleles) #Remove all before and up to "/"
  MinorAllele$ALT <- ifelse(MinorAllele$Ref == MinorAllele$allele1, MinorAllele$allele2, MinorAllele$allele1)
  
  df1 <- MinorAllele %>% filter_all(any_vars(. %in% c(1))) #list of monomorphic snps
  
  if(nrow(df1) == 0){
    df2 <- MinorAllele
    df2[df2 == "0"] <- "NA" #replace 0 with NA so that they aren't seen as then min number
    df2$minor_allele <- names(df2)[apply(df2, MARGIN = 1, FUN = which.min)]
    ma <- as.data.frame(gsub("allelic_frequency.","",df2$minor_allele)) #remove character string leaving just the allele
    colnames(ma)[1]<- "ma"
    df2 <- cbind(df2,ma)
    df2 <- df2[,c("SNP", "ALT", "allelic_frequency.A", "allelic_frequency.T", "allelic_frequency.G", "allelic_frequency.C", "ma")]
    bready <- merge(bread, df2[,c("SNP", "ALT", "ma")], by="SNP", all=TRUE)
    colnames(bready)[2] <- "PhentypicMeanMinorHomozygous"
    colnames(bready)[3] <- "PhentypicMeanHeterozygous"
    colnames(bready)[4] <- "PhentypicMeanMajorHomozygous"
    colnames(bready)[12] <- "PVE"
    colnames(bready)[42] <- "MinorAllele"
    
    bready <- bready[,c("SNP", "Chr", "Pos", "alleles", "Ref", "ALT", "MinorAllele", "trait", "plot", "P.value", "MAF", "nobs", "Effect", "H.B.P.Value", "PVE", "GG", "AG", "AA", "CG", "CC", "TT", "GT", "CT", "AT", "AC", "genotypic_frequency.CT", "genotypic_frequency.CC", "genotypic_frequency.TT", "genotypic_frequency.AT", "genotypic_frequency.AA", "genotypic_frequency.AG", "genotypic_frequency.GG", "genotypic_frequency.CG", "genotypic_frequency.AC", "genotypic_frequency.GT", "allelic_frequency.A", "allelic_frequency.T", "allelic_frequency.G", "allelic_frequency.C", "PhentypicMeanMinorHomozygous", "PhentypicMeanHeterozygous", "PhentypicMeanMajorHomozygous")] # reorder
    
    assign(paste0(table_name), bready, envir=parent.frame())
    
  } else {
    df2 <- setdiff(MinorAllele, df1) #list of normal snps
    
    df1[df1 == "0"] <- "NA" #replace 0 with NA so that they aren't seen as then min number
    df1$minor_allele <- names(df1)[apply(df1, MARGIN = 1, FUN = which.min)]
    df1$minor_allele <- gsub("allelic_frequency.","",df1$minor_allele) #Remove all before and up to "/"
    df1$ma <- ifelse(df1$minor_allele == df1$Ref, df1$ALT, df1$Ref)
    df1 <- df1[,c("SNP", "ALT", "allelic_frequency.A", "allelic_frequency.T", "allelic_frequency.G", "allelic_frequency.C", "ma")]
    
    df2[df2 == "0"] <- "NA" #replace 0 with NA so that they aren't seen as then min number
    df2$minor_allele <- names(df2)[apply(df2, MARGIN = 1, FUN = which.min)]
    ma <- as.data.frame(gsub("allelic_frequency.","",df2$minor_allele)) #remove character string leaving just the allele
    colnames(ma)[1]<- "ma"
    df2 <- cbind(df2,ma)
    df2 <- df2[,c("SNP", "ALT", "allelic_frequency.A", "allelic_frequency.T", "allelic_frequency.G", "allelic_frequency.C", "ma")]
    
    df <- rbind(df2, if(exists("df1")) df1)
    
    bready <- merge(bread, df[,c("SNP","ALT", "ma")], by="SNP", all=TRUE)
    colnames(bready)[2] <- "PhentypicMeanMinorHomozygous"
    colnames(bready)[3] <- "PhentypicMeanHeterozygous"
    colnames(bready)[4] <- "PhentypicMeanMajorHomozygous"
    colnames(bready)[12] <- "PVE"
    colnames(bready)[42] <- "MinorAllele"
    
    bready <- bready[,c("SNP", "Chr", "Pos", "alleles", "Ref", "ALT", "MinorAllele", "trait", "plot", "P.value", "MAF", "nobs", "Effect", "H.B.P.Value", "PVE", "GG", "AG", "AA", "CG", "CC", "TT", "GT", "CT", "AT", "AC", "genotypic_frequency.CT", "genotypic_frequency.CC", "genotypic_frequency.TT", "genotypic_frequency.AT", "genotypic_frequency.AA", "genotypic_frequency.AG", "genotypic_frequency.GG", "genotypic_frequency.CG", "genotypic_frequency.AC", "genotypic_frequency.GT", "allelic_frequency.A", "allelic_frequency.T", "allelic_frequency.G", "allelic_frequency.C", "PhentypicMeanMinorHomozygous", "PhentypicMeanHeterozygous", "PhentypicMeanMajorHomozygous")] # reorder
    
    assign(paste0(table_name), bready, envir=parent.frame())
  }
  
  
}

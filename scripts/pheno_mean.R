# get phenotypic mean of each genotype (0,1,2)
# used in 2.4_Phenotypic_Means.Rmd and 2.7_Phenotypic_Means_2025
#function
pheno_mean <- function(trait_table,snp_table,output_name){
  traitt <- data.frame(matrix(ncol=4,nrow=0)) #empty dataframe
  colnames(traitt) <- c('SNP', 'mean0', 'mean1', 'mean2') #column header names
  for (l in snp_table){
    Hom0 <- trait_table[trait_table[, l] == '0', ]
    mean0 <- as.data.frame(mean(Hom0[, 2], na.rm=TRUE))
    row.names(mean0) <- l
    colnames(mean0) <- 'mean0'
    mean0 <- tibble::rownames_to_column(mean0, "SNP")
    
    Het1 <- trait_table[trait_table[, l] == '1', ]
    mean1 <- as.data.frame(mean(Het1[, 2], na.rm=TRUE))
    row.names(mean1) <- l
    colnames(mean1) <- 'mean1'
    mean1 <- tibble::rownames_to_column(mean1, "SNP")
    
    Hom2 <- trait_table[trait_table[, l] == '2', ]
    mean2 <- as.data.frame(mean(Hom2[, 2], na.rm=TRUE))
    row.names(mean2) <- l
    colnames(mean2) <- 'mean2'
    mean2 <- tibble::rownames_to_column(mean2, "SNP")
    
    means <- (list(mean0,mean1,mean2) %>% reduce(inner_join, by='SNP'))
    traitt <- rbind(traitt, means)
  }
  assign(paste0(output_name),traitt,envir=parent.frame())
}
#requires: library(dplyr) 
#columns you MUST have: SNP, Plot, CHR, Pos
  #SNP: character column identification for accessions (ex.markers/snp), this is the column that duplicates are looked for in
  #Plot: character column for observation group/factor (ex. cycle, field site, year) used when arranging table, don't need to have more than one but needs to have data in column 
  #CHR: character column with the chromosome name of the chromosome each marker is a part of; table is subset by each of these groups and positions within this group are compared to each other
  #Pos: numeric column of position of each marker/snp
#additional columns are fine to have, but will not be looked at by the function
#these columns can be in any order
#optional arguments run when =TRUE:
  #order: order/arrange output list  by CHR
  #space: add a blank row after each pair for visual aid in seeing matched rows
#output will have a last column called 'pair'. this is a visual way to find the groups that met the distance requirement 

#example table to test function: 
#list4 = data.table(
  #SNP = c("SJ01_320309248","SJ01_320309278","SJ01_320109248","SS05_83112789","SS05_83112769","SS05_83912769","SV07_590255485","SV07_590255488","SV07_590255488","SV07_590255482","SV07_590255482","SV07_590255482","SV07_590255482","SV07_590255482"),
  #traits = c("FarmCPU.FHBDON","FarmCPU.FHBDON","FarmCPU.PTHT","FarmCPU.FHBD3G.DON","FarmCPU.FHBD3G.DON","FarmCPU.STMANG","FarmCPU.FHBD3G","FarmCPU.FHBD3G","FarmCPU.FHBD3G","FarmCPU.FHBSPKINC","FarmCPU.FHBSPKINC","FarmCPU.FHBSPKINC","FarmCPU.FHBSPKINC","FarmCPU.FHBSPKINC"),
  #Plot = c("C07_SAL_2018","C07_SAL_2019","C07_SAL_2020","C07_SAL_2018","C07_SAL_2019","C07_SAL_2020","C07_SAL_2018","C07_SAL_2019","C07_SAL_2019","C07_SAL_2020","C07_SAL_2020","C07_SAL_2020","C07_SAL_2020","C07_SAL_2020"),
  #CHR = c("J01","J01","J01","S05","S05","S05","V07","V07","V07","V07","V07","V07","V07","V07"),
  #Pos = c("320309248","320309278","320109248","83112789","83112769","83912769","590255485","590255488","590255488","590255482","590255482","590255482","590255482","590255482"))
#list4$Pos <- as.numeric(list4$Pos)

similar.region <- function(snp_table, range){
  snp_list <- split(snp_table, snp_table$CHR) # Split data frame to get lists for each chr
  chr <- unique(snp_table$CHR) #list of all 21 chromosomes (J01-J07,S01-S07,V01-V07) 
  z = range #assign the numeric desired range (distance +/- from SNP) to the character z 
  
  for (i in chr){
    match_list <- data.frame() #empty table for each chromosome to add any rows with snps within desired range of each other   
    temp = snp_list[[i]] #table of one chromosome info
    
    exact_list <- temp[duplicated(temp$SNP)|duplicated(temp$SNP, fromLast=TRUE),] #Extract duplicate records in dataframe
    exact_list <- arrange(exact_list, SNP, Plot) #order/arrange list first by SNP then Plot
    
    for(j in 1:nrow(temp)){
      for(k in 1:nrow(temp)){
        u = temp[j,] #one row 
        v = temp[k,] #second row
      
      comp <- rbind(u,v) #bind two rows 
      
      f = u$Pos - v$Pos #distance between two rows
      
      if (f > 0 && f <= z){
        match_list <- rbind(match_list,comp) #bind compared rows into saved table if the difference between their position was within the desired range 
      }
      }
    }
    
    b = as.numeric(nrow(exact_list)) #get number of rows 
    if (b > 0) {
      exact_list <- exact_list %>% group_by(Pos) %>% mutate(id=cur_group_id()) #assign numeric id for groups of identical numbers in Pos column of exact_list
      exact_list$pair = chartr("123456789", "abcdefghi", exact_list$id) #character to letter, assign letter to numeric
      exact_list <- select(exact_list, -c(id)) #remove id column 
    } else {
      #columns = c("SNP","traits","Plot","CHR","Pos","pair") 
      #exact_list = data.frame(matrix(nrow = 0, ncol = length(columns))) 
      #colnames(exact_list) = columns
      exact_list <- data.frame()
    }
    
    c = as.numeric(nrow(match_list)) #get number of rows
    if (c > 0){
      cc = c/2 #divide number of rows in half
      d = as.numeric(c(1:cc)) #make list numeric 
      e = chartr("123456789", "ABCDEFGHI", d) #character to letter, assign letter to numeric from the list 'd'
      match_list$pair <- rep(e,each=2) #add column to output that gives letter to each match set which are binded in the list in pairs already
    } else {
      #columns = c("SNP","traits","Plot","CHR","Pos","pair") 
      #match_list = data.frame(matrix(nrow = 0, ncol = length(columns))) 
      #colnames(match_list) = columns
      match_list <- data.frame()
    }
    final_list <- rbind(exact_list,match_list)
    assign(paste0(i,'_qtl_list'), final_list, envir=parent.frame()) #save table for chromosome to environment 
  }
}

#similar.region(list4, 100) # run function


#similar.region <- function(snp_table, snp_list, lower_bound, upper_bound){
#y = lower_bound #assign the numeric lower bound of desired range to the character y 
#if (f < 0 && f >= y || f > 0 && f <= z){
#similar.region <- function(list4, my_splits3, -100, 100){

qtl.comparison <- function(snp_table, range, order=NULL, space=NULL){
  exact_list <- snp_table[duplicated(snp_table$SNP)|duplicated(snp_table$SNP, fromLast=TRUE),] #Extract duplicate records in dataframe
  exact_list <- arrange(exact_list, SNP, Plot) #order/arrange list first by SNP then Plot
    b = as.numeric(nrow(exact_list)) #get number of rows 
  if (b > 0) {
    exact_list <- exact_list %>% group_by(SNP) %>% mutate(id=cur_group_id()) #assign numeric id for groups of identical numbers in SNP column of exact_list
    exact_list$p = chartr("123456789", "abcdefghi", exact_list$id) #character to letter, assign letter to numeric
    exact_list <- select(exact_list, -c(id)) #remove id column 
  } else {
    exact_list <- data.frame()
  }

  snp_table$p <- NA #add empty column so that this table can be bound to the exact match table
  snp_list <- split(snp_table, snp_table$CHR) # Split data frame to get lists for each chr
  
    chr <- unique(snp_table$CHR) #list of all chromosomes in table (J01-J07,S01-S07,V01-V07) 
    z = range #assign the numeric desired range (distance +/- from SNP) to the character z 
    
    for (i in chr){
      match_list <- data.frame() #empty table for each chromosome to add any rows with snps within desired range of each other   
      temp = snp_list[[i]] #table of one chromosome info
      
      for(j in 1:nrow(temp)){
        for(k in 1:nrow(temp)){
          first_row = temp[j,] #one row 
          second_row = temp[k,] #second row
          
          comp <- rbind(first_row,second_row) #bind two rows 
          
          f = first_row$Pos - second_row$Pos #distance between two rows
          
          if (f > 0 && f <= z){
            match_list <- rbind(match_list,comp) #bind compared rows into saved table if the difference between their position was within the desired range 
          }
        }
      }
      exact_list <- rbind(exact_list,match_list) #add the match list of the chromosome to the exact list, which then becomes a running list of all snps as more lists from other chromosomes are added each iteration
    }
    perf_list <- exact_list[!is.na(exact_list$p),] #extract rows with a value other than NA
    colnames(perf_list)[colnames(perf_list) == "p"] ="pair" #rename p column to pair
    
    close_list <- exact_list[is.na(exact_list$p),] #extract rows with NA
    c = as.numeric(nrow(close_list)) #get number of rows 
    if (c > 0) {
      cc = c/2 #divide number of rows in half
      d = as.numeric(c(1:cc)) #make list numeric 
      e = chartr("123456789", "ABCDEFGHI", d) #character to letter, assign letter to numeric from the list 'd'
      close_list$pair <- rep(e,each=2) #add column to output that gives letter to each match set which are binded in the list in pairs already
      close_list <- select(close_list, -c(p)) #remove id column 
    } else {
      colnames(close_list)[colnames(close_list) == "p"] ="pair"
    }
    full_list <- rbind(perf_list,close_list) #rebind table
    if (order == TRUE) {full_list <- arrange(full_list, CHR, pair)} #optional argument order/arrange list first by CHR then pair, default is FALSE
    if (space == TRUE) {data_new <- as.data.frame(lapply(full_list, as.character), stringsAsFactors = FALSE)
    data_new <- head(do.call(rbind, by(data_new, full_list$pair, rbind, "")), -1)
    full_list <- data_new} #optional argument to add a blank row after each pair for visual aid in seeing matched rows, default is FALSE
    if (order == TRUE && space == TRUE){
      full_list$rip <- paste0(full_list$CHR,'_',full_list$pair)
      full_list <- arrange(full_list, rip)
      data_new <- as.data.frame(lapply(full_list, as.character), stringsAsFactors = FALSE)
      data_new <- head(do.call(rbind, by(data_new, full_list$rip, rbind, "")), -1)
      full_list <- data_new[data_new$rip != "_", ]
      full_list <- subset(full_list, select=-c(rip))
      rownames(full_list) <- NULL
      full_list <- full_list[-c(1),]
    }
    assign(paste0('qtl_list'), full_list, envir=parent.frame()) #save table of exact matches and snps within desired range to environment as table called 'qtl_list'
}

#qtl.comparison(list4, 100, order=TRUE, space=TRUE)

library(tidyr)
library(dplyr)
library(stringr)
library(graphics)
library(stats)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(janitor)
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
install.packages("Hmisc")
library(Hmisc)
install.packages("corrplot")
library(corrplot)
library(ggplot2)
library(psych)
library(data.table)
library(reshape2)

#Histograms
# Make a histogram for each trait column in a table
# input: phenotype file (table) with taxa in first column followed by columns of for each trait
# input: list of traits to loop through 


# Histogram with Normal Curve
# input: list of traits to loop through 
# input: phenotype file (table) of columns of for each trait
Histogram <- function(trait_list, pheno_table) {
  for (i in trait_list){
    i = paste0(i)
    
    trait <- na.omit(pheno_table[i])
    trait <- sapply(trait, as.numeric)
    
    h <- hist(trait, main=i, breaks=40)
    xfit<-seq(min(trait),max(trait),length=40) 
    yfit<-dnorm(xfit,mean=mean(trait),sd=sd(trait)) 
    yfit <- yfit*diff(h$mids[1:2])*length(trait) 
    lines(xfit, yfit, col="blue", lwd=2)
  }
  
}

# Histogram with Normal Curve
# Use for comparing IWG to Wheat in same graph
# input: directory path for output pdf
# input: list of traits to loop through 
# input: phenotype file (table) of columns of for each trait
HistogramOverlap <- function(output_name, trait_list, pheno_table) {
  plot_list = list()
  for (i in trait_list){
    i = paste0(i)
    
    table <- pheno_table[,c('entry',i)]
    trait <- na.omit(table)
    trait <- as.data.frame(trait)
    colnames(trait)[2] <- "values"
  
  #plot <- ggplot(trait,aes(x=values)) + 
    #geom_histogram(data=subset(trait,entry == 'IWG'),aes(fill = "orange"), alpha = 0.7) +
    #geom_histogram(data=subset(trait,entry == 'Wheat'),aes(fill = "cornflowerblue"), alpha = 0.7)+
  #scale_fill_manual(name = "entry", values = c("orange", "cornflowerblue"), labels=c("IWG", "Wheat"))+
    #ggtitle(i) # for the main title
    
    plot <- ggplot(trait,aes(x=values)) + 
    geom_histogram(data=subset(trait,entry == 'IWG'), fill = "orange", alpha = 0.7) +
    geom_histogram(data=subset(trait,entry == 'Wheat'), fill = "cornflowerblue", alpha = 0.7)+
    scale_color_manual(values = c(IWG = "orange", Wheat = "cornflowerblue"))+
      ggtitle(i) # for the main title
    

    plot_list[[i]] = plot
    
    pdf(output_name)
    for (k in trait_list) {
      print(plot_list[[k]])
    }
    dev.off()
  }
}

# Histograms + Kernel Density Plots
# input : path to save to that includes name of the table 
# input: list of traits to loop through 
# input: phenotype file (table) of columns of for each trait 
Histogram_Kernel_Density <- function(table_name, trait_list, pheno_table){
  plot_list = list()
  for (j in trait_list){
    
    plot <- ggplot(pheno_table, aes_string(x=pheno_table[[j]])) +
      geom_histogram(aes(y=..density..), colour="black", fill="white") +
      geom_density(alpha=.2, fill="#FF6666") +
      ggtitle(j)
    
    plot_list[[j]] = plot
    
    pdf(table_name)
    for (k in trait_list) {
      print(plot_list[[k]])
    }
    dev.off()
  }
}

#Normalized Histograms + Kernel Density Plots
# input: path to save to that includes name of the table 
# input: list of traits to loop through 
# input: phenotype file (table) of columns of for each trait 
N_Histogram_Kernel_Density <- function(table_name, trait_list, pheno_table){
  plot_list = list()
  for (l in trait_list){
    
    plot <- ggplot(pheno_table, aes_string(x=pheno_table[[l]])) +
      geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = .15) +
      geom_density(alpha=.2, fill="#FF6666") +
      ggtitle(l)
    
    plot_list[[l]] = plot
    
    pdf(table_name)
    for (m in trait_list) {
      print(plot_list[[m]])
    }
    dev.off()
  }
}

#Compare Raw to Normalized Plots
# input: path to save to that includes name of the table 
# input: list of traits to loop through 
# input: phenotype file (table) of columns of for each trait, raw ratings
# input: phenotype file (table) of columns of for each trait, log transformed ratings
Raw_Normalized <- function(table_name, trait_list, r_pheno, n_pheno){
  plot_list = list()
  for (n in trait_list){
    n <- paste0(n)
    
    one <- as.data.frame(r_pheno[[n]])
    two <- as.data.frame(n_pheno[[n]])
    one['marker'] = 'raw'
    two['marker'] = 'normalizied'
    colnames(one)[1] <- "score"
    colnames(two)[1] <- "score"
    three <- rbind(one,two) #bind into new table
    three <- three[!(is.na(three$score) | three$score==""), ] #remove rows with any blank cell
    three$score <- as.numeric(three$score) #as numeric 
    
    plot <- ggplot(three, aes(x=score)) +
      geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = .15) +
      geom_density(alpha=.2, fill="#FF6666")+
      facet_wrap(~marker, scales = "free", ncol = 2) +
      labs(title=n) + 
      theme(strip.text.x = element_text(size = 6))
    
    plot_list[[n]] = plot
    
    pdf(table_name)
    for (o in trait_list) {
      print(plot_list[[o]])
    }
    dev.off()
    
  }
}


#Correlation Plots
#http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#Box Plots
# input: phenotype file (table) with taxa in first column followed by columns of for each trait
# input: list of traits to loop through 
BoxPlots <- function(table_name, trait_list, pheno_table){
  plot_list = list()
  for (j in trait_list){
    
    j = paste0(j)
    main = paste0('Distribution of',' ',j,' ','Ratings')
    
    pheno = pheno_table[,c("Year","plant_id","Cycle","Site","Time",j)]
    colnames(pheno)[6] <- "Value"
    pheno <- na.omit(pheno) #remove rows with NA for Value
    pheno$Values <- pheno$Value
    pheno$Year= as.character(pheno$Year)
    pheno$Cycle= as.factor(as.character(pheno$Cycle))
    #pheno$Values= as.integer(as.factor(pheno$Values)) 
    
    plot <- ggplot(pheno, aes(x=Year, y=Values, fill=Cycle))+
      geom_boxplot()+
      labs(title=main)
    
    plot_list[[j]] = plot
    
    pdf(table_name)
    for (k in trait_list) {
      print(plot_list[[k]])
    }
    dev.off()
  }
}

#Needs work if we want to use it
BoxPlot <- function(table_name, trait_list, pheno_table){
  plot_list = list()
  for (j in trait_list){
    
    j = paste0(j)
    main = paste0('Average Cycle Differences for',' ',j)
    
    pheno = pheno_table[,c("Year","plant_id","Cycle","Site","Time",j)]
    colnames(pheno)[6] <- "Value"
    pheno <- na.omit(pheno) #remove rows with NA for Value
    pheno$Values <- pheno$Value
    pheno$Year= as.character(pheno$Year)
    pheno$Cycle= as.factor(as.character(pheno$Cycle))
    pheno$Values= as.integer(as.factor(pheno$Values)) 
    
    plot <- ggplot(table, aes(x=Cycle, y=Values, fill=Cycle))+
      geom_boxplot()+
      labs(title=main)
    
    plot_list[[j]] = plot
    
    pdf(table_name)
    for (k in trait_list) {
      print(plot_list[[k]])
    }
    dev.off()
  }
}

QQplot.save <- function(data.scan, trait.list, plot.title){
  plot_list = list()
  for (i in trait.list){
    i <- paste0(i)
    
    title <- paste0(plot.title,'_',i)
    table_path <- paste0("output/rrBLUP/QQplot_", plot.title)
    
    plot <- qq.plot(data.scan,trait=i) + ggtitle(label=title)
    plot_list[[i]] = plot
    
    pdf(table_path, width = 12, height = 8.91)
    for (k in trait.list) {
      print(plot_list[[k]])
    }
    dev.off()
  }
}

pdf('output/rrBLUP/QQplot_Original_C07.pdf', width = 12, height = 8.91)
QQplot.save(data.original.scan, trait.list, "Original_C07")
dev.off()

QQplot.save <- function(data.scan, trait.list, plot.title) {
  for (i in trait.list){
    i = paste0(i)

    title <- paste0(plot.title,'_',i)
    
    plot <- qq.plot(data.scan,trait=i) + ggtitle(label=title)
    
  }
  
}


# Plots <- c("C07_SAL_2018", "C07_SAL_2019", "C07_SAL_2020", "C10_SAL_2021", "C10_OLA_2021", "C10_SAL_2022", "C10_OLA_2022")
# Traits <- c("FHBDON", "FHBD3G", "FHBDISIND", "FHBSPKINC")
# List <- list(C07_SAL_2018 = C7_2018_S, C07_SAL_2019 = C7_2019_S, C07_SAL_2020 =C7_2020_S, C10_SAL_2021 = C10_2021_S, C10_OLA_2021 = C10_2021_O, C10_SAL_2022 = C10_2022_S, C10_OLA_2022 = C10_2022_O)
# path <- 'data/Figures/BarPlot_FHB.pdf'
# bar_plot_fhb(Plot,Traits,List,path)
#Plots: vector of the names of each plot 
#Traits: vector of the names of traits to find and use from each table
#List: list containing each year's table, named with same names as in the plots vector, needs to be in long format and contain plant ids, trait ids, and trait values. there can be more trait ids than there is in Traits vector but those that are in both must be exact matches
#path: vector of output path, used to save plots as a pdf

#if adding traits, will have to add in other ifelse statements to customize the scale
####################
bar_plot_fhb <- function(Plots,Traits,List,path){
  plot_list = list()
  for (b in Traits){
    b = paste0(b)
    bdf <- data.frame()
    for (a in Plots){
      temp = List[[a]] 
      test <- as_tibble(temp[c("plant_id", "trait_id", "phenotype_value")]) #tibble table with three columns
      test$phenotype_value <- as.numeric(test$phenotype_value) #make one column numeric 
      test <- distinct(test, plant_id, trait_id, .keep_all= TRUE) #remove duplicate rows based on multiple columns
      
      test <- test %>% 
        pivot_wider(
          names_from = trait_id,
          values_from = phenotype_value
        )
      
      # Using multiple conditions on DataFrame
      df <- test[,c("plant_id", b)]
      if (b == "FHBDON" || b == "FHBD3G" || b == "FHBZEA"){
        colnames(df)[2] <- "phenotype_value"
        df$select <- df$phenotype_value
        df$select <- ifelse(df$phenotype_value <= 0, "0", df$select) 
        df$select <- ifelse(df$phenotype_value > 0 & df$phenotype_value <= 0.5, "0.01-0.5", df$select) 
        df$select <- ifelse(df$phenotype_value > 0.5 & df$phenotype_value <= 1, "0.51-1", df$select) 
        df$select <- ifelse(df$phenotype_value > 1 & df$phenotype_value <= 5, "1.01-5", df$select) 
        df$select <- ifelse(df$phenotype_value > 5 & df$phenotype_value <= 10, "5.01-10", df$select) 
        df$select <- ifelse(df$phenotype_value > 10, ">10", df$select) 
        
        x = as.numeric(nrow(df)) #number of plants in population
        df.total = x - (sum(as.numeric(is.na(df$phenotype_value)))) #number of plants with data for trait
        df.one = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.two = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.three = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.four = ((as.numeric(nrow(df[df$select == "1.01-5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.five = ((as.numeric(nrow(df[df$select == "5.01-10",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.six = ((as.numeric(nrow(df[df$select == ">10",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.1 = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select)))))
        df.2 = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select)))))
        df.3 = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select)))))
        df.4 = ((as.numeric(nrow(df[df$select == "1.01-5",]))) - (sum(as.numeric(is.na(df$select)))))
        df.5 = ((as.numeric(nrow(df[df$select == "5.01-10",]))) - (sum(as.numeric(is.na(df$select)))))
        df.6 = ((as.numeric(nrow(df[df$select == ">10",]))) - (sum(as.numeric(is.na(df$select)))))
        
        sel = c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm")
        p = c(df.one, df.two, df.three, df.four, df.five, df.six)
        c = c(df.1, df.2, df.3, df.4, df.5, df.6)
        dff <- data.frame(sel,c,p)
        colnames(dff)[1] <- "toxin_level"
        colnames(dff)[2] <- "count"
        colnames(dff)[3] <- "percentage"
        dff$labels <- dff$percentage * 100
        dff$labels <- round(dff$labels,2)
        dff$labels <- paste0(dff$labels,"%")
        dff$count <- as.numeric(dff$count)
        dff$Plot <- paste0(a)
        
        dff$score <- factor(dff$toxin_level, levels=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))
        bdf <- rbind(bdf, dff)
        rm(temp,test,df,dff,x,df.one,df.two,df.three,df.four,df.five,df.six,df.1,df.2,df.3,df.4,df.5,df.6)
        bdf <- na.omit(bdf)
        assign(paste0(b,"table"), bdf, envir=parent.frame())
      } else if (b == "FHBDISIND"){
        colnames(df)[2] <- "phenotype_value"
        df$select <- df$phenotype_value
        df$select <- ifelse(df$phenotype_value <= 0, "0", df$select) 
        df$select <- ifelse(df$phenotype_value > 0 & df$phenotype_value <= 0.5, "0.01-0.5", df$select) 
        df$select <- ifelse(df$phenotype_value > 0.5 & df$phenotype_value <= 1, "0.51-1", df$select) 
        df$select <- ifelse(df$phenotype_value > 1 & df$phenotype_value <= 2, "1.01-2", df$select) 
        df$select <- ifelse(df$phenotype_value > 2 & df$phenotype_value <= 3, "2.01-3", df$select) 
        df$select <- ifelse(df$phenotype_value > 3 & df$phenotype_value <= 4, "3.01-4", df$select) 
        df$select <- ifelse(df$phenotype_value > 4 & df$phenotype_value <= 5, "4.01-5", df$select) 
        
        x = as.numeric(nrow(df)) #number of plants in population
        df.total = x - (sum(as.numeric(is.na(df$phenotype_value)))) #number of plants with data for trait
        df.one = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.two = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.three = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.four = ((as.numeric(nrow(df[df$select == "1.01-2",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.five = ((as.numeric(nrow(df[df$select == "2.01-3",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.six = ((as.numeric(nrow(df[df$select == "3.01-4",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.seven = ((as.numeric(nrow(df[df$select == "4.01-5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.1 = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select)))))
        df.2 = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select)))))
        df.3 = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select)))))
        df.4 = ((as.numeric(nrow(df[df$select == "1.01-2",]))) - (sum(as.numeric(is.na(df$select)))))
        df.5 = ((as.numeric(nrow(df[df$select == "2.01-3",]))) - (sum(as.numeric(is.na(df$select)))))
        df.6 = ((as.numeric(nrow(df[df$select == "3.01-4",]))) - (sum(as.numeric(is.na(df$select)))))
        df.7 = ((as.numeric(nrow(df[df$select == "4.01-5",]))) - (sum(as.numeric(is.na(df$select)))))
        
        sel = c("0", "0.01-0.5", "0.51-1", "1.01-2", "2.01-3", "3.01-4", "4.01-5")
        p = c(df.one, df.two, df.three, df.four, df.five, df.six, df.seven)
        c = c(df.1, df.2, df.3, df.4, df.5, df.6, df.7)
        dff <- data.frame(sel,c,p)
        colnames(dff)[1] <- "value"
        colnames(dff)[2] <- "count"
        colnames(dff)[3] <- "percentage"
        dff$labels <- dff$percentage * 100
        dff$labels <- round(dff$labels,2)
        dff$labels <- paste0(dff$labels,"%")
        dff$count <- as.numeric(dff$count)
        dff$Plot <- paste0(a)
        
        dff$score <- factor(dff$value, levels=c("0", "0.01-0.5", "0.51-1", "1.01-2", "2.01-3", "3.01-4", "4.01-5"))
        bdf <- rbind(bdf, dff)
        rm(df,dff,x,df.one,df.two,df.three,df.four,df.five,df.six,df.seven,df.1,df.2,df.3,df.4,df.5,df.6,df.7)
        bdf <- na.omit(bdf)
        assign(paste0(b,"table"), bdf, envir=parent.frame())  
      } else if (b == "FHBSPKINC") {
        colnames(df)[2] <- "phenotype_value"
        df$select <- df$phenotype_value
        df$select <- ifelse(df$phenotype_value <= 0, "0", df$select) 
        df$select <- ifelse(df$phenotype_value > 0 & df$phenotype_value <= 10, "1-10", df$select) 
        df$select <- ifelse(df$phenotype_value > 10 & df$phenotype_value <= 20, "11-20", df$select) 
        df$select <- ifelse(df$phenotype_value > 20 & df$phenotype_value <= 30, "21-30", df$select) 
        df$select <- ifelse(df$phenotype_value > 30 & df$phenotype_value <= 40, "31-40", df$select) 
        df$select <- ifelse(df$phenotype_value > 40 & df$phenotype_value <= 50, "41-50", df$select) 
        df$select <- ifelse(df$phenotype_value > 50 & df$phenotype_value <= 60, "51-60", df$select) 
        df$select <- ifelse(df$phenotype_value > 60 & df$phenotype_value <= 70, "61-70", df$select) 
        df$select <- ifelse(df$phenotype_value > 70 & df$phenotype_value <= 80, "71-80", df$select) 
        df$select <- ifelse(df$phenotype_value > 80 & df$phenotype_value <= 90, "81-90", df$select) 
        df$select <- ifelse(df$phenotype_value > 90 & df$phenotype_value <= 100, "91-100", df$select) 
        
        x = as.numeric(nrow(df)) #number of plants in population
        df.total = x - (sum(as.numeric(is.na(df$phenotype_value)))) #number of plants with data for trait
        df.one = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.two = ((as.numeric(nrow(df[df$select == "1-10",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.three = ((as.numeric(nrow(df[df$select == "11-20",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.four = ((as.numeric(nrow(df[df$select == "21-30",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.five = ((as.numeric(nrow(df[df$select == "31-40",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.six = ((as.numeric(nrow(df[df$select == "41-50",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.seven = ((as.numeric(nrow(df[df$select == "51-60",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.eight = ((as.numeric(nrow(df[df$select == "61-70",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.nine = ((as.numeric(nrow(df[df$select == "71-80",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.ten = ((as.numeric(nrow(df[df$select == "81-90",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.eleven = ((as.numeric(nrow(df[df$select == "91-100",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.1 = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select)))))
        df.2 = ((as.numeric(nrow(df[df$select == "1-10",]))) - (sum(as.numeric(is.na(df$select)))))
        df.3 = ((as.numeric(nrow(df[df$select == "11-20",]))) - (sum(as.numeric(is.na(df$select)))))
        df.4 = ((as.numeric(nrow(df[df$select == "21-30",]))) - (sum(as.numeric(is.na(df$select)))))
        df.5 = ((as.numeric(nrow(df[df$select == "31-40",]))) - (sum(as.numeric(is.na(df$select)))))
        df.6 = ((as.numeric(nrow(df[df$select == "41-50",]))) - (sum(as.numeric(is.na(df$select)))))
        df.7 = ((as.numeric(nrow(df[df$select == "51-60",]))) - (sum(as.numeric(is.na(df$select)))))
        df.8 = ((as.numeric(nrow(df[df$select == "61-70",]))) - (sum(as.numeric(is.na(df$select)))))
        df.9 = ((as.numeric(nrow(df[df$select == "71-80",]))) - (sum(as.numeric(is.na(df$select)))))
        df.10 = ((as.numeric(nrow(df[df$select == "81-90",]))) - (sum(as.numeric(is.na(df$select)))))
        df.11 = ((as.numeric(nrow(df[df$select == "91-100",]))) - (sum(as.numeric(is.na(df$select)))))
        
        sel = c("0", "1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100")
        p = c(df.one, df.two, df.three, df.four, df.five, df.six, df.seven, df.eight, df.nine, df.ten, df.eleven)
        c = c(df.1, df.2, df.3, df.4, df.5, df.6, df.7, df.8, df.9, df.10, df.11)
        dff <- data.frame(sel,c,p)
        colnames(dff)[1] <- "value"
        colnames(dff)[2] <- "count"
        colnames(dff)[3] <- "percentage"
        dff$labels <- dff$percentage * 100
        dff$labels <- round(dff$labels,2)
        dff$labels <- paste0(dff$labels,"%")
        dff$count <- as.numeric(dff$count)
        dff$Plot <- paste0(a)
        
        dff$score <- factor(dff$value, levels= c("0", "1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100"))
        bdf <- rbind(bdf, dff)
        rm(df,dff,x,df.one,df.two,df.three,df.four,df.five,df.six,df.seven,df.eight,df.nine,df.ten,df.eleven,df.1,df.2,df.3,df.4,df.5,df.6,df.7,df.8,df.9,df.10,df.11)
        bdf <- na.omit(bdf)
        assign(paste0(b,"table"), bdf, envir=parent.frame())  
      } else {print("if did not work")}
      rm(temp,test)
    }
    
    main_title = paste0(b)
    # FHBZEAtable, FHBDONtable, FHBD3Gtable, FHBDISINDtable, FHBSPKINCtable
    if (b == "FHBZEA") {
      dat <- FHBZEAtable
      plot <- ggplot(dat, aes(x = score, y = count, fill = score))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#F6AE2D", "#F26419", "#999933", "#999999"), labels= c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))+
        scale_x_discrete(limits=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm")) +
        geom_text(aes(label = labels), vjust = 2, colour = "black", size= 3)+
        facet_wrap(~Plot, scales = "fixed", ncol = 3)+
        labs(title= main_title, color=NULL) + 
        theme(strip.text.x = element_text(size = 10), legend.title=element_blank())
      plot_list[[b]] = plot
    } else if (b == "FHBDON") {
      dat <- FHBDONtable
      plot <- ggplot(dat, aes(x = score, y = count, fill = score))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#F6AE2D", "#F26419", "#999933", "#999999"), labels= c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))+
        scale_x_discrete(limits=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm")) +
        geom_text(aes(label = labels), vjust = 2, colour = "black", size= 3)+
        facet_wrap(~Plot, scales = "fixed", ncol = 3)+
        labs(title= main_title, color=NULL) + 
        theme(strip.text.x = element_text(size = 10), legend.title=element_blank())
      plot_list[[b]] = plot
    } else if (b == "FHBD3G") {
      dat <- FHBD3Gtable
      plot <- ggplot(dat, aes(x = score, y = count, fill = score))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#F6AE2D", "#F26419", "#999933", "#999999"), labels= c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))+
        scale_x_discrete(limits=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm")) +
        geom_text(aes(label = labels), vjust = 2, colour = "black", size= 3)+
        facet_wrap(~Plot, scales = "fixed", ncol = 3)+
        labs(title= main_title) + 
        theme(strip.text.x = element_text(size = 10), legend.title=element_blank())
      plot_list[[b]] = plot
    } else if (b == "FHBDISIND") {
      dat <- FHBDISINDtable
      plot <- ggplot(dat, aes(x = score, y = count, fill = score))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#44AA99", "#F6AE2D", "#F26419", "#999933", "#999999"), labels= c("0", "0.01-0.5", "0.51-1", "1.01-2", "2.01-3", "3.01-4", "4.01-5"))+
        scale_x_discrete(limits=c("0", "0.01-0.5", "0.51-1", "1.01-2", "2.01-3", "3.01-4", "4.01-5")) +
        geom_text(aes(label = labels), vjust = 2, colour = "black", size= 3)+
        facet_wrap(~Plot, scales = "fixed", ncol = 3)+
        labs(title= main_title) + 
        theme(strip.text.x = element_text(size = 10), legend.title=element_blank())
      plot_list[[b]] = plot
    } else if (b == "FHBSPKINC") {
      dat <- FHBSPKINCtable
      plot <- ggplot(dat, aes(x = score, y = count, fill = score))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#44AA99", "#F6AE2D", "#DDCC77", "#999933", "#F26419", "#CC6677", "#882255","#332288", "#999999"), labels= c("0", "1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100"))+
        scale_x_discrete(limits=c("0", "1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100")) +
        geom_text(aes(label = labels), vjust = 2, colour = "black", size= 3)+
        facet_wrap(~Plot, scales = "fixed", ncol = 3)+
        labs(title= main_title) + 
        theme(strip.text.x = element_text(size = 10), legend.title=element_blank())
      plot_list[[b]] = plot
    } else {print('help me, plots are not working')}
    
    trait_list <- names(plot_list)
    
    pdf(path, height = 15, width = 18)
    for (c in trait_list) {
      print(plot_list[[c]])
    }
    dev.off()
  }
}

#
#plot <- ggplot(dat, aes(x = score, y = count, fill = score))+
#geom_bar(stat = "identity")+
#scale_fill_manual(values = c("#55DDE0", "#33658A", "#F6AE2D", "#F26419", "#999933", "#999999"), labels= sel)+
#scale_x_discrete(limits=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm")) +
#geom_text(aes(label = labels), vjust = 2, colour = "black", size= 3)+
#facet_wrap(~Plot, scales = "fixed", ncol = 3)+
#labs(title= main_title) + 
#theme(strip.text.x = element_text(size = 10))

#plot_list[[b]] = plot




# Plots <- c("C07_SAL_2018", "C07_SAL_2019", "C07_SAL_2020", "C10_SAL_2021", "C10_OLA_2021", "C10_SAL_2022", "C10_OLA_2022")
# Traits <- c("FHBDON", "FHBD3G", "FHBDISIND", "FHBSPKINC")
# List <- list(C07_SAL_2018 = C7_2018_S, C07_SAL_2019 = C7_2019_S, C07_SAL_2020 =C7_2020_S, C10_SAL_2021 = C10_2021_S, C10_OLA_2021 = C10_2021_O, C10_SAL_2022 = C10_2022_S, C10_OLA_2022 = C10_2022_O)
# path <- 'data/Figures/PiePlot_FHB.pdf'
# pie_plot_fhb(Plot,Traits,List,path)
#Plots: vector of the names of each plot 
#Traits: vector of the names of traits to find and use from each table
#List: list containing each year's table, named with same names as in the plots vector, needs to be in long format and contain plant ids, trait ids, and trait values. there can be more trait ids than there is in Traits vector but those that are in both must be exact matches
#path: vector of output path, used to save plots as a pdf

pie_plot_fhb <- function(Plots,Traits,List,path){
  library(ggrepel)
  plot_list = list()
  for (b in Traits){
    b = paste0(b)
    bdf <- data.frame()
    for (a in Plots){
      temp = List[[a]] 
      test <- as_tibble(temp[c("plant_id", "trait_id", "phenotype_value")]) #tibble table with three columns
      test$phenotype_value <- as.numeric(test$phenotype_value) #make one column numeric 
      test <- distinct(test, plant_id, trait_id, .keep_all= TRUE) #remove duplicate rows based on multiple columns
      
      test <- test %>% 
        pivot_wider(
          names_from = trait_id,
          values_from = phenotype_value
        )
      
      # Using multiple conditions on DataFrame
      df <- test[,c("plant_id", b)]
      if (b == "FHBDON" || b == "FHBD3G" || b == "FHBZEA"){
        colnames(df)[2] <- "phenotype_value"
        df$select <- df$phenotype_value
        df$select <- ifelse(df$phenotype_value <= 0, "0", df$select) 
        df$select <- ifelse(df$phenotype_value > 0 & df$phenotype_value <= 0.5, "0.01-0.5", df$select) 
        df$select <- ifelse(df$phenotype_value > 0.5 & df$phenotype_value <= 1, "0.51-1", df$select) 
        df$select <- ifelse(df$phenotype_value > 1 & df$phenotype_value <= 5, "1.01-5", df$select) 
        df$select <- ifelse(df$phenotype_value > 5 & df$phenotype_value <= 10, "5.01-10", df$select) 
        df$select <- ifelse(df$phenotype_value > 10, ">10", df$select) 
        
        x = as.numeric(nrow(df)) #number of plants in population
        df.total = x - (sum(as.numeric(is.na(df$phenotype_value)))) #number of plants with data for trait
        df.one = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.two = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.three = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.four = ((as.numeric(nrow(df[df$select == "1.01-5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.five = ((as.numeric(nrow(df[df$select == "5.01-10",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.six = ((as.numeric(nrow(df[df$select == ">10",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.1 = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select)))))
        df.2 = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select)))))
        df.3 = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select)))))
        df.4 = ((as.numeric(nrow(df[df$select == "1.01-5",]))) - (sum(as.numeric(is.na(df$select)))))
        df.5 = ((as.numeric(nrow(df[df$select == "5.01-10",]))) - (sum(as.numeric(is.na(df$select)))))
        df.6 = ((as.numeric(nrow(df[df$select == ">10",]))) - (sum(as.numeric(is.na(df$select)))))
        
        sel = c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm")
        p = c(df.one, df.two, df.three, df.four, df.five, df.six)
        c = c(df.1, df.2, df.3, df.4, df.5, df.6)
        dff <- data.frame(sel,c,p)
        colnames(dff)[1] <- "toxin_level"
        colnames(dff)[2] <- "count"
        colnames(dff)[3] <- "percentage"
        dff$labels <- dff$percentage * 100
        dff$labels <- round(dff$labels,2)
        dff$labels <- paste0(dff$labels,"%")
        dff$count <- as.numeric(dff$count)
        dff$Plot <- paste0(a)
        
        dff$score <- factor(dff$toxin_level, levels=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))
        bdf <- rbind(bdf, dff)
        rm(temp,test,df,dff,x,df.one,df.two,df.three,df.four,df.five,df.six,df.1,df.2,df.3,df.4,df.5,df.6)
        bdf <- na.omit(bdf)
        assign(paste0(b,"table"), bdf, envir=parent.frame())
      } else if (b == "FHBDISIND"){
        colnames(df)[2] <- "phenotype_value"
        df$select <- df$phenotype_value
        df$select <- ifelse(df$phenotype_value <= 0, "0", df$select) 
        df$select <- ifelse(df$phenotype_value > 0 & df$phenotype_value <= 0.5, "0.01-0.5", df$select) 
        df$select <- ifelse(df$phenotype_value > 0.5 & df$phenotype_value <= 1, "0.51-1", df$select) 
        df$select <- ifelse(df$phenotype_value > 1 & df$phenotype_value <= 2, "1.01-2", df$select) 
        df$select <- ifelse(df$phenotype_value > 2 & df$phenotype_value <= 3, "2.01-3", df$select) 
        df$select <- ifelse(df$phenotype_value > 3 & df$phenotype_value <= 4, "3.01-4", df$select) 
        df$select <- ifelse(df$phenotype_value > 4 & df$phenotype_value <= 5, "4.01-5", df$select) 
        
        x = as.numeric(nrow(df)) #number of plants in population
        df.total = x - (sum(as.numeric(is.na(df$phenotype_value)))) #number of plants with data for trait
        df.one = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.two = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.three = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.four = ((as.numeric(nrow(df[df$select == "1.01-2",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.five = ((as.numeric(nrow(df[df$select == "2.01-3",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.six = ((as.numeric(nrow(df[df$select == "3.01-4",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.seven = ((as.numeric(nrow(df[df$select == "4.01-5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.1 = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select)))))
        df.2 = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select)))))
        df.3 = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select)))))
        df.4 = ((as.numeric(nrow(df[df$select == "1.01-2",]))) - (sum(as.numeric(is.na(df$select)))))
        df.5 = ((as.numeric(nrow(df[df$select == "2.01-3",]))) - (sum(as.numeric(is.na(df$select)))))
        df.6 = ((as.numeric(nrow(df[df$select == "3.01-4",]))) - (sum(as.numeric(is.na(df$select)))))
        df.7 = ((as.numeric(nrow(df[df$select == "4.01-5",]))) - (sum(as.numeric(is.na(df$select)))))
        
        sel = c("0", "0.01-0.5", "0.51-1", "1.01-2", "2.01-3", "3.01-4", "4.01-5")
        p = c(df.one, df.two, df.three, df.four, df.five, df.six, df.seven)
        c = c(df.1, df.2, df.3, df.4, df.5, df.6, df.7)
        dff <- data.frame(sel,c,p)
        colnames(dff)[1] <- "value"
        colnames(dff)[2] <- "count"
        colnames(dff)[3] <- "percentage"
        dff$labels <- dff$percentage * 100
        dff$labels <- round(dff$labels,2)
        dff$labels <- paste0(dff$labels,"%")
        dff$count <- as.numeric(dff$count)
        dff$Plot <- paste0(a)
        
        dff$score <- factor(dff$value, levels=c("0", "0.01-0.5", "0.51-1", "1.01-2", "2.01-3", "3.01-4", "4.01-5"))
        bdf <- rbind(bdf, dff)
        rm(df,dff,x,df.one,df.two,df.three,df.four,df.five,df.six,df.seven,df.1,df.2,df.3,df.4,df.5,df.6,df.7)
        bdf <- na.omit(bdf)
        assign(paste0(b,"table"), bdf, envir=parent.frame())  
      } else if (b == "FHBSPKINC") {
        colnames(df)[2] <- "phenotype_value"
        df$select <- df$phenotype_value
        df$select <- ifelse(df$phenotype_value <= 0, "0", df$select) 
        df$select <- ifelse(df$phenotype_value > 0 & df$phenotype_value <= 10, "1-10", df$select) 
        df$select <- ifelse(df$phenotype_value > 10 & df$phenotype_value <= 20, "11-20", df$select) 
        df$select <- ifelse(df$phenotype_value > 20 & df$phenotype_value <= 30, "21-30", df$select) 
        df$select <- ifelse(df$phenotype_value > 30 & df$phenotype_value <= 40, "31-40", df$select) 
        df$select <- ifelse(df$phenotype_value > 40 & df$phenotype_value <= 50, "41-50", df$select) 
        df$select <- ifelse(df$phenotype_value > 50 & df$phenotype_value <= 60, "51-60", df$select) 
        df$select <- ifelse(df$phenotype_value > 60 & df$phenotype_value <= 70, "61-70", df$select) 
        df$select <- ifelse(df$phenotype_value > 70 & df$phenotype_value <= 80, "71-80", df$select) 
        df$select <- ifelse(df$phenotype_value > 80 & df$phenotype_value <= 90, "81-90", df$select) 
        df$select <- ifelse(df$phenotype_value > 90 & df$phenotype_value <= 100, "91-100", df$select) 
        
        x = as.numeric(nrow(df)) #number of plants in population
        df.total = x - (sum(as.numeric(is.na(df$phenotype_value)))) #number of plants with data for trait
        df.one = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.two = ((as.numeric(nrow(df[df$select == "1-10",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.three = ((as.numeric(nrow(df[df$select == "11-20",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.four = ((as.numeric(nrow(df[df$select == "21-30",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.five = ((as.numeric(nrow(df[df$select == "31-40",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.six = ((as.numeric(nrow(df[df$select == "41-50",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.seven = ((as.numeric(nrow(df[df$select == "51-60",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.eight = ((as.numeric(nrow(df[df$select == "61-70",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.nine = ((as.numeric(nrow(df[df$select == "71-80",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.ten = ((as.numeric(nrow(df[df$select == "81-90",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.eleven = ((as.numeric(nrow(df[df$select == "91-100",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.1 = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select)))))
        df.2 = ((as.numeric(nrow(df[df$select == "1-10",]))) - (sum(as.numeric(is.na(df$select)))))
        df.3 = ((as.numeric(nrow(df[df$select == "11-20",]))) - (sum(as.numeric(is.na(df$select)))))
        df.4 = ((as.numeric(nrow(df[df$select == "21-30",]))) - (sum(as.numeric(is.na(df$select)))))
        df.5 = ((as.numeric(nrow(df[df$select == "31-40",]))) - (sum(as.numeric(is.na(df$select)))))
        df.6 = ((as.numeric(nrow(df[df$select == "41-50",]))) - (sum(as.numeric(is.na(df$select)))))
        df.7 = ((as.numeric(nrow(df[df$select == "51-60",]))) - (sum(as.numeric(is.na(df$select)))))
        df.8 = ((as.numeric(nrow(df[df$select == "61-70",]))) - (sum(as.numeric(is.na(df$select)))))
        df.9 = ((as.numeric(nrow(df[df$select == "71-80",]))) - (sum(as.numeric(is.na(df$select)))))
        df.10 = ((as.numeric(nrow(df[df$select == "81-90",]))) - (sum(as.numeric(is.na(df$select)))))
        df.11 = ((as.numeric(nrow(df[df$select == "91-100",]))) - (sum(as.numeric(is.na(df$select)))))
        
        sel = c("0", "1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100")
        p = c(df.one, df.two, df.three, df.four, df.five, df.six, df.seven, df.eight, df.nine, df.ten, df.eleven)
        c = c(df.1, df.2, df.3, df.4, df.5, df.6, df.7, df.8, df.9, df.10, df.11)
        dff <- data.frame(sel,c,p)
        colnames(dff)[1] <- "value"
        colnames(dff)[2] <- "count"
        colnames(dff)[3] <- "percentage"
        dff$labels <- dff$percentage * 100
        dff$labels <- round(dff$labels,2)
        dff$labels <- paste0(dff$labels,"%")
        dff$count <- as.numeric(dff$count)
        dff$Plot <- paste0(a)
        
        dff$score <- factor(dff$value, levels= c("0", "1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100"))
        bdf <- rbind(bdf, dff)
        rm(df,dff,x,df.one,df.two,df.three,df.four,df.five,df.six,df.seven,df.eight,df.nine,df.ten,df.eleven,df.1,df.2,df.3,df.4,df.5,df.6,df.7,df.8,df.9,df.10,df.11)
        bdf <- na.omit(bdf)
        assign(paste0(b,"table"), bdf, envir=parent.frame())  
      } else {print("it did not work")}
      rm(temp,test)
    }
    
    main_title = paste0(b)
    # FHBZEAtable, FHBDONtable, FHBD3Gtable, FHBDISINDtable, FHBSPKINCtable
    if (b == "FHBZEA") {
      dat <- FHBZEAtable
      dat$labels <- ifelse(dat$labels == "0%", NA, dat$labels)
      dat$scorer <- ifelse(dat$count == 0, NA, dat$score)
      dat$scorer <- ifelse(dat$scorer == 1, "0 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 2, "0.01-0.5 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 3, "0.51-1 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 4, "1.01-5 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 5, "5.01-10 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 6, ">10 ppm", dat$scorer)
      dat$scorer <- factor(dat$scorer, levels=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))
      plot <- ggplot(dat, aes(x = "", y = percentage, fill = score)) +
        geom_col(stat="identity", width=1, color = "black") +
        coord_polar("y", start=0, direction=-1)+
        geom_text(aes(x = 1.35, label = labels), color="white", size= 3, angle = 50,position = position_stack(vjust = .5)) +
        labs(x = NULL, y = NULL, fill = NULL)+
        theme_void()+
        theme(axis.ticks=element_blank(),
              axis.title=element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              panel.grid  = element_blank())+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#999933", "#F6AE2D","#F26419", "#999999"))+
        facet_wrap(~Plot, ncol = 3)+
        geom_text(aes(x = 1.7, label = scorer), color="black", size= 3, angle = 60, position = position_stack(vjust = .5)) +
        labs(title= main_title, color=NULL)
      plot_list[[b]] = plot
    } else if (b == "FHBDON") {
      dat <- FHBDONtable
      dat$labels <- ifelse(dat$labels == "0%", NA, dat$labels)
      dat$scorer <- ifelse(dat$count == 0, NA, dat$score)
      dat$scorer <- ifelse(dat$scorer == 1, "0 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 2, "0.01-0.5 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 3, "0.51-1 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 4, "1.01-5 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 5, "5.01-10 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 6, ">10 ppm", dat$scorer)
      dat$scorer <- factor(dat$scorer, levels=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))
      plot <- ggplot(dat, aes(x = "", y = percentage, fill = score)) +
        geom_col(stat="identity", width=1, color = "black") +
        coord_polar("y", start=0, direction=-1)+
        geom_text(aes(x = 1.35, label = labels), color="white", size= 3, angle = 50,position = position_stack(vjust = .5)) +
        labs(x = NULL, y = NULL, fill = NULL)+
        theme_void()+
        theme(axis.ticks=element_blank(),
              axis.title=element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              panel.grid  = element_blank())+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#999933", "#F6AE2D","#F26419", "#999999"))+
        facet_wrap(~Plot, ncol = 3)+
        geom_text(aes(x = 1.7, label = scorer), color="black", size= 3, angle = 60, position = position_stack(vjust = .5)) +
        labs(title= main_title, color=NULL)
      plot_list[[b]] = plot
    } else if (b == "FHBD3G") {
      dat <- FHBD3Gtable
      dat$labels <- ifelse(dat$labels == "0%", NA, dat$labels)
      dat$scorer <- ifelse(dat$count == 0, NA, dat$score)
      dat$scorer <- ifelse(dat$scorer == 1, "0 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 2, "0.01-0.5 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 3, "0.51-1 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 4, "1.01-5 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 5, "5.01-10 ppm", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 6, ">10 ppm", dat$scorer)
      dat$scorer <- factor(dat$scorer, levels=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))
      plot <- ggplot(dat, aes(x = "", y = percentage, fill = score)) +
        geom_col(stat="identity", width=1, color = "black") +
        coord_polar("y", start=0, direction=-1)+
        geom_text(aes(x = 1.3, label = labels), color="white", size= 3, angle = 45, position = position_stack(vjust = 0.5)) +
        labs(x = NULL, y = NULL, fill = NULL)+
        theme_void()+
        theme(axis.ticks=element_blank(),
              axis.title=element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              panel.grid  = element_blank())+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#999933", "#F6AE2D", "#F26419", "#999999"))+
        facet_wrap(~Plot, ncol = 3)+
        geom_text(aes(x = 1.7, label = scorer), color="black", size= 3, angle = 60, position = position_stack(vjust = 0.5)) +
        labs(title= main_title, color=NULL)
      plot_list[[b]] = plot
    } else if (b == "FHBDISIND") {
      dat <- FHBDISINDtable
      dat$labels <- ifelse(dat$labels == "0%", NA, dat$labels)
      dat$scorer <- ifelse(dat$count == 0, NA, dat$score)
      dat$scorer <- ifelse(dat$scorer == 1, "0", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 2, "0.01-0.5", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 3, "0.51-1", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 4, "1.01-2", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 5, "2.01-3", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 6, "3.01-4", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 7, "4.01-5", dat$scorer)
      dat$scorer <- factor(dat$scorer, levels=c("0", "0.01-0.5", "0.51-1", "1.01-2", "2.01-3", "3.01-4", "4.01-5"))
      plot <- ggplot(dat, aes(x = "", y = percentage, fill = score)) +
        geom_col(stat="identity", width=1, color = "black") +
        coord_polar("y", start=0, direction=-1)+
        geom_text(aes(x = 1.3, label = labels), color="white", size= 3, angle = 45, position = position_stack(vjust = .5)) +
        labs(x = NULL, y = NULL, fill = NULL)+
        theme_void()+
        theme(axis.ticks=element_blank(),
              axis.title=element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              panel.grid  = element_blank())+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#44AA99", "#999933", "#F6AE2D", "#F26419", "#999999"))+
        facet_wrap(~Plot, ncol = 3)+
        geom_text(aes(x = 1.7, label = scorer), color="black", size= 3, angle = 60, position = position_stack(vjust = .5)) +
        labs(title= main_title, color=NULL)
      plot_list[[b]] = plot
    } else if (b == "FHBSPKINC") {
      dat <- FHBSPKINCtable
      dat$labels <- ifelse(dat$labels == "0%", NA, dat$labels)
      dat$scorer <- ifelse(dat$count == 0, NA, dat$score)
      dat$scorer <- ifelse(dat$scorer == 1, "0", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 2, "1-10", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 3, "11-20", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 4, "21-30", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 5, "31-40", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 6, "41-50", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 7, "51-60", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 8, "61-70", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 9, "71-80", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 10, "81-90", dat$scorer)
      dat$scorer <- ifelse(dat$scorer == 11, "91-100", dat$scorer)
      dat$scorer <- factor(dat$scorer, levels=c("0", "1-10", "11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100"))
      plot <- ggplot(dat, aes(x = "", y = percentage, fill = score)) +
        geom_col(stat="identity", width=1, color = "black") +
        coord_polar("y", start=0, direction=-1)+
        geom_text(aes(x = 1.3, label = labels), color="white", size= 3, angle = 45, position = position_stack(vjust = .5)) +
        labs(x = NULL, y = NULL, fill = NULL)+
        theme_void()+
        theme(axis.ticks=element_blank(),
              axis.title=element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              panel.grid  = element_blank())+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#44AA99", "#999933", "#DDCC77", "#F6AE2D","#F26419", "#CC6677", "#882255","#332288", "#999999"))+
        facet_wrap(~Plot, ncol = 3)+
        geom_text(aes(x = 1.7, label = scorer), color="black", size= 3, angle = 60, position = position_stack(vjust = .5)) +
        labs(title= main_title, color=NULL)
      plot_list[[b]] = plot
    } else {print('help me, plots are not working')}
    
    trait_list <- names(plot_list)
    
    pdf(path, height = 15, width = 18)
    for (c in trait_list) {
      print(plot_list[[c]])
    }
    dev.off()
  }
}


#pie_plot_fhb(Plots,Traits,List,path)


#ggplot(db, aes(x = "", y = percentage, fill = toxin_level)) +
#geom_col(stat="identity", width=1, color = "black") +
#coord_polar("y", start=0)+
#geom_text(aes(x = 1.3, label = labels), color="white", size= 3, position = position_stack(vjust = .5)) +
#labs(x = NULL, y = NULL, fill = NULL)+
#theme_void()+
#theme(axis.ticks=element_blank(),
#axis.title=element_blank(),
#axis.text.y = element_blank(),
#axis.text.x = element_blank(),
#panel.grid  = element_blank())+
#scale_fill_manual(values = c("#55DDE0", "#33658A", "#F6AE2D", "#F26419", "#999933", "#999999"))+
#facet_wrap(~Plot, ncol = 2)+
#labs(title= main_title, color=NULL)


# Plots <- c("C07_SAL_2018", "C07_SAL_2019", "C07_SAL_2020", "C10_SAL_2021", "C10_OLA_2021", "C10_SAL_2022", "C10_OLA_2022")
# Traits <- c("FHBDON", "FHBD3G", "FHBDISIND", "FHBSPKINC")
# List <- list(C07_SAL_2018 = C7_2018_S, C07_SAL_2019 = C7_2019_S, C07_SAL_2020 =C7_2020_S, C10_SAL_2021 = C10_2021_S, C10_OLA_2021 = C10_2021_O, C10_SAL_2022 = C10_2022_S, C10_OLA_2022 = C10_2022_O)
# path <- 'data/Figures/BarPlot_FHB.pdf'
# bar_plot_fhb(Plot,Traits,List,path)
#Plots: vector of the names of each plot 
#Traits: vector of the names of traits to find and use from each table
#List: list containing each year's table, named with same names as in the plots vector, needs to be in long format and contain plant ids, trait ids, and trait values. there can be more trait ids than there is in Traits vector but those that are in both must be exact matches
#path: vector of output path, used to save plots as a pdf

#if adding traits, will have to add in other ifelse statements to customize the scale
####################
bar_plot_fhb2 <- function(Plots,Traits,List,List3,path){
  plot_list = list()
  for (b in Traits){
    b = paste0(b)
    bdf <- data.frame()
    for (a in Plots){
      temp = List[[a]] 
      test <- as_tibble(temp[c("plant_id", "trait_id", "phenotype_value")]) #tibble table with three columns
      test$phenotype_value <- as.numeric(test$phenotype_value) #make one column numeric 
      test <- distinct(test, plant_id, trait_id, .keep_all= TRUE) #remove duplicate rows based on multiple columns
      
      test <- test %>% 
        pivot_wider(
          names_from = trait_id,
          values_from = phenotype_value
        )
      
      # Using multiple conditions on DataFrame
      df <- test[,c("plant_id", b)]
      if (b == "FHBDON" || b == "FHBD3G" || b == "FHBZEA"){
        colnames(df)[2] <- "phenotype_value"
        df$select <- df$phenotype_value
        df$select <- ifelse(df$phenotype_value <= 0, "0", df$select) 
        df$select <- ifelse(df$phenotype_value > 0 & df$phenotype_value <= 0.5, "0.01-0.5", df$select) 
        df$select <- ifelse(df$phenotype_value > 0.5 & df$phenotype_value <= 1, "0.51-1", df$select) 
        df$select <- ifelse(df$phenotype_value > 1 & df$phenotype_value <= 5, "1.01-5", df$select) 
        df$select <- ifelse(df$phenotype_value > 5 & df$phenotype_value <= 10, "5.01-10", df$select) 
        df$select <- ifelse(df$phenotype_value > 10, ">10", df$select) 
        
        x = as.numeric(nrow(df)) #number of plants in population
        df.total = x - (sum(as.numeric(is.na(df$phenotype_value)))) #number of plants with data for trait
        df.one = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.two = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.three = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.four = ((as.numeric(nrow(df[df$select == "1.01-5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.five = ((as.numeric(nrow(df[df$select == "5.01-10",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.six = ((as.numeric(nrow(df[df$select == ">10",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.1 = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select)))))
        df.2 = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select)))))
        df.3 = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select)))))
        df.4 = ((as.numeric(nrow(df[df$select == "1.01-5",]))) - (sum(as.numeric(is.na(df$select)))))
        df.5 = ((as.numeric(nrow(df[df$select == "5.01-10",]))) - (sum(as.numeric(is.na(df$select)))))
        df.6 = ((as.numeric(nrow(df[df$select == ">10",]))) - (sum(as.numeric(is.na(df$select)))))
        
        sel = c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm")
        p = c(df.one, df.two, df.three, df.four, df.five, df.six)
        c = c(df.1, df.2, df.3, df.4, df.5, df.6)
        dff <- data.frame(sel,c,p)
        colnames(dff)[1] <- "toxin_level"
        colnames(dff)[2] <- "count"
        colnames(dff)[3] <- "percentage"
        dff$labels <- dff$percentage * 100
        dff$labels <- round(dff$labels,2)
        dff$labels <- paste0(dff$labels,"%")
        dff$count <- as.numeric(dff$count)
        dff$Plot <- paste0(a)
        
        dff$score <- factor(dff$toxin_level, levels=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))
        bdf <- rbind(bdf, dff)
        rm(temp,test,df,dff,x,df.one,df.two,df.three,df.four,df.five,df.six,df.1,df.2,df.3,df.4,df.5,df.6)
        bdf <- na.omit(bdf)
        assign(paste0(b,"table"), bdf, envir=parent.frame())
      } else if (b == "FHBDISIND"){
        colnames(df)[2] <- "phenotype_value"
        df$select <- df$phenotype_value
        df$select <- ifelse(df$phenotype_value <= 0, "0", df$select) 
        df$select <- ifelse(df$phenotype_value > 0 & df$phenotype_value <= 0.5, "0.01-0.5", df$select) 
        df$select <- ifelse(df$phenotype_value > 0.5 & df$phenotype_value <= 1, "0.51-1", df$select) 
        df$select <- ifelse(df$phenotype_value > 1 & df$phenotype_value <= 2, "1.01-2", df$select) 
        df$select <- ifelse(df$phenotype_value > 2 & df$phenotype_value <= 3, "2.01-3", df$select) 
        df$select <- ifelse(df$phenotype_value > 3 & df$phenotype_value <= 4, "3.01-4", df$select) 
        df$select <- ifelse(df$phenotype_value > 4 & df$phenotype_value <= 5, "4.01-5", df$select) 
        
        x = as.numeric(nrow(df)) #number of plants in population
        df.total = x - (sum(as.numeric(is.na(df$phenotype_value)))) #number of plants with data for trait
        df.one = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.two = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.three = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.four = ((as.numeric(nrow(df[df$select == "1.01-2",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.five = ((as.numeric(nrow(df[df$select == "2.01-3",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.six = ((as.numeric(nrow(df[df$select == "3.01-4",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.seven = ((as.numeric(nrow(df[df$select == "4.01-5",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.1 = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select)))))
        df.2 = ((as.numeric(nrow(df[df$select == "0.01-0.5",]))) - (sum(as.numeric(is.na(df$select)))))
        df.3 = ((as.numeric(nrow(df[df$select == "0.51-1",]))) - (sum(as.numeric(is.na(df$select)))))
        df.4 = ((as.numeric(nrow(df[df$select == "1.01-2",]))) - (sum(as.numeric(is.na(df$select)))))
        df.5 = ((as.numeric(nrow(df[df$select == "2.01-3",]))) - (sum(as.numeric(is.na(df$select)))))
        df.6 = ((as.numeric(nrow(df[df$select == "3.01-4",]))) - (sum(as.numeric(is.na(df$select)))))
        df.7 = ((as.numeric(nrow(df[df$select == "4.01-5",]))) - (sum(as.numeric(is.na(df$select)))))
        
        sel = c("0", "0.01-0.5", "0.51-1", "1.01-2", "2.01-3", "3.01-4", "4.01-5")
        p = c(df.one, df.two, df.three, df.four, df.five, df.six, df.seven)
        c = c(df.1, df.2, df.3, df.4, df.5, df.6, df.7)
        dff <- data.frame(sel,c,p)
        colnames(dff)[1] <- "value"
        colnames(dff)[2] <- "count"
        colnames(dff)[3] <- "percentage"
        dff$labels <- dff$percentage * 100
        dff$labels <- round(dff$labels,2)
        dff$labels <- paste0(dff$labels,"%")
        dff$count <- as.numeric(dff$count)
        dff$Plot <- paste0(a)
        
        dff$score <- factor(dff$value, levels=c("0", "0.01-0.5", "0.51-1", "1.01-2", "2.01-3", "3.01-4", "4.01-5"))
        bdf <- rbind(bdf, dff)
        rm(df,dff,x,df.one,df.two,df.three,df.four,df.five,df.six,df.seven,df.1,df.2,df.3,df.4,df.5,df.6,df.7)
        bdf <- na.omit(bdf)
        assign(paste0(b,"table"), bdf, envir=parent.frame())  
      } else if (b == "FHBSPKINC") {
        colnames(df)[2] <- "phenotype_value"
        df$select <- df$phenotype_value
        df$select <- ifelse(df$phenotype_value <= 0, "0", df$select) 
        df$select <- ifelse(df$phenotype_value > 0 & df$phenotype_value <= 10, "1-10", df$select) 
        df$select <- ifelse(df$phenotype_value > 10 & df$phenotype_value <= 20, "11-20", df$select) 
        df$select <- ifelse(df$phenotype_value > 20 & df$phenotype_value <= 30, "21-30", df$select) 
        df$select <- ifelse(df$phenotype_value > 30 & df$phenotype_value <= 40, "31-40", df$select) 
        df$select <- ifelse(df$phenotype_value > 40 & df$phenotype_value <= 50, "41-50", df$select) 
        df$select <- ifelse(df$phenotype_value > 50 & df$phenotype_value <= 60, "51-60", df$select) 
        df$select <- ifelse(df$phenotype_value > 60 & df$phenotype_value <= 70, "61-70", df$select) 
        df$select <- ifelse(df$phenotype_value > 70 & df$phenotype_value <= 80, "71-80", df$select) 
        df$select <- ifelse(df$phenotype_value > 80 & df$phenotype_value <= 90, "81-90", df$select) 
        df$select <- ifelse(df$phenotype_value > 90 & df$phenotype_value <= 100, "91-100", df$select) 
        
        x = as.numeric(nrow(df)) #number of plants in population
        df.total = x - (sum(as.numeric(is.na(df$phenotype_value)))) #number of plants with data for trait
        df.one = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.two = ((as.numeric(nrow(df[df$select == "1-10",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.three = ((as.numeric(nrow(df[df$select == "11-20",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.four = ((as.numeric(nrow(df[df$select == "21-30",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.five = ((as.numeric(nrow(df[df$select == "31-40",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.six = ((as.numeric(nrow(df[df$select == "41-50",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.seven = ((as.numeric(nrow(df[df$select == "51-60",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.eight = ((as.numeric(nrow(df[df$select == "61-70",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.nine = ((as.numeric(nrow(df[df$select == "71-80",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.ten = ((as.numeric(nrow(df[df$select == "81-90",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.eleven = ((as.numeric(nrow(df[df$select == "91-100",]))) - (sum(as.numeric(is.na(df$select))))) / df.total
        df.1 = ((as.numeric(nrow(df[df$select == "0",]))) - (sum(as.numeric(is.na(df$select)))))
        df.2 = ((as.numeric(nrow(df[df$select == "1-10",]))) - (sum(as.numeric(is.na(df$select)))))
        df.3 = ((as.numeric(nrow(df[df$select == "11-20",]))) - (sum(as.numeric(is.na(df$select)))))
        df.4 = ((as.numeric(nrow(df[df$select == "21-30",]))) - (sum(as.numeric(is.na(df$select)))))
        df.5 = ((as.numeric(nrow(df[df$select == "31-40",]))) - (sum(as.numeric(is.na(df$select)))))
        df.6 = ((as.numeric(nrow(df[df$select == "41-50",]))) - (sum(as.numeric(is.na(df$select)))))
        df.7 = ((as.numeric(nrow(df[df$select == "51-60",]))) - (sum(as.numeric(is.na(df$select)))))
        df.8 = ((as.numeric(nrow(df[df$select == "61-70",]))) - (sum(as.numeric(is.na(df$select)))))
        df.9 = ((as.numeric(nrow(df[df$select == "71-80",]))) - (sum(as.numeric(is.na(df$select)))))
        df.10 = ((as.numeric(nrow(df[df$select == "81-90",]))) - (sum(as.numeric(is.na(df$select)))))
        df.11 = ((as.numeric(nrow(df[df$select == "91-100",]))) - (sum(as.numeric(is.na(df$select)))))
        
        sel = c("0", "1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100")
        p = c(df.one, df.two, df.three, df.four, df.five, df.six, df.seven, df.eight, df.nine, df.ten, df.eleven)
        c = c(df.1, df.2, df.3, df.4, df.5, df.6, df.7, df.8, df.9, df.10, df.11)
        dff <- data.frame(sel,c,p)
        colnames(dff)[1] <- "value"
        colnames(dff)[2] <- "count"
        colnames(dff)[3] <- "percentage"
        dff$labels <- dff$percentage * 100
        dff$labels <- round(dff$labels,2)
        dff$labels <- paste0(dff$labels,"%")
        dff$count <- as.numeric(dff$count)
        dff$Plot <- paste0(a)
        
        dff$score <- factor(dff$value, levels= c("0", "1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100"))
        bdf <- rbind(bdf, dff)
        rm(df,dff,x,df.one,df.two,df.three,df.four,df.five,df.six,df.seven,df.eight,df.nine,df.ten,df.eleven,df.1,df.2,df.3,df.4,df.5,df.6,df.7,df.8,df.9,df.10,df.11)
        bdf <- na.omit(bdf)
        assign(paste0(b,"table"), bdf, envir=parent.frame())  
      } else {print("if did not work")}
      rm(temp,test)
    }
    
    main_title = paste0(b)
    # FHBZEAtable, FHBDONtable, FHBD3Gtable, FHBDISINDtable, FHBSPKINCtable
    if (b == "FHBZEA") {
      dat <- FHBZEAtable
      plot <- ggplot(dat, aes(x = score, y = count, fill = score))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#F6AE2D", "#F26419", "#999933", "#999999"), labels= c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))+
        scale_x_discrete(limits=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm")) +
        geom_text(aes(label = labels), vjust = 2, colour = "black", size= 3)+
        facet_wrap(~Plot, scales = "fixed", ncol = 3)+
        labs(title= main_title, color=NULL) + 
        theme(strip.text.x = element_text(size = 10), legend.title=element_blank())+
        theme(panel.background = element_blank())+ 
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.text = element_text(size=1))
      plot_list[[b]] = plot
    } else if (b == "FHBDON") {
      dat <- FHBDONtable
      A = List[["C07_SAL_2018"]]
      B = List[["C07_SAL_2019"]]
      C = List[["C07_SAL_2020"]]
      D = List[["C07_SAL_2021"]]
      E = List[["C07_OLA_2021"]]
      G = List[["C07_SAL_2022"]]
      H = List[["C07_OLA_2022"]]
      des <- "
        ABC
        GD#
        HE#
      "
      plot <- ggplot(dat, aes(x = score, y = count, fill = score))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#F6AE2D", "#F26419", "#999933", "#999999"), labels= c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))+
        scale_x_discrete(limits=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm")) +
        geom_text(aes(label = labels), vjust = 2, colour = "black", size= 3)+
        facet_manual(~Plot, scales = "fixed", design = des)+
        labs(title= main_title, color=NULL) + 
        theme(strip.text.x = element_text(size = 10), legend.title=element_blank())+
        theme(panel.background = element_blank())+ 
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.text = element_text(size=10))
      plot_list[[b]] = plot
    } else if (b == "FHBD3G") {
      dat <- FHBD3Gtable
      plot <- ggplot(dat, aes(x = score, y = count, fill = score))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#F6AE2D", "#F26419", "#999933", "#999999"), labels= c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm"))+
        scale_x_discrete(limits=c("0 ppm", "0.01-0.5 ppm", "0.51-1 ppm", "1.01-5 ppm", "5.01-10 ppm", ">10 ppm")) +
        geom_text(aes(label = labels), vjust = 2, colour = "black", size= 3)+
        facet_wrap(~Plot, scales = "fixed", ncol = 3)+
        labs(title= main_title) + 
        theme(strip.text.x = element_text(size = 10), legend.title=element_blank())+
        theme(panel.background = element_blank())+ 
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.text = element_text(size=10))
      plot_list[[b]] = plot
    } else if (b == "FHBDISIND") {
      dat <- FHBDISINDtable
      A = List[["C07_SAL_2018"]]
      B = List[["C07_SAL_2019"]]
      C = List[["C07_SAL_2020"]]
      D = List[["C07_SAL_2021"]]
      E = List[["C07_OLA_2021"]]
      G = List[["C07_SAL_2022"]]
      H = List[["C07_OLA_2022"]]
      des <- "
        ABC
        GD#
        HE#
      "
      plot <- ggplot(dat, aes(x = score, y = count, fill = score))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#44AA99", "#F6AE2D", "#F26419", "#999933", "#999999"), labels= c("0", "0.01-0.5", "0.51-1", "1.01-2", "2.01-3", "3.01-4", "4.01-5"))+
        scale_x_discrete(limits=c("0", "0.01-0.5", "0.51-1", "1.01-2", "2.01-3", "3.01-4", "4.01-5")) +
        geom_text(aes(label = labels), vjust = 2, colour = "black", size= 3)+
        facet_manual(~Plot, scales = "fixed", design=des)+
        geom_vline(data = List3$FHBDISIND, aes(xintercept=select), color= 'black')+
        labs(title= main_title) + 
        theme(strip.text.x = element_text(size = 10), legend.title=element_blank())+
        theme(panel.background = element_blank())+ 
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.text = element_text(size=10))
      plot_list[[b]] = plot
    } else if (b == "FHBSPKINC") {
      dat <- FHBSPKINCtable
      A = List[["C07_SAL_2018"]]
      B = List[["C07_SAL_2019"]]
      C = List[["C07_SAL_2020"]]
      D = List[["C07_SAL_2021"]]
      E = List[["C07_OLA_2021"]]
      G = List[["C07_SAL_2022"]]
      H = List[["C07_OLA_2022"]]
      des <- "
        ABC
        GD#
        HE#
      "
      plot <- ggplot(dat, aes(x = score, y = count, fill = score))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = c("#55DDE0", "#33658A", "#44AA99", "#F6AE2D", "#DDCC77", "#999933", "#F26419", "#CC6677", "#882255","#332288", "#999999"), labels= c("0", "1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100"))+
        scale_x_discrete(limits=c("0", "1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100")) +
        geom_text(aes(label = labels), vjust = 2, colour = "black", size= 3)+
        facet_manual(~Plot, scales = "fixed", design=des)+
        geom_vline(data = List3$FHBSPKINC, aes(xintercept=select), color= 'black')+
        labs(title= main_title) + 
        theme(strip.text.x = element_text(size = 10), legend.title=element_blank())+
        theme(panel.background = element_blank())+ 
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.text = element_text(size=10))
      plot_list[[b]] = plot
    } else {print('help me, plots are not working')}
    
    trait_list <- names(plot_list)
    
    pdf(path, height = 15, width = 18)
    for (c in trait_list) {
      print(plot_list[[c]])
    }
    dev.off()
  }
}


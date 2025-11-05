data <- read.csv('~/Desktop/output.csv', header = T)

###
# Example p-values and effect sizes for SNPs (replace with your actual data)
p_values <- c(0.01, 0.005, 0.02, 0.1)  # Example p-values
effect_sizes <- c(0.1, 0.2, 0.15, -0.08)  # Example effect sizes

# Example total variance of the phenotype
total_variance_phenotype <- 100

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - p_values) * p_values) * effect_sizes^2

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype

# Print PVE for each SNP
print(pve_snp)

# Assuming you have vectors of p-values, effect sizes, and MAFs for each SNP
p_values <- c(0.001, 0.002, 0.003)  # Example p-values
effect_sizes <- c(0.1, 0.2, 0.3)  # Example effect sizes
maf <- c(0.2, 0.3, 0.4)  # Example MAFs

# Calculate PVE for each SNP
pve <- 2 * (1 - maf) * maf * effect_sizes^2

# Print the PVE values
print(pve)
  
  
###
# Example phenotype data (replace with your actual phenotype data)
##phenotype <- c(20, 25, 30, 35, 40)  # Example phenotype measurements
# My phenotypes
#setwd('/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS')
#myY <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C07_SAL_2018_pheno.txt", head=TRUE)

# Calculate the total variance of the phenotype
##total_variance <- var(phenotype)
#total_variance <- var(myY$FHBDON, na.rm = T)
###


#Replace p_values and effect_sizes with your actual p-values and effect sizes. 
#Also, adjust total_variance to reflect the total variance of your phenotype.

#This code will calculate PVE using the formula you provided, 
#where alpha^2 is the sum of (2(1 - p)p * beta^2) over all significant SNPs 
#(where p is the p-value and beta is the effect size). 
#This formula accounts for the contribution of each SNP to the variance explained by the genetic variants.

#p_values <- data$P.value
#effect_sizes <- data$effect

#alpha_squared <- (2 * (1 - p_values) * p_values) * effect_sizes^2



#####
### FHBDON 2018
# p-values and effect sizes

subset_snps <- subset(data, Trait == "FarmCPU.FHBDON" & Plot == "C07_SAL_2018")

p_values <- subset_snps$P.value
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C07_SAL_2018_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY$FHBDON, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - p_values) * p_values) * effect_sizes^2

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype

# Print PVE for each SNP
print(pve_snp)

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBDON18 <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 3) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBDON 2018 SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBDON18)

#####
### FHBDON 2019
# p-values and effect sizes
subset_snps <- subset(data, Trait == "FarmCPU.FHBDON" & Plot == "C07_SAL_2019")
subset_snps

p_values <- subset_snps$P.value
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C07_SAL_2019_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY$FHBDON, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - p_values) * p_values) * effect_sizes^2

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype

# Print PVE for each SNP
print(pve_snp)

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBDON19 <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 3) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBDON 2019 SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBDON19)
#####
### FHBDON 2020
# p-values and effect sizes
subset_snps <- subset(data, Trait == "FarmCPU.FHBDON" & Plot == "C07_SAL_2020")
subset_snps

p_values <- subset_snps$P.value
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C07_SAL_2020_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY$FHBDON, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - p_values) * p_values) * effect_sizes^2

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype

# Print PVE for each SNP
print(pve_snp)

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBDON20 <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 3) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBDON 2020 SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBDON20)





# Load the required libraries
library(ggplot2)
library(gridExtra)

# Combine the plots into a single page
combined_page <- grid.arrange(FHBDON18, FHBDON19, FHBDON20, ncol = 2)

# Print the combined page
print(combined_page)


##########
data_examp <- data[c(22:29),c(1,2,10,11,12,14)]
write.csv(data_examp, file="~/Desktop/table.csv")
##########







### need to use MAF not pval

### FHBDON 2018
# maf and effect sizes

subset_snps <- subset(data, Trait == "FarmCPU.FHBDON" & Plot == "C07_SAL_2018")

maf <- subset_snps$MAF
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY18 <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C07_SAL_2018_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY18$FHBDON, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - maf) * maf) * effect_sizes^2

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype

# Print PVE for each SNP
print(pve_snp)

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBDON18 <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 2.5) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBDON 2018 SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBDON18)

rm(subset_snps, maf, effect_sizes, total_variance_phenotype, contribution_snp, pve_snp, original_pve, manual_pve, snps, pve_data)
#####

### FHBDON 2019
# p-values and effect sizes
subset_snps <- subset(data, Trait == "FarmCPU.FHBDON" & Plot == "C07_SAL_2019")

maf <- subset_snps$MAF
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY19 <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C07_SAL_2019_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY19$FHBDON, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - maf) * maf) * effect_sizes^2

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype

# Print PVE for each SNP
print(pve_snp)

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBDON19 <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 3) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBDON 2019 SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBDON19)

rm(subset_snps, maf, effect_sizes, total_variance_phenotype, contribution_snp, pve_snp, original_pve, manual_pve, snps, pve_data)
#####
### FHBDON 2020
# p-values and effect sizes
subset_snps <- subset(data, Trait == "FarmCPU.FHBDON" & Plot == "C07_SAL_2020")

maf <- subset_snps$MAF
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY20 <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C07_SAL_2020_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY20$FHBDON, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - maf) * maf) * effect_sizes^2

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype

# Print PVE for each SNP
print(pve_snp)

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBDON20 <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 2.2) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBDON 2020 SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBDON20)
rm(subset_snps, maf, effect_sizes, total_variance_phenotype, contribution_snp, pve_snp, original_pve, manual_pve, snps, pve_data)
#####
### FHBDON 2021 SAL
# p-values and effect sizes
subset_snps <- subset(data, Trait == "FarmCPU.FHBDON" & Plot == "C10_SAL_2021")

maf <- subset_snps$MAF
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY21S <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C10_SAL_2021_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY21S$FHBDON, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - maf) * maf) * effect_sizes^2

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype

# Print PVE for each SNP
print(pve_snp)

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBDON21S <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 2.2) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBDON 2021 Salina SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBDON21S)
rm(subset_snps, maf, effect_sizes, total_variance_phenotype, contribution_snp, pve_snp, original_pve, manual_pve, snps, pve_data)
#####
### FHBDON 2021 OLA
# p-values and effect sizes
subset_snps <- subset(data, Trait == "FarmCPU.FHBDON" & Plot == "C10_OLA_2021")

maf <- subset_snps$MAF
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY21O <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C10_OLA_2021_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY21O$FHBDON, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - maf) * maf) * effect_sizes^2

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype

# Print PVE for each SNP
print(pve_snp)

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBDON21O <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 2.2) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBDON 2021 Olathe SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBDON21O)
rm(subset_snps, maf, effect_sizes, total_variance_phenotype, contribution_snp, pve_snp, original_pve, manual_pve, snps, pve_data)
#####
### FHBDON 2022 SAL
# p-values and effect sizes
subset_snps <- subset(data, Trait == "FarmCPU.FHBDON" & Plot == "C10_SAL_2022")

maf <- subset_snps$MAF
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY22S <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C10_SAL_2022_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY21S$FHBDON, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - maf) * maf) * effect_sizes^2

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype

# Print PVE for each SNP
print(pve_snp)

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBDON22S <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 2.2) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBDON 2022 Salina SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBDON22S)
rm(subset_snps, maf, effect_sizes, total_variance_phenotype, contribution_snp, pve_snp, original_pve, manual_pve, snps, pve_data)
#####
### FHBDON 2022 OLA
# p-values and effect sizes
subset_snps <- subset(data, Trait == "FarmCPU.FHBDON" & Plot == "C10_OLA_2022")

maf <- subset_snps$MAF
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY22O <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C10_OLA_2022_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY22O$FHBDON, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - maf) * maf) * effect_sizes^2 #proportion of variance explained

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype #percentage of the total phenotypic variance

# Print PVE for each SNP
print(pve_snp)

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBDON22O <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 2.2) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBDON 2022 Olathe SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBDON22O)
rm(subset_snps, maf, effect_sizes, total_variance_phenotype, contribution_snp, pve_snp, original_pve, manual_pve, snps, pve_data)




#####
# Combine the plots into a single page
combined_page <- grid.arrange(FHBDON18, FHBDON19, FHBDON20, FHBDON21S, FHBDON21O, FHBDON22S, FHBDON22O, ncol = 2)

# Print the combined page
print(combined_page)




### FHBSPKINC 2018 : no sig snps
### FHBSPKINC 2019 : no sig snps
### FHBSPKINC 2020 : no sig snps
### FHBSPKINC 2021 SAL : no sig snps

### FHBSPKINC 2021 OLA
# p-values and effect sizes
subset_snps <- subset(data, Trait == "FarmCPU.FHBSPKINC" & Plot == "C10_OLA_2021")

maf <- subset_snps$MAF
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY21O <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C10_OLA_2021_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY21O$FHBSPKINC, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - maf) * maf) * effect_sizes^2

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBSPKINC21O <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 2.2) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBSPKINC 2021 Olathe SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBSPKINC21O)
rm(subset_snps, maf, effect_sizes, total_variance_phenotype, contribution_snp, pve_snp, original_pve, manual_pve, snps, pve_data)
#####
### FHBSPKINC 2022 SAL : no sig snps

### FHBSPKINC 2022 OLA
# p-values and effect sizes
subset_snps <- subset(data, Trait == "FarmCPU.FHBSPKINC" & Plot == "C10_OLA_2022")

maf <- subset_snps$MAF
effect_sizes <- subset_snps$effect

# Example total variance of the phenotype
myY22O <- read.table("/Volumes/Backup_Plus/AM/New_GWAS/2023.01.27/IWG_GWAS/data/Intermediate_Files/C10_OLA_2022_pheno.txt", head=TRUE)
total_variance_phenotype <- var(myY22O$FHBSPKINC, na.rm = T)

# Calculate the contribution of each SNP to variance explained
contribution_snp <- (2 * (1 - maf) * maf) * effect_sizes^2 #proportion of variance explained

# Calculate the PVE for each SNP
pve_snp <- contribution_snp / total_variance_phenotype #percentage of the total phenotypic variance

# Assuming you have vectors of originally calculated PVE and manually calculated PVE
original_pve <- subset_snps$Phenotype_Variance_Explained...  #  original PVE values
manual_pve <- pve_snp    #  manually calculated PVE values

# Assuming you have a vector of SNP names
snps <- paste0("SNP", seq_along(original_pve))

# Create a data frame with SNP names, original PVE, and transformed PVE
pve_data <- data.frame(SNP = snps, GAPIT_PVE = original_pve, Manual_PVE = manual_pve)

# Create a line plot with points
FHBSPKINC22O <- ggplot(pve_data, aes(x = SNP)) +
  # Original PVE points with labels
  geom_point(aes(y = GAPIT_PVE, color = "GAPIT PVE"), size = 3) +
  geom_text(aes(y = GAPIT_PVE, label = round(GAPIT_PVE, 2)), 
            hjust = -0.2, vjust = 0.5, color = "blue", size = 3) +
  # Transformed PVE points with labels
  geom_point(aes(y = Manual_PVE, color = "Manual PVE"), size = 3, shape = 2) +
  geom_text(aes(y = Manual_PVE, label = round(Manual_PVE, 12)), 
            hjust = -0.2, vjust = 1, color = "red", size = 2.2) +
  # Difference segments
  geom_segment(aes(y = GAPIT_PVE, yend = Manual_PVE), 
               color = "grey", size = 0.8, 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Axes and labels
  labs(x = "SNP", y = "PVE", title = "Original vs Manual PVE for Each FHBSPKINC 2022 Olathe SNP") +
  # Color scale
  scale_color_manual(values = c("GAPIT PVE" = "blue", "Manual PVE" = "red")) +
  # Rotate x-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank())

# Print the plot
print(FHBSPKINC22O)
rm(subset_snps, maf, effect_sizes, total_variance_phenotype, contribution_snp, pve_snp, original_pve, manual_pve, snps, pve_data)

#####
# Combine the plots into a single page
combined_page <- grid.arrange(FHBSPKINC21O, FHBSPKINC22O, ncol = 2)






###########################################
pve <- 2 * (1 - maf) * maf * effect_sizes^2
###########################################

Var(X) = E(x^2) - E(x)^2

average phenotype value for individuals with a certain genotype (aa, Aa, AA)

E(X^2) = frequency(aa) * average phenotype value for indivials with aa + frequency(Aa) * average phenotype value for individuals with Aa + f(AA) * average pheno value for indivuduals with AA

frequencies could be taken from hardy weinburg or actual frequncy of the genotype at the alleles


Var(X) = 2p(1-p)a^2


E(X) = f(aa) * -a + f(Aa) * 0 + f(AA) * a

a(f(A) - f(a))

a^2(F(A)+f(a))

[a(f(A) - f(a))]^2

#this is for additive effect, if GAPIT can give dominent effect (d) that can be used as the middle effect instead of 0
  #look at allele box plots for this 

E(X) = f(aa) * -a + f(AA) * a
E(X) = f(AA) * a - f(aa) * a  
E(X) = a(f(AA) - f(aa))
E(X)^2 = [a(f(AA) - f(aa))]^2 = a^2 * (f(AA) - f(aa))^2
E(X^2) = a^2(f(AA) - f(aa)) 

Var(X) = E(x^2) - E(x)^2
Var(X) = a^2(f(AA) - f(aa)) - [a(f(AA) - f(aa))]^2 



cols <- c("genotypic_frequency.AA", "genotypic_frequency.CC", "genotypic_frequency.TT", "genotypic_frequency.GG")
  
library(dplyr)

df <- data %>%
  mutate(
    max_column = apply(select(., all_of(cols)), 1, function(x) {
      if (all(x == 0)) NA  # handle all-zero case
      else {
        names(x)[which.max(x * (x > 0))]  # find the name of the max non-zero value
      }
    }),
    max_value = apply(select(., all_of(cols)), 1, function(x) {
      if (all(x == 0)) NA  # handle all-zero case
      else {
        max(x[x > 0])  # find the max non-zero value
      }
    }),
    min_column = apply(select(., all_of(cols)), 1, function(x) {
      if (all(x == 0)) NA  # handle all-zero case
      else {
        names(x)[which.min(x * (x > 0) + (x == 0) * max(x) * 100)]  # find the name of the min non-zero value
      }
    }),
    min_value = apply(select(., all_of(cols)), 1, function(x) {
      if (all(x == 0)) NA  # handle all-zero case
      else {
        min(x[x > 0])  # find the min non-zero value
      }
    })
  )

# now i have a max value (AA), min value (aa), and effect (a)
# Var(X) = a^2(f(AA) - f(aa)) - [a(f(AA) - f(aa))]^2 

df$var <- ((df$effect)^2 * (df$max_value - df$min_value)) - (df$effect * (df$max_value - df$min_value))^2

write.csv(df,'~/Desktop/output2.csv')

###
 
# Var(SNP) = beta^2_i * 2 * p_i * (1-p_i)
  # where: beta_i is the effect size of the ith SNP, and p_i is the allele frequency of the ith SNP

# Total Genetic Variance = SUM{Var(SNP_i)}
  # sum the variances explained by all significant SNPs to get total genetic variance explained by these SNPs

# Phenotypic Variance = Var(Phenotype)
  # sample variance of the phenotypic values

# PVE = Total Genetic Variance /  Phenotypic Variance
  # proportion of variance explained (PVE) by the significant SNPs

# 1. Calculate variance explained by each SNP
test <- data #copy of original data to play with
test$var <- (test$effect)^2 * 2 * test$MAF * (1-test$MAF)

# 2. Sum genetic variances (Total Genetic Variance)
TGV <- sum(test$var)

# 3. Calculate Phenotypic Variance












###
# Load necessary libraries
library(lme4)   # For mixed models
library(Matrix) # For matrix operations

# Assuming `y` is the phenotypic vector, `X` is the fixed effects matrix,
# `Z1` is the random effects design matrix for candidate regions,
# and `Z2` is the random effects design matrix for the rest of the genome

# Example data (replace with your actual data)
# y <- ...
y <- phenotype_data$y #phenotype values
# X <- ...
X <- model.matrix(~ PC1 + PC2 + PC3 + sex + farm, data = phenotype_data) 
# Z1 <- ...
# Z2 <- ...
candidate_indices <- ... # This should be a vector of column indices or logical vector indicating SNPs in candidate regions
Z1 <- as.matrix(genotype_data[, candidate_indices, with = FALSE])
Z2 <- as.matrix(genotype_data[, !candidate_indices, with = FALSE])

# Fit the mixed linear model
model <- lmer(y ~ X + (1|Z1) + (1|Z2), REML=TRUE)

# Extract variance components
var_components <- VarCorr(model)
sigma_u1 <- attr(var_components$Z1, "stddev")^2
sigma_u2 <- attr(var_components$Z2, "stddev")^2
sigma_e <- attr(var_components$Residual, "stddev")^2

# Calculate total genetic variance
sigma_u <- sigma_u1 + sigma_u2

# Calculate the proportion of genetic variance explained by candidate regions
PVE_genetic <- sigma_u1 / sigma_u

# Calculate the proportion of phenotypic variance explained by candidate regions
phenotypic_variance <- sigma_u + sigma_e
PVE_phenotypic <- sigma_u1 / phenotypic_variance

# Print results
cat("Proportion of Genetic Variance Explained by Candidate Regions:", PVE_genetic, "\n")
cat("Proportion of Phenotypic Variance Explained by Candidate Regions:", PVE_phenotypic, "\n")



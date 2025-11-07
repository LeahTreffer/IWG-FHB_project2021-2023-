library(flextable)
library(dplyr)
library(readxl)
library(officer)

Figure_3 <- read_excel("data/Intermediate_Files/Figure 3_Sig Loci.xlsx", 
                       sheet = "Table 1")
colnames(Figure_3)[6] <- "PVAL"
Figure_3$p.value <- formatC(Figure_3$PVAL, format = "e", digits = 3)
Figure_3$MAF <- formatC(Figure_3$MAF, digits = 4)
Figure_3$PVE <- formatC(Figure_3$PVE, digits = 6)
Figure_3$Effect <- formatC(Figure_3$Effect, digits = 6)
Figure_3$`NCBI Per.Ident` <- as.numeric(Figure_3$`NCBI Per.Ident`)
Figure_3$`NCBI Per.Ident` <- as.character(Figure_3$`NCBI Per.Ident`)
Figure_3[c(1:2,9:10), 12:14] <- "No Significant Hits"
Figure_3[1:2, 10:11] <- "No Proximal IWG Genes"
Figure_3 <- Figure_3[,c(1:5,7,15,8:14)]

myft <- flextable(Figure_3) 
myft <- theme_box(myft)  

std_border <- fp_border(color = "black", width = 1)

myft <- flextable(Figure_3) %>% 
  theme_box()%>%
  merge_at(
    i = 1:2, j = 1) %>%
  merge_at(
    i = 3:4, j = 1)%>%
  merge_at(
    i = 5:6, j = 1)%>%
  merge_at(
    i = 3:4, j = 2)%>%
  merge_at(
    i = 5:6, j = 3)%>%
  merge_at(
    i = 1:2, j = 4)%>%
  merge_at(
    i = 1:2, j = 5)%>%
  merge_at(
    i = 1:2, j = 6)%>%
  merge_at(
    i = 3:4, j = 4)%>%
  merge_at(
    i = 3:4, j = 5)%>%
  merge_at(
    i = 3:4, j = 6)%>%
  merge_at(
    i = 5:6, j = 4)%>%
  merge_at(
    i = 5:6, j = 5)%>%
  merge_at(
    i = 5:6, j = 6)%>%
  merge_at(
    i = 1:8, j = 10:11)%>%
  #merge_at(
    #i = 1:8, j = 11)%>%
  merge_at(
    i = 1:8, j = 12:14)%>%
  #merge_at(
    #i = 1:8, j = 13)%>%
  #merge_at(
    #i = 1:8, j = 14)%>%
  merge_at(
    i = 9, j=12:14)%>%
  merge_at(
    i = 10, j=12:14)%>%
  merge_at(
    i = 12:13, j = 1) %>%
  merge_at(
    i = 12:13, j = 3) %>%
  merge_at(
    i = 12:13, j = 4) %>%
  merge_at(
    i = 12:13, j = 5) %>%
  merge_at(
    i = 12:13, j = 6) %>%
  merge_at(
    i = 12:13, j = 10) %>%
  merge_at(
    i = 12:13, j = 11) %>%
  merge_at(
    i = 12:13, j = 12) %>%
  merge_at(
    i = 12:13, j = 13) %>%
  merge_at(
    i = 12:13, j = 14) %>%
  align(j = 1:14, align = "left", part = "all") %>%  # Align headers to the left
  padding(padding.top = 0, padding.bottom = 0, padding.left = 2.5, padding.right=0,part = "all")%>%
  plot
#export as image 1450x550



Figure_5 <- read_excel("data/Intermediate_Files/Figure 5_ Adjusted P.value.xlsx", 
                       sheet = "significant loci")

Figure_5$Effect <- formatC(Figure_5$Effect, digits = 5)
Figure_5$p.value <- formatC(Figure_5$P.value, format = "e", digits = 3)

Figure_5 <- Figure_5[,c(1,2,10,4:9)]

ft <- flextable(Figure_5)%>%
  theme_box()%>%
  merge_at(
    i = 2, j = 5:6) %>%
  merge_at(
    i = 5, j = 5:6) %>%
  merge_at(
    i = 1, j = 7:9) %>%
  merge_at(
    i = 5, j = 7:9) %>%
  align(j = 1:9, align = "left", part = "all") %>%  # Align headers to the left
  padding(padding.top = 0, padding.bottom = 0, padding.left=2.5, padding.right=0, part = "header")
  
set_table_properties(ft, width = 1, layout = "autofit")

ft_2 <- autofit(ft, add_w = 1, add_h = 0)
ft_2
# the dimensions are 'optimized here, so use the as guidlines to adjust

ft_3 <- ft_2 
ft_3[["body"]][["colwidths"]][["T. intermedium Description"]] <- 11.25
ft_3[["header"]][["colwidths"]][["T. intermedium Description"]] <- 11.25 
ft_3[["body"]][["colwidths"]][["NCBI Query Hit Description"]] <- 6.75 
ft_3[["header"]][["colwidths"]][["NCBI Query Hit Description"]] <- 6.75 
ft_3

#export as image 1400x351

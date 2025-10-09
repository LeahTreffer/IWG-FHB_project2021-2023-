library(ggplot2)

locs <- c("loc_a", "loc_b")
years <- 2017:2021
my_data <- expand.grid(loc  = locs, year = years) # x = site, y = year
my_data <- my_data[my_data$loc=="loc_a"|my_data$loc=="loc_b"&my_data$year>2019,]
measured_variable <- rnorm(28)
my_data2 <- my_data[rep(seq_len(nrow(my_data)), each = 4), ] # site, year
my_data3 <- cbind(my_data2, measured_variable) # site, year, value

ggplot(my_data3, aes(x=as.character(year), y=measured_variable, fill=loc)) + 
  geom_boxplot() + 
  #subset of data so we're just getting a couple of boxplots, no overplotting previous plots
  geom_boxplot(data = my_data3[my_data3$year>2019,],aes(x=as.character(year), 
                                                        y=measured_variable, fill=as.character(year)), lwd = 0.2, alpha = 0.4)

#line width adjusted with lwd argument, transparency adjuted with alpha



Leah_data <- data.frame (Site  = c("Salina", "Salina", "Salina", "Salina", "Olathe", "Salina", "Olathe"),
                  Year = c("2018", "2019", "2020", "2021", "2021", "2022", "2022")
)
#OR
Site <- c("Salina", "Olathe")
Year <- 2018:2022
Leah_data <- expand.grid(Site  = Site, Year = Year) # x = site, y = year
Leah_data <- Leah_data[Leah_data$Site=="Salina"|Leah_data$Site=="Olathe"&Leah_data$Year>2020,]

Leah_data2 <- pheno_table[,c("Site","Year")]

Leah_data3 <- pheno_table[,c("Site","Year","Cycle","FHBDON")]

 ggplot(Leah_data3, aes(x=as.character(Year), y=FHBDON, color=Site, fill= Cycle)) + 
  geom_boxplot() +
  scale_color_manual(name = "Site", values = c("Salina" = "#882255", "Olathe" = "#DDCC77")) +
  scale_fill_manual(name = "Cycle", values = c("C07" = "#44AA99", "C10" = "#CC6677"))+
  geom_segment(aes(x = as.numeric(as.factor(Year)), xend = as.numeric(as.factor(Year)),
        y = red_line_values, yend = red_line_values),
    color = "red", size = 1
  )

 
d=data.frame(x=c(0.5,1,1.5,2,2,2.5,2.5), y=c(20,50,45,100,100,75,75), vx=c(1,1.5,2,2.5,2.5,3,3), vy=c(20,50,45,100,100,75,75)) # these would be mean annual wheat scores

red_line_values <- c(20, 50, 45, 100, 100, 75, 75)

g + geom_segment(aes(x = c(0,0.2,0.4,0.6,0.6,0.8,0.8), y = c(20,50,45,100,100,75,75), xend = c(0.2,0.4,0.6,0.8,0.8,1,1), yend = c(20,50,45,100,100,75,75), colour = "red"))



ggplot(Leah_data3, aes(x = as.character(Year), y = FHBDON, color = Site, fill = Cycle)) + 
  geom_boxplot() +
  scale_color_manual(name = "Site", values = c("Salina" = "#882255", "Olathe" = "#DDCC77")) +
  scale_fill_manual(name = "Cycle", values = c("C07" = "#44AA99", "C10" = "#CC6677")) +
  geom_segment(
    aes(x = as.numeric(as.factor(Year)), xend = as.numeric(as.factor(Year)),
        y = red_line_values, yend = red_line_values),
    color = "red", size = 1
  )





#+ 
  #subset of data so we're just getting a couple of boxplots, no overplotting previous plots
  geom_boxplot(data = Leah_data3[Leah_data3$Year>2020,],aes(x=as.character(Year), 
                                                        y=FHBDON, fill=as.character(Year)), lwd = 0.2, alpha = 0.4)+
  scale_color_manual(name = "Site", values = c("Salina" = "#882255", "Olathe" = "#DDCC77"))
  

Leah_data3


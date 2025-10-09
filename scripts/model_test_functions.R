#Generalized linear mixed model (GLMM)
# Treat Time as random: 
#in R lmer and glmer (1|'factor') is a limiting factor
# poisson family = log, identity, or sqrt
# use Values, not Value (Values is an integer)

AIC_model_selection_glmer<- function(table) {
# Site and Cycle as fixed, time and site:time interaction as random effect
glmm_poisson1=glmer(Values ~Site + Cycle + (1|Time) + (1|Site:Time),
                    data = table, family = poisson(link = "sqrt"))  # Poisson 
# Cycle and Site as fixed, count and Time as random effects
glmm_poisson2=glmer(Values ~ Cycle + Site + (1|count)+ (1|Time), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle  as fixed, Site, count, and Time as random effects
glmm_poisson3=glmer(Values ~ Cycle + (1|Site) + (1|count) + (1|Time), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson4=glmer(Values ~ Cycle + Time + Site + (1|count), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson5=glmer(Values ~ Cycle + Site + Time + (1|count), 
                    data = table, family = poisson(link = "sqrt"))  
# Site and Cycle and Time as fixed; site:time interaction as random effect
glmm_poisson6=glmer(Values ~Site + Cycle + Time + (1|Site:Time),
                    data = table, family = poisson(link = "sqrt"))  # Poisson 
# Cycle, Site, and Time as fixed, count and site:time as random
glmm_poisson7=glmer(Values ~ Cycle + Time + Site + (1|count) + (1|Site:Time), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson8=glmer(Values ~ Cycle + Site + Time + (1|count) + (1|Site:Time), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson9=glmer(Values ~Site + Cycle + Time + (1|count), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson10=glmer(Values ~ Site + Time + Cycle + (1|count), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson11=glmer(Values ~ Time + Cycle + Site + (1|count), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson12=glmer(Values ~ Time + Site + Cycle + (1|count), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson13=glmer(Values ~  Site + Cycle + Time + (1|count) + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson14=glmer(Values ~ Site + Time + Cycle + (1|count) + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson15=glmer(Values ~ Time + Cycle + Site + (1|count) + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson16=glmer(Values ~ Time + Site + Cycle + (1|count) + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count and site:time as random
glmm_poisson17=glmer(Values ~ Cycle + Time + Site + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson18=glmer(Values ~ Cycle + Site + Time + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson19=glmer(Values ~  Site + Cycle + Time + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson20=glmer(Values ~ Site + Time + Cycle + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson21=glmer(Values ~ Time + Cycle + Site + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson22=glmer(Values ~ Time + Site + Cycle + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Models Added 1/4/2022
glmm_poisson23=glmer(Values ~Site*Cycle*(1|Time)*(1|Site:Time),
                     data = table, family = poisson(link = "sqrt"))  # Poisson 
# Cycle and Site as fixed, count and Time as random effects
glmm_poisson24=glmer(Values ~ Cycle*Site*(1|count)*(1|Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle  as fixed, Site, count, and Time as random effects
glmm_poisson25=glmer(Values ~ Cycle*(1|Site)*(1|count)*(1|Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson26=glmer(Values ~ Cycle*Time*Site*(1|count), 
                     data = table, family = poisson(link = "sqrt"))  
# Site and Cycle as fixed, time and site:time interaction as random effect
glmm_poisson27=glmer(Values ~Site*Cycle + (1|Time) + (1|Site:Time),
                     data = table, family = poisson(link = "sqrt"))  # Poisson 
# Cycle and Site as fixed, count and Time as random effects
glmm_poisson28=glmer(Values ~ Cycle*Site + (1|count) + (1|Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Site and Cycle and Time as fixed; site:time interaction as random effect
glmm_poisson29=glmer(Values ~Site*Cycle*Time + (1|Site:Time),
                     data = table, family = poisson(link = "sqrt"))  # Poisson 
# Site and Cycle and Time as fixed; site:time interaction as random effect
glmm_poisson30=glmer(Values ~Site*Cycle*Time*(1|Site:Time),
                     data = table, family = poisson(link = "sqrt"))  # Poisson 
# Cycle, Site, and Time as fixed, count and site:time as random
glmm_poisson31=glmer(Values ~ Cycle*Time*Site*(1|count)*(1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count and site:time as random
glmm_poisson32=glmer(Values ~ Cycle*Time*Site + (1|count) + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))
glmm_poisson33=glmer(Values ~ Cycle*Site*(1|Time), 
                     data = table, family = poisson(link = "sqrt"))

# Find the best-fit model # https://www.scribbr.com/statistics/anova-in-r/
model.set <- list(glmm_poisson1, glmm_poisson2, glmm_poisson3, glmm_poisson4, glmm_poisson5, glmm_poisson6, glmm_poisson7, glmm_poisson8, glmm_poisson9, glmm_poisson10, glmm_poisson11, glmm_poisson12, glmm_poisson13, glmm_poisson14, glmm_poisson15, glmm_poisson16, glmm_poisson17, glmm_poisson18, glmm_poisson19, glmm_poisson20, glmm_poisson21, glmm_poisson22, glmm_poisson23, glmm_poisson24, glmm_poisson25, glmm_poisson26, glmm_poisson27, glmm_poisson28, glmm_poisson29, glmm_poisson30, glmm_poisson31, glmm_poisson32, glmm_poisson33)
model.names <- c("glmm_poisson1", "glmm_poisson2", "glmm_poisson3", "glmm_poisson4", "glmm_poisson5", "glmm_poisson6", "glmm_poisson7","glmm_poisson8", "glmm_poisson9", "glmm_poisson10", "glmm_poisson11", "glmm_poisson12", "glmm_poisson13", "glmm_poisson14", "glmm_poisson15", "glmm_poisson16", "glmm_poisson17", "glmm_poisson18", "glmm_poisson19", "glmm_poisson20", "glmm_poisson21", "glmm_poisson22", "glmm_poisson23", "glmm_poisson24", "glmm_poisson25", "glmm_poisson26", "glmm_poisson27", "glmm_poisson28", "glmm_poisson29", "glmm_poisson30", "glmm_poisson31", "glmm_poisson32", "glmm_poisson33")
aic_table <- aictab(model.set, modnames = model.names)

ind = which(aic_table[,'AICc'] == min(aic_table[,'AICc'])) # which in AIC table are the lowest scores
try <- aic_table[rownames(aic_table)[ind],] # subset AIC table to just have the rows of models with top AICc scores
day = rownames(aic_table)[ind] 
names = try$Modnames # list of model names 
return(day)

}

# can be run separately using aic_table<-AIC_model_selection_glmer(table) then best_model(aic_table)
best_model <- function(aic_table){

ind = which(aic_table[,'AICc'] == min(aic_table[,'AICc'])) # which in AIC table are the lowest scores
try <- aic_table[rownames(aic_table)[ind],] # subset AIC table to just have the rows of models with top AICc scores
day = rownames(aic_table)[ind] 
names = try$Modnames # list of model names 
return(day)
}

best_models_glmer <- function(a){
# Site and Cycle as fixed, time and site:time interaction as random effect
glmm_poisson1=glmer(Values ~Site + Cycle + (1|Time) + (1|Site:Time),
                    data = table, family = poisson(link = "sqrt"))  # Poisson 
# Cycle and Site as fixed, count and Time as random effects
glmm_poisson2=glmer(Values ~ Cycle + Site + (1|count)+ (1|Time), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle  as fixed, Site, count, and Time as random effects
glmm_poisson3=glmer(Values ~ Cycle + (1|Site) + (1|count) + (1|Time), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson4=glmer(Values ~ Cycle + Time + Site + (1|count), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson5=glmer(Values ~ Cycle + Site + Time + (1|count), 
                    data = table, family = poisson(link = "sqrt"))  
# Site and Cycle and Time as fixed; site:time interaction as random effect
glmm_poisson6=glmer(Values ~Site + Cycle + Time + (1|Site:Time),
                    data = table, family = poisson(link = "sqrt"))  # Poisson 
# Cycle, Site, and Time as fixed, count and site:time as random
glmm_poisson7=glmer(Values ~ Cycle + Time + Site + (1|count) + (1|Site:Time), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson8=glmer(Values ~ Cycle + Site + Time + (1|count) + (1|Site:Time), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson9=glmer(Values ~Site + Cycle + Time + (1|count), 
                    data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson10=glmer(Values ~ Site + Time + Cycle + (1|count), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson11=glmer(Values ~ Time + Cycle + Site + (1|count), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson12=glmer(Values ~ Time + Site + Cycle + (1|count), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson13=glmer(Values ~  Site + Cycle + Time + (1|count) + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson14=glmer(Values ~ Site + Time + Cycle + (1|count) + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson15=glmer(Values ~ Time + Cycle + Site + (1|count) + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson16=glmer(Values ~ Time + Site + Cycle + (1|count) + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count and site:time as random
glmm_poisson17=glmer(Values ~ Cycle + Time + Site + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson18=glmer(Values ~ Cycle + Site + Time + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson19=glmer(Values ~  Site + Cycle + Time + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson20=glmer(Values ~ Site + Time + Cycle + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson21=glmer(Values ~ Time + Cycle + Site + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count  and site:time as random
glmm_poisson22=glmer(Values ~ Time + Site + Cycle + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Models Added 1/4/2022
glmm_poisson23=glmer(Values ~Site*Cycle*(1|Time)*(1|Site:Time),
                     data = table, family = poisson(link = "sqrt"))  # Poisson 
# Cycle and Site as fixed, count and Time as random effects
glmm_poisson24=glmer(Values ~ Cycle*Site*(1|count)*(1|Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle  as fixed, Site, count, and Time as random effects
glmm_poisson25=glmer(Values ~ Cycle*(1|Site)*(1|count)*(1|Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count as random
glmm_poisson26=glmer(Values ~ Cycle*Time*Site*(1|count), 
                     data = table, family = poisson(link = "sqrt"))  
# Site and Cycle as fixed, time and site:time interaction as random effect
glmm_poisson27=glmer(Values ~Site*Cycle + (1|Time) + (1|Site:Time),
                     data = table, family = poisson(link = "sqrt"))  # Poisson 
# Cycle and Site as fixed, count and Time as random effects
glmm_poisson28=glmer(Values ~ Cycle*Site + (1|count) + (1|Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Site and Cycle and Time as fixed; site:time interaction as random effect
glmm_poisson29=glmer(Values ~Site*Cycle*Time + (1|Site:Time),
                     data = table, family = poisson(link = "sqrt"))  # Poisson 
# Site and Cycle and Time as fixed; site:time interaction as random effect
glmm_poisson30=glmer(Values ~Site*Cycle*Time*(1|Site:Time),
                     data = table, family = poisson(link = "sqrt"))  # Poisson 
# Cycle, Site, and Time as fixed, count and site:time as random
glmm_poisson31=glmer(Values ~ Cycle*Time*Site*(1|count)*(1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))  
# Cycle, Site, and Time as fixed, count and site:time as random
glmm_poisson32=glmer(Values ~ Cycle*Time*Site + (1|count) + (1|Site:Time), 
                     data = table, family = poisson(link = "sqrt"))
glmm_poisson33=glmer(Values ~ Cycle*Site*(1|Time), 
                     data = table, family = poisson(link = "sqrt"))


model.set <- list(glmm_poisson1, glmm_poisson2, glmm_poisson3, glmm_poisson4, glmm_poisson5, glmm_poisson6, glmm_poisson7, glmm_poisson8, glmm_poisson9, glmm_poisson10, glmm_poisson11, glmm_poisson12, glmm_poisson13, glmm_poisson14, glmm_poisson15, glmm_poisson16, glmm_poisson17, glmm_poisson18, glmm_poisson19, glmm_poisson20, glmm_poisson21, glmm_poisson22, glmm_poisson23, glmm_poisson24, glmm_poisson25, glmm_poisson26, glmm_poisson27, glmm_poisson28, glmm_poisson29, glmm_poisson30, glmm_poisson31, glmm_poisson32, glmm_poisson33)

aa <- (model.set[a])
aaa <- (summary(model.set[[a]]))

my_list <- list(aa,aaa)
return(my_list)

}


#For Analysis of Cycle07
AIC_model_selection_C07<- function(table) {
  one.way <- aov(Value ~ Time, data = table)
  two.way <- aov(Value ~ Time + count, data = table)
  interaction <- aov(Value ~ Time*count, data = table)  
  
  model.set <- list(one.way, two.way, interaction)
  model.names <- c("one.way", "two.way", "interaction")
  aic_table <- aictab(model.set, modnames = model.names)
  
  ind = which(aic_table[,'AICc'] == min(aic_table[,'AICc'])) # which in AIC table are the lowest scores
  try <- aic_table[rownames(aic_table)[ind],] # subset AIC table to just have the rows of models with top AICc scores
  day = rownames(aic_table)[ind] 
  names = try$Modnames # list of model names 
  return(day)
  
}

best_models_C07 <- function(a){
  one.way <- aov(Value ~ Time, data = table)
  two.way <- aov(Value ~ Time + count, data = table)
  interaction <- aov(Value ~ Time*count, data = table) 
  
  model.set <- list(one.way, two.way, interaction)
  
  aa <- (model.set[a])
  aaa <- (summary(model.set[[a]]))
  
  my_list <- list(aa,aaa)
  return(my_list)
  
}


AIC_model_selection_C10 <- function(table){
  one.way <- aov(Value ~ Time, data = table)
  one.way2 <- aov(Value ~ Site, data = table)
  one.way3 <- aov(Value ~ count, data = table)
  two.way <- aov(Value ~ Time + Site, data = table)
  interaction <- aov(Value ~ Time*Site, data = table)
  blocking <- aov(Value ~ Time + Site + count, data = table)
  blocking2 <- aov(Value ~ Time + count + Site, data = table)
  blocking3 <- aov(Value ~ Site + Time + count, data = table)
  blocking4 <- aov(Value ~ Site + count + Time, data = table)
  blocking5 <- aov(Value ~ count + Time + Site, data = table)
  blocking6 <- aov(Value ~ count + Site + Time, data = table)
  unknown <- aov(Value ~ Time*Site + count, data = table)
  
  model.set <- list(one.way, one.way2, one.way3, two.way, interaction, blocking, blocking2, blocking3, blocking4, blocking5, blocking6, unknown)
  model.names <- c("one.way", "one.way2", "one.way3", "two.way", "interaction", "blocking", "blocking2", "blocking3", "blocking4", "blocking5", "blocking6", "unknown")
  aic_table <- aictab(model.set, modnames = model.names)
  
  ind = which(aic_table[,'AICc'] == min(aic_table[,'AICc'])) # which in AIC table are the lowest scores
  try <- aic_table[rownames(aic_table)[ind],] # subset AIC table to just have the rows of models with top AICc scores
  day = rownames(aic_table)[ind] 
  names = try$Modnames # list of model names 
  return(day)
}

best_models_C10 <- function(a){
  one.way <- aov(Value ~ Time, data = table)
  one.way2 <- aov(Value ~ Site, data = table)
  one.way3 <- aov(Value ~ count, data = table)
  two.way <- aov(Value ~ Time + Site, data = table)
  interaction <- aov(Value ~ Time*Site, data = table)
  blocking <- aov(Value ~ Time + Site + count, data = table)
  blocking2 <- aov(Value ~ Time + count + Site, data = table)
  blocking3 <- aov(Value ~ Site + Time + count, data = table)
  blocking4 <- aov(Value ~ Site + count + Time, data = table)
  blocking5 <- aov(Value ~ count + Time + Site, data = table)
  blocking6 <- aov(Value ~ count + Site + Time, data = table)
  unknown <- aov(Value ~ Time*Site + count, data = table)
  
  model.set <- list(one.way, one.way2, one.way3, two.way, interaction, blocking, blocking2, blocking3, blocking4, blocking5, blocking6, unknown)
  
  aa <- (model.set[a])
  aaa <- (summary(model.set[[a]]))
  
  my_list <- list(aa,aaa)
  return(my_list)

}

#library(sjPlot)
#library(sjmisc)
#library(sjlabelled)
#tab_model(model)

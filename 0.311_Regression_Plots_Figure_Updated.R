library(tidyverse)
library(ggpubr)


j <- ggscatter(
  C07, x = "FHBDISIND", y = "FHBDON", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 30) +
  stat_regline_equation(label.y = 27)+
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

k <- ggscatter(
  C07, x = "FHBSPKINC", y = "FHBDON", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 30) +
  stat_regline_equation(label.y = 28) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))  # This line removes the legend

l <- ggscatter(
  C07, x = "FHBDISIND", y = "FHBSPKINC", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 200) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
m <- ggscatter(
  C07, x = "FHBD3G", y = "FHBDON", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 30) +
  stat_regline_equation(label.y = 28)+
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
o <- ggscatter(
  C10, x = "FHBDISIND", y = "FHBDON", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 150) +
  stat_regline_equation(label.y = 139)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  labs(color = NULL) +  # Remove the legend title for color
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 20))
p <- ggscatter(
  C10, x = "FHBSPKINC", y = "FHBDON", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 150) +
  stat_regline_equation(label.y = 140)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
q <- ggscatter(
  C10, x = "FHBDISIND", y = "FHBSPKINC", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 150) +
  stat_regline_equation(label.y = 140)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
r <- ggscatter(
  C10, x = "FHBD3G", y = "FHBDON", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 175) +
  stat_regline_equation(label.y = 165)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))


figure <- ggarrange(j, o, k, p, l, q, m, r,
                    ncol = 2, nrow = 4, 
                    common.legend = FALSE)
figure
#export as image:
#2010x1395

####################################################################

j <- ggscatter(
  C07, x = "FHBDISIND", y = "FHBDON", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 30) +
  stat_regline_equation(label.y = 27)+
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

k <- ggscatter(
  C07, x = "FHBSPKINC", y = "FHBDON", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 30) +
  stat_regline_equation(label.y = 28) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))  # This line removes the legend

l <- ggscatter(
  C07, x = "FHBDISIND", y = "FHBSPKINC", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 200) +
  stat_regline_equation(label.y = 185)+
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
m <- ggscatter(
  C07, x = "FHBD3G", y = "FHBDON", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 30) +
  stat_regline_equation(label.y = 28)+
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
n <- ggscatter(
  C07, x = "FHBDISIND", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 165) +
  stat_regline_equation(label.y = 160)+
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
t <- ggscatter(
  C07, x = "FHBDISIND", y = "SPKYLD", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 5) +
  stat_regline_equation(label.y = 4)+
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
o <- ggscatter(
  C10, x = "FHBDISIND", y = "FHBDON", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 150) +
  stat_regline_equation(label.y = 139)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  labs(color = NULL) +  # Remove the legend title for color
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 20))
p <- ggscatter(
  C10, x = "FHBSPKINC", y = "FHBDON", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 150) +
  stat_regline_equation(label.y = 140)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
q <- ggscatter(
  C10, x = "FHBDISIND", y = "FHBSPKINC", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 150) +
  stat_regline_equation(label.y = 140)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
r <- ggscatter(
  C10, x = "FHBD3G", y = "FHBDON", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 175) +
  stat_regline_equation(label.y = 165)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
s <- ggscatter(
  C10, x = "FHBDISIND", y = "PTHT", 
  color = "phenotype_year_site", palette = c("#490092","lightgrey","#b66dff","#999999"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 150) +
  stat_regline_equation(label.y = 140)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))


figure <- ggarrange(j, o, k, p, l, q, m, r, n, s, t,
                    ncol = 2, nrow = 6, 
                    common.legend = FALSE)
figure
#export as image:
#2010x1395


########################################################################

j <- ggscatter(
  C07, x = "FHBDON", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 212) +
  stat_regline_equation(label.y = 195)+
  theme(axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22))

k <- ggscatter(
  C07, x = "FHBD3G", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 205) +
  stat_regline_equation(label.y = 193) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))  # This line removes the legend
m <- ggscatter(
  C07, x = "FHBDISIND", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 205) +
  stat_regline_equation(label.y = 193)+
  theme(legend.position = "none", 
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))  # This line removes the legend
n <- ggscatter(
  C07, x = "FHBSPKINC", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 205) +
  stat_regline_equation(label.y = 193)+
  theme(legend.position = "none", 
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))  # This line removes the legend
o <- ggscatter(
  C10, x = "FHBDON", y = "PTHT", 
  color = "phenotype_year_site", palette = c("#004949","#999999","#490092","#b66dff"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 198) +
  stat_regline_equation(label.y = 180)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  labs(color = NULL) +  # Remove the legend title for color
  theme(axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 22))
p <- ggscatter(
  C10, x = "FHBD3G", y = "PTHT", 
  color = "phenotype_year_site", palette = c("#004949","#999999","#490092","#b66dff"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 156) +
  stat_regline_equation(label.y = 144)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))
r <- ggscatter(
  C10, x = "FHBDISIND", y = "PTHT", 
  color = "phenotype_year_site", palette = c("#004949","#999999","#490092","#b66dff"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 172) +
  stat_regline_equation(label.y = 160)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))
s <- ggscatter(
  C10, x = "FHBSPKINC", y = "PTHT", 
  color = "phenotype_year_site", palette = c("#004949","#999999","#490092","#b66dff"),
  add = "reg.line"
  
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 160) +
  stat_regline_equation(label.y = 148)+
  labs(x = NULL, y = NULL)+  # This line removes the axis labels
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))
t <- ggscatter(
  C07, x = "SPKYLD", y = "PTHT", 
  color = "phenotype_year", palette = c("#920000", "#DDCC77", "#006ddb"),
  add = "reg.line"
) +
  facet_wrap(~phenotype_year) +
  stat_cor(label.y = 195) +
  stat_regline_equation(label.y = 182)+
  theme(legend.position = "none")+  # This line removes the legend
  theme(legend.position = "none",  # This line removes the legend
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))

figure <- ggarrange(j, o, k, p, m, r, n, s, t, NULL,
                    #labels = c("A", "B", "C"),
                    ncol = 2, nrow = 5, 
                    common.legend = FALSE)
figure
#2680x1856

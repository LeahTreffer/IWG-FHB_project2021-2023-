library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(purrr)
library(broom)

# made a NOAA request for weather data from 2018-2022 from Salina and Olathe KS
# https://www.noaa.gov/tools-and-resources/weather-and-climate-resources#historic
# https://www.ncdc.noaa.gov/cdo-web/datatools/lcd

# read in raw NOAA data
file <- read.csv('data/KS_weather.csv')
colnames(file)

# columns that sound promising
# "STATION", "DATE", "REPORT_TYPE", "SOURCE", "DailyAverageDewPointTemperature", 
# "DailyAverageDryBulbTemperature", "DailyAverageRelativeHumidity", "DailyAverageWetBulbTemperature",
# "DailyMaximumDryBulbTemperature", "DailyMinimumDryBulbTemperature", "DailyPrecipitation",
# "DailyWeather", "MonthlyAverageRH", "MonthlyDaysWithGT001Precip", "MonthlyDaysWithGT010Precip",
# "MonthlyDaysWithGT32Temp", "MonthlyDaysWithGT90Temp", "MonthlyDaysWithLT0Temp",
# "MonthlyDaysWithLT32Temp", "MonthlyMeanTemperature", "MonthlyMinimumTemperature", 
# "MonthlyTotalLiquidPrecipitation", "REPORT_TYPE.1", "SOURCE.1"

# rename station identifiers to something recognizable  
file$STATION[file$STATION == "72458603919"] <- "Salina"
file$STATION[file$STATION == "72447593909"] <- "Olathe"
file$STATION <- as.factor(file$STATION)

# Salina = WBAN:03919 = 72458603919
Salina <- file[file$STATION == 'Salina', c("STATION", "DATE", "REPORT_TYPE", "SOURCE", "REPORT_TYPE.1", "SOURCE.1", 
                                      "DailyAverageDryBulbTemperature",
                                      "DailyMaximumDryBulbTemperature", "DailyMinimumDryBulbTemperature", "DailyPrecipitation",
                                      "DailyWeather", "MonthlyDaysWithGT001Precip", "MonthlyDaysWithGT010Precip",
                                      "MonthlyDaysWithGT32Temp", "MonthlyDaysWithGT90Temp", "MonthlyDaysWithLT0Temp",
                                      "MonthlyDaysWithLT32Temp", "MonthlyMeanTemperature", "MonthlyMaximumTemperature",
                                      "MonthlyMinimumTemperature", "MonthlyTotalLiquidPrecipitation")]

# Olathe = WBAN:93909 = 72447593909
Olathe <- file[file$STATION == 'Olathe', c("STATION", "DATE", "REPORT_TYPE", "SOURCE", "REPORT_TYPE.1", "SOURCE.1", 
                                     "DailyAverageDryBulbTemperature",
                                     "DailyMaximumDryBulbTemperature", "DailyMinimumDryBulbTemperature", "DailyPrecipitation",
                                     "DailyWeather", "MonthlyDaysWithGT001Precip", "MonthlyDaysWithGT010Precip",
                                     "MonthlyDaysWithGT32Temp", "MonthlyDaysWithGT90Temp", "MonthlyDaysWithLT0Temp",
                                     "MonthlyDaysWithLT32Temp", "MonthlyMeanTemperature", "MonthlyMaximumTemperature",
                                     "MonthlyMinimumTemperature", "MonthlyTotalLiquidPrecipitation")]

### Monthly Temperature

MeanTemp <- rbind(Salina[,c("STATION", "DATE", "REPORT_TYPE", "SOURCE", "MonthlyMeanTemperature")], Olathe[,c("STATION", "DATE", "REPORT_TYPE", "SOURCE", "MonthlyMeanTemperature")]) %>% filter(!is.na(MonthlyMeanTemperature)) 

# Ensure DATE is a Date type (if it's not already)
MeanTemp$DATE <- as.Date(MeanTemp$DATE)

ggplot(MeanTemp, aes(x = DATE, y = MonthlyMeanTemperature, color = STATION)) +
  geom_line(size = 1) +
  labs(title = "Monthly Mean Temperature by Location",
       x = "Month",
       y = "Mean Temperature (°F)",
       color = "Station") +
  theme_minimal()


### Monthly Precipitation

MeanPrecip <- rbind(Salina[,c("STATION", "DATE", "REPORT_TYPE", "SOURCE", "MonthlyTotalLiquidPrecipitation")], Olathe[,c("STATION", "DATE", "REPORT_TYPE", "SOURCE", "MonthlyTotalLiquidPrecipitation")]) %>% filter(!is.na(MonthlyTotalLiquidPrecipitation)) 

# Ensure DATE is a Date type (if it's not already)
MeanPrecip$DATE <- as.Date(MeanPrecip$DATE)

MeanPrecip$MonthlyTotalLiquidPrecipitation <- gsub("s", "", MeanPrecip$MonthlyTotalLiquidPrecipitation)
MeanPrecip$MonthlyTotalLiquidPrecipitation <- gsub("T", "", MeanPrecip$MonthlyTotalLiquidPrecipitation)
MeanPrecip$MonthlyTotalLiquidPrecipitation <- as.numeric(MeanPrecip$MonthlyTotalLiquidPrecipitation)
MeanPrecip <- na.omit(MeanPrecip) 

ggplot(MeanPrecip, aes(x = DATE, y = MonthlyTotalLiquidPrecipitation, color = STATION)) +
  geom_line(size = 1) +
  labs(title = "Monthly Mean Precipitation by Location",
       x = "Month",
       y = "Mean Precipitation",
       color = "Station") +
  theme_minimal()

# Create a data frame of shaded date ranges (May 1 to Sept 30 of each year)
# Get min and max y-values for shaded area height
y_min <- min(MeanPrecip$MonthlyTotalLiquidPrecipitation, na.rm = TRUE)
y_max <- max(MeanPrecip$MonthlyTotalLiquidPrecipitation, na.rm = TRUE)

shaded_ranges <- MeanPrecip %>%
  mutate(Year = year(DATE)) %>%
  distinct(Year) %>%
  mutate(
    xmin = as.Date(paste0(Year, "-05-01")),
    xmax = as.Date(paste0(Year, "-09-30")),
    ymin = y_min,
    ymax = y_max
  )

ggplot() +
  geom_rect(
    data = shaded_ranges,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "gray80", alpha = 0.4
  ) +
  geom_line(
    data = MeanPrecip,
    aes(x = DATE, y = MonthlyTotalLiquidPrecipitation, color = STATION),
    size = 1
  ) +
  labs(
    title = "Monthly Mean Precipitation by Station (May–Sep Highlighted)",
    x = "Date",
    y = "Mean Precipitation",
    color = "Station"
  ) +
  theme_minimal()

### DailyPrecipitation

DailyPrecip <- rbind(Salina[,c("STATION", "DATE", "REPORT_TYPE", "SOURCE", "DailyPrecipitation")], Olathe[,c("STATION", "DATE", "REPORT_TYPE", "SOURCE", "DailyPrecipitation")]) %>% filter(!is.na(DailyPrecipitation)) 

# Ensure DATE is a Date type (if it's not already)
DailyPrecip$DATE <- as.Date(DailyPrecip$DATE)

unique(DailyPrecip$DailyPrecipitation)

DailyPrecip$DailyPrecipitation <- gsub("s", "", DailyPrecip$DailyPrecipitation)
DailyPrecip$DailyPrecipitation <- gsub("T", "", DailyPrecip$DailyPrecipitation)
DailyPrecip$DailyPrecipitation <- as.numeric(DailyPrecip$DailyPrecipitation)
DailyPrecip <- na.omit(DailyPrecip) 

ggplot(DailyPrecip, aes(x = DATE, y = DailyPrecipitation, color = STATION)) +
  geom_line(size = 1) +
  labs(title = "Daily Precipitation by Location",
       x = "Month",
       y = "Precipitation",
       color = "Station") +
  theme_minimal()

# Create a data frame of shaded date ranges (May 1 to Sept 30 of each year)
# Get min and max y-values for shaded area height
y_min <- min(DailyPrecip$DailyPrecipitation, na.rm = TRUE)
y_max <- max(DailyPrecip$DailyPrecipitation, na.rm = TRUE)

shaded_ranges <- DailyPrecip %>%
  mutate(Year = year(DATE)) %>%
  distinct(Year) %>%
  mutate(
    xmin = as.Date(paste0(Year, "-05-01")),
    xmax = as.Date(paste0(Year, "-09-30")),
    ymin = y_min,
    ymax = y_max
  )

ggplot() +
  geom_rect(
    data = shaded_ranges,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "gray80", alpha = 0.4
  ) +
  geom_line(
    data = DailyPrecip,
    aes(x = DATE, y = DailyPrecipitation, color = STATION),
    size = 1
  ) +
  labs(
    title = "Daily Precipitation by Location (May–Sep Highlighted)",
    x = "Date",
    y = "Precipitation",
    color = "Station"
  ) +
  theme_minimal()


### Summary Stats
## Function for value difference between sites in each month
difference_table_function <- function(data, Min_col = NULL, Max_col = NULL, Mean_col = NULL, table_name) {
  
  # Function to compute difference and identify "winning" station
  run_diff_tests <- function(data, value_col, comparison = c("higher", "lower")) {
    comparison <- match.arg(comparison)
    
    data %>%
      select(DATE, STATION, value = !!sym(value_col)) %>%
      filter(!is.na(value)) %>%
      group_by(DATE) %>%
      filter(n_distinct(STATION) == 2) %>%
      summarise(
        diff = abs(diff(value)),
        Station = STATION[which.max(value)],
        variable = value_col,
        .groups = "drop"
      )
  }
  
  # Initialize an empty list to collect results
  results <- list()
  
  # Run for each column if it's provided
  if (!is.null(Min_col) && Min_col %in% names(data)) {
    results[["min_diff"]] <- run_diff_tests(data, Min_col, comparison = "lower")
  }
  
  if (!is.null(Max_col) && Max_col %in% names(data)) {
    results[["max_diff"]] <- run_diff_tests(data, Max_col, comparison = "higher")
  }
  
  if (!is.null(Mean_col) && Mean_col %in% names(data)) {
    results[["mean_diff"]] <- run_diff_tests(data, Mean_col, comparison = "higher")
  }
  
  # Combine any non-empty results
  if (length(results) == 0) {
    warning("No valid columns provided. Nothing to compute.")
    return(invisible(NULL))
  }
  
  # Combine all results into one table
  test_results <- bind_rows(results) %>%
    arrange(DATE, variable)
  
  # Save to global environment using the desired name
  assign(table_name, test_results, envir = .GlobalEnv)
}



# Monthly Temperature
MeanTemp_summary <- file %>%
  filter(
    !is.na(MonthlyMinimumTemperature),
    !is.na(MonthlyMaximumTemperature),
    !is.na(MonthlyMeanTemperature)
  ) %>%
  group_by(STATION, DATE) %>%
  mutate(
    Min = MonthlyMinimumTemperature,
    Max = MonthlyMaximumTemperature,
    Mean = MonthlyMeanTemperature,
  ) %>%
  ungroup() %>%
  select(STATION,DATE,Min, Max, Mean) %>%  
  mutate(DATE = as.Date(DATE),
         Min  = as.numeric(gsub("T|s", "", Min)),
         Max  = as.numeric(gsub("T|s", "", Max)),
         Mean = as.numeric(gsub("T|s", "", Mean))
  ) %>%
  arrange(DATE, STATION)

difference_table_function(
  data = MeanTemp_summary,
  Min_col = "Min",
  Max_col = "Max",
  Mean_col = "Mean",
  table_name = "MonthlyTemp_Summary"
)


# Monthly Precipitation
MeanPrecip_summary <- file %>%
  filter(
    !is.na(MonthlyTotalLiquidPrecipitation)
  ) %>%
  group_by(STATION, DATE) %>%
  mutate(
    Mean = MonthlyTotalLiquidPrecipitation,
  ) %>%
  ungroup() %>%
  select(STATION,DATE,Mean) %>%  
  mutate(DATE = as.Date(DATE),
         Mean = as.numeric(gsub("T|s", "", Mean))
  ) %>%
  arrange(DATE, STATION)

difference_table_function(
  data = MeanPrecip_summary,
  Mean_col = "Mean",
  table_name = "MonthlyPrecip_summary"
)


## Function for t-test between sites in each month
ttest_table_function <- function(data, value_col, table_name) {
  
  run_monthly_tests <- function(data, value_col) {
    data %>%
      mutate(
        YearMonth = floor_date(as.Date(DATE), "month") # set every day (row) within a month to the same value (first day of the month)
      ) %>%
      select(YearMonth, STATION, value = !!sym(value_col)) %>% # rename the value_col to a generic value colname
      filter(!is.na(value)) %>%
      group_by(YearMonth) %>% # group by month, all days within a month group together
      filter(n_distinct(STATION) == 2) %>% # only keeps rows where values for a day exist in both sites 
      summarise(
        p_value = tryCatch(
          t.test(value ~ STATION)$p.value, # t-test comparing sites each month since data is currently grouped by month 
          error = function(e) NA
        ),
        mean_diff = tryCatch(
          diff(tapply(value, STATION, mean)), # average difference between the sites
          error = function(e) NA
        ),
        effect_size = tryCatch({
          vals <- split(value, STATION) # splits by station 
          m1 <- mean(vals[[1]]) # mean for station 1
          m2 <- mean(vals[[2]]) # mean for station 2
          s1 <- sd(vals[[1]]) # sd for station 1
          s2 <- sd(vals[[2]]) # sd for station 2
          n1 <- length(vals[[1]]) # sample size for station 1
          n2 <- length(vals[[2]]) # sample size for station 2
          
          # Pooled standard deviation
          # denominator in Cohen’s d
          spooled <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)) # combines the variability of both groups into one standard deviation measure, weighting by degrees of freedom
          # Cohen's d effect size
          (m1 - m2) / spooled # Difference between means, divided by pooled standard deviation
        }, error = function(e) NA),
        .groups = "drop"
      )
  }
  
  test_results <- run_monthly_tests(data, value_col) %>%
    arrange(YearMonth)
  
  assign(table_name, test_results, envir = .GlobalEnv)
}



# Monthly Precipitation (using daily data points)
DailyPrecip_summary <- file %>%
  mutate(
    DATE = as.Date(DATE),
    DailyPrecipitation = as.numeric(gsub("T|s", "", DailyPrecipitation))  # Convert and clean
  ) %>%
  filter(!is.na(DailyPrecipitation)) %>%  # Remove all NA (including converted "T", "s")
  select(STATION, DATE, DailyPrecipitation) %>%
  arrange(DATE, STATION)

ttest_table_function(
  data = DailyPrecip_summary,
  value_col = "DailyPrecipitation",
  table_name = "MonthlyPrecip_ttest"
)

# Monthly Temperature (using DailyAverageDryBulbTemperature)

DailyAvTemp_summary <- file %>%
  mutate(
    DATE = as.Date(DATE),
    DailyAvTemp = as.numeric(gsub("T|s", "", DailyAverageDryBulbTemperature))  # Convert and clean
  ) %>%
  filter(!is.na(DailyAvTemp)) %>%  # Remove all NA (including converted "T", "s")
  select(STATION, DATE, DailyAvTemp) %>%
  arrange(DATE, STATION)

ttest_table_function(
  data = DailyAvTemp_summary,
  value_col = "DailyAvTemp",
  table_name = "MonthlyAvTemp_ttest"
)



# Is there differences between years
# precipitation
DailyPrecip_summary2 <- DailyPrecip_summary %>% 
  mutate(
    YEAR = substr(DATE, 1, 4),
    YEARMONTH = substr(DATE, 1, 7)) %>%
  group_by(YEAR, YEARMONTH, STATION) %>%
  summarise(mean = mean(DailyPrecipitation))%>%
  ungroup()

# repeated measure anova
library(lmerTest)
library(lme4)
library(emmeans)
library(kableExtra)
anova.model <- anova(lmer(mean ~ YEARMONTH + (1|STATION),data = DailyPrecip_summary2))
anova.model


# https://www.statology.org/repeated-measures-anova-in-r/
aov_model <- aov(mean ~ factor(YEARMONTH) + Error(factor(STATION)), data = DailyPrecip_summary2)
summary(aov_model)

# https://www.r-bloggers.com/2025/02/one-way-repeated-measure-anova-in-r/
model <- aov(formula = mean ~ YEARMONTH + Error(STATION/YEARMONTH), data = DailyPrecip_summary2)
summary(model)
pttest <- pairwise.t.test(x=DailyPrecip_summary2$mean,
                          g = DailyPrecip_summary2$YEARMONTH,
                          p.adjust.method = 'bonferroni')

pub_table <- as.data.frame(as.table(pttest$p.value)) %>%
  drop_na() %>%
  rename(
    Group_1 = Var1,
    Group_2 = Var2,
    P_value = Freq
  ) %>%
  mutate(
    P_value = round(P_value, 4),
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )%>%filter(P_value < 0.05,
             str_sub(Group_1, -2, -1) == str_sub(Group_2, -2, -1)
  )

library(gt)

pub_table %>%
  gt() %>%
  tab_header(
    title = "Pairwise t-test Results",
    subtitle = "Bonferroni-adjusted p-values"
  ) %>%
  cols_label(
    Group_1 = "Group 1",
    Group_2 = "Group 2",
    P_value = "Adjusted p-value",
    Significance = "Significance"
  ) %>%
  fmt_number(
    columns = P_value,
    decimals = 4
  ) %>%
  tab_source_note(
    source_note = "Significance codes: *** p < 0.001, ** p < 0.01, * p < 0.05"
  )


# temperature 
DailyTemp_summary2 <- DailyAvTemp_summary %>% 
  mutate(
    YEAR = substr(DATE, 1, 4),
    YEARMONTH = substr(DATE, 1, 7)) %>%
  group_by(YEAR, YEARMONTH, STATION) %>%
  summarise(mean = mean(DailyAvTemp))%>%
  ungroup()

model2 <- aov(formula = mean ~ YEARMONTH + Error(STATION/YEARMONTH), data = DailyTemp_summary2)
summary(model2)
pttest2 <- pairwise.t.test(x=DailyTemp_summary2$mean,
                           g = DailyTemp_summary2$YEARMONTH,
                           p.adjust.method = 'bonferroni')

pub_table2 <- as.data.frame(as.table(pttest2$p.value)) %>%
  drop_na() %>%
  rename(
    Group_1 = Var1,
    Group_2 = Var2,
    P_value = Freq
  ) %>%
  mutate(
    P_value = round(P_value, 4),
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )%>%filter(P_value < 0.05,
             str_sub(Group_1, -2, -1) == str_sub(Group_2, -2, -1)
  )%>%
  arrange(Group_1)

pub_table2 %>%
  gt() %>%
  tab_header(
    title = "Pairwise t-test Results",
    subtitle = "Bonferroni-adjusted p-values"
  ) %>%
  cols_label(
    Group_1 = "Group 1",
    Group_2 = "Group 2",
    P_value = "Adjusted p-value",
    Significance = "Significance"
  ) %>%
  fmt_number(
    columns = P_value,
    decimals = 4
  ) %>%
  tab_source_note(
    source_note = "Significance codes: *** p < 0.001, ** p < 0.01, * p < 0.05"
  )




# Sig differences in months of interest (comparing June, July, Aug to same month in each year)
#Precipitation 
#2020-08 2019-08 : 2019 > 2020
#2021-07 2020-07 : 2020 > 2021

#Temperature
#2019-06 2018-06 : 2018 > 2019

interest <- DailyPrecip_summary2[
  DailyPrecip_summary2$YEARMONTH %in% c("2020-08", "2019-08"),
]
arrange(interest,YEARMONTH)

interest <- DailyPrecip_summary2[
  DailyPrecip_summary2$YEARMONTH %in% c("2021-07", "2020-07"),
]
arrange(interest,YEARMONTH)

interest <- DailyTemp_summary2[
  DailyTemp_summary2$YEARMONTH %in% c("2019-06","2018-06"),
]
arrange(interest,YEARMONTH)




# Sig differences in June-Aug
DailyPrecip_summary3 <- DailyPrecip_summary2 %>% 
  filter(str_sub(YEARMONTH, -2, -1) %in% c("06", "07", "08"))%>%
  group_by(YEAR, STATION) %>%
  summarise(mean = mean(mean))%>%
  ungroup()

model3 <- aov(formula = mean ~ YEAR + Error(STATION/YEAR), data = DailyPrecip_summary3)
summary(model3)
# not significant (Pr(>F) = 0.0858)


DailyTemp_summary3 <- DailyTemp_summary2 %>% 
  filter(str_sub(YEARMONTH, -2, -1) %in% c("06", "07", "08"))%>%
  group_by(YEAR, STATION) %>%
  summarise(mean = mean(mean))%>%
  ungroup()

model4 <- aov(formula = mean ~ YEAR + Error(STATION/YEAR), data = DailyTemp_summary3)
summary(model4)
pttest4 <- pairwise.t.test(x=DailyTemp_summary3$mean,
                           g = DailyTemp_summary3$YEAR,
                           p.adjust.method = 'bonferroni')

emm <- emmeans(model4, list(pairwise ~ YEAR), adjust = "tukey")

pub_table4 <- as.data.frame(emm[["pairwise differences of YEAR"]]) %>%
  filter(p.value < 0.05) %>%
  arrange(p.value)
pub_table4

# 2018- 2019 : 2018 > 2019
# 2019 - 2022 : 2022 > 2019
# 2019 - 2021 : 2021 > 2019

interest <- DailyTemp_summary3[
  DailyTemp_summary3$YEAR %in% c("2018","2019"),
]
arrange(interest,YEAR)

interest <- DailyTemp_summary3[
  DailyTemp_summary3$YEAR %in% c("2019","2022"),
]
arrange(interest,YEAR)

interest <- DailyTemp_summary3[
  DailyTemp_summary3$YEAR %in% c("2019","2021"),
]
arrange(interest,YEAR)




# Sig differences in June-July
DailyPrecip_summary4 <- DailyPrecip_summary2 %>% 
  filter(str_sub(YEARMONTH, -2, -1) %in% c("06", "07"))%>%
  group_by(YEAR, STATION) %>%
  summarise(mean = mean(mean))%>%
  ungroup()

model5 <- aov(formula = mean ~ YEAR + Error(STATION/YEAR), data = DailyPrecip_summary4)
summary(model5)
# not significant (Pr(>F) = 0.0836)

DailyTemp_summary4 <- DailyTemp_summary2 %>% 
  filter(str_sub(YEARMONTH, -2, -1) %in% c("06", "07"))%>%
  group_by(YEAR, STATION) %>%
  summarise(mean = mean(mean))%>%
  ungroup()

model6 <- aov(formula = mean ~ YEAR + Error(STATION/YEAR), data = DailyTemp_summary4)
summary(model6)
pttest5 <- pairwise.t.test(x=DailyTemp_summary4$mean,
                           g = DailyTemp_summary4$YEAR,
                           p.adjust.method = 'bonferroni')

emm <- emmeans(model6, list(pairwise ~ YEAR), adjust = "tukey")

pub_table5 <- as.data.frame(emm[["pairwise differences of YEAR"]]) %>%
  filter(p.value < 0.05) %>%
  arrange(p.value)
pub_table5

# 2018 - 2019 : 2018 > 2019
# 2019 - 2020 : 2020 > 2019
# 2019 - 2022 : 2022 > 2019
# 2018 - 2021 : 2018 > 2021
# 2019 - 2021 : 2021 > 2019
# 2018 - 2022 : 2018 > 2022
# 2018 - 2020 : 2018 > 2020

interest <- DailyTemp_summary4[
  DailyTemp_summary4$YEAR %in% c("2018","2019"),
]
arrange(interest,YEAR)

interest <- DailyTemp_summary4[
  DailyTemp_summary4$YEAR %in% c("2019","2020"),
]
arrange(interest,YEAR)

interest <- DailyTemp_summary4[
  DailyTemp_summary4$YEAR %in% c("2019","2022"),
]
arrange(interest,YEAR)

interest <- DailyTemp_summary4[
  DailyTemp_summary4$YEAR %in% c("2018","2021"),
]
arrange(interest,YEAR)

interest <- DailyTemp_summary4[
  DailyTemp_summary4$YEAR %in% c("2019","2021"),
]
arrange(interest,YEAR)

interest <- DailyTemp_summary4[
  DailyTemp_summary4$YEAR %in% c("2018","2022"),
]
arrange(interest,YEAR)

interest <- DailyTemp_summary4[
  DailyTemp_summary4$YEAR %in% c("2018","2020"),
]
arrange(interest,YEAR)

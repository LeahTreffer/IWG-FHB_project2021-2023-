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


library(dplyr)
library(fastDummies)

prep_data <- function(data, scale=FALSE, ...){
  # first want to drop all 0 electricity columns
  # data <- filter(data, bldgarea > 0 & electricity_mwh > 0)
  
  # now projecting it into sq meters
  data$bldgarea <- data$bldgarea * 0.092903
  
  drop_cols = c(
    "...1",
    "bbl",
    "bin",
    "Property Name",
    "Month",
    "cnstrct_yr",
    "coun_dist",
    "index_right"
  )
  data$month <- as.factor(data$month)
  data$year <- as.factor(data$year)
  
  data_dropped <- data[,!(names(data) %in% drop_cols)]
  
  # now want to normalize for the volume
  building_volume <- data_dropped$heightroof * data_dropped$area
  building_area <- data_dropped$bldgarea
  electric_daily <- data_dropped$electricity_mwh / data_dropped$month_days
  gas_daily <- data_dropped$gas_mwh / data_dropped$month_days
  
  # so here we can swap in and out the volume and the area
  electric_norm <- log(electric_daily / building_area)
  gas_norm <- log(gas_daily / building_area)
  
  major_terms_drop = c(
    "heightroof",
    "area",
    "month_days",
    "electricity_mwh",
    "gas_mwh",
    "bldgarea"
  )
  
  data_terms <- data_dropped[,!(names(data_dropped) %in% major_terms_drop)]
  
  # now into processing the numeric data
  data_numeric <- select_if(data_terms, is.numeric)
  
  # if (hasArg(data_mean) && hasArg(data_std)){
  #   data_numeric <- (data_numeric - data_mean) / data_std
  # } else {
  #   # here we just want to standardize the variables to avoid multicollinearity issues
  #   data_mean <- apply(data_numeric, 2, mean)
  #   # data_std <- apply(data_numeric, 2, sd)
  #   data_numeric <- (data_numeric - data_mean) # / data_std
  # }
  
  # now want to pull out the mean terms to reduce multicolinearity issues again but maintain effect sizes
  data_numeric <- scale(data_numeric, scale=scale) # / data_std

  ## res <- cor(data_numeric)
  ## correlated_idx <- findCorrelation(res, cutoff = 0.99)
  ## noncorr <- data_numeric[,-c(correlated_idx)]
  
  # now thinking about the non numeric
  data_nonnumeric <- data_terms %>% select_if(~!is.numeric(.x))
  data_num_processed <- cbind(data_numeric, data_nonnumeric)
  data_prepped <- dummy_cols(data_num_processed)
  
  returning_values <- list(electric_norm, gas_norm, data_numeric)
  return(returning_values)
}
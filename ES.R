packages <- c(
  "readxl", "dplyr", "tidyverse", "tableone", "epitools", "gtsummary", "dsr","survival","survminer")

lapply(packages, require, character.only = TRUE)

ES <- read.csv("D:/Sarcoma/Data/EwingSarcoma_DATA_LABELS_2026-02-10_1703.csv") 

age_levels <- c(
  "0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
  "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79",
  "80-84","85-89","90-94")

ES_data <- ES %>%
  mutate(Presentation.age = as.numeric(Presentation.age),
         Presentation.year = as.numeric(Presentation.year),
        year_period = cut(as.integer(Presentation.year),
                           c(1970,1981,1992,2003,2014,2025),
                           right = FALSE,
                           labels = c("1970-1980","1981-1991","1992-2002","2003-2013","2014-2024")),
         Rurality=case_when(GCH %in% c("U1","U2")~"Urban",GCH %in% c("R1","R2")~"Rural",TRUE~"Unknown"),
         Ethnicity1=case_when(Ethnicity.category %in% c("MELAA","Missing") ~"Other/Unknown",TRUE~Ethnicity.category),
         Location = case_when(
           Extraskeletal=="Yes" ~ "Extraskeletal",
           Extraskeletal!="Yes" & Site %in% c("Femur","Tibia","Fibula","Foot","Thigh","Humerus","Scapula","Clavicle") ~ "Appendicular",
           Extraskeletal!="Yes" & Site %in% c("Pelvis","Ilium","Sacrum",
                       "Cervical spine","Thoracic spine","Lumbar spine","Chest wall","Rib") ~ "Axial",
          Site == "Missing" ~ "Missing",TRUE ~ "Other"),
        age_group = case_when(
          Presentation.age >= 0 & Presentation.age <= 90 ~ as.character(cut(
            Presentation.age,
            breaks = seq(0, 90, by = 5), 
            right = FALSE,
            labels = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
              "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89"))),
          TRUE ~ "Unknown"),
        age_group = factor(age_group, levels = c(age_levels, "Unknown")))

#Table 1
Table1 <- ES_data %>%
  select(Presentation.age,Gender,Ethnicity1,year_period,Rurality,Laterality,Location,Extraskeletal,Metastasis.at.diagnosis,
         Surgery,Chemotherapy,Radiotherapy,FISH.for.EWSR1) %>%
  mutate(Ethnicity1 = fct_relevel(Ethnicity1, c("Maori", "Pacific", "Asian", "European", "Other/Unknown"))) %>% 
  tbl_summary(
    percent = "column", 
    type = list(c(Ethnicity1,Gender,year_period,Rurality,Laterality,Location,Extraskeletal,Metastasis.at.diagnosis,
                  Surgery,Chemotherapy,Radiotherapy,FISH.for.EWSR1) ~ "categorical"),
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})"),
    digits = list(all_categorical() ~ c(0, 1)))

Table1 %>% as_flex_table() %>% flextable::save_as_docx(path = "D:/Sarcoma/Result/Table1.docx")

# age distribution
ES_data %>%
  filter(!is.na(Presentation.age)) %>%
  ggplot(aes(x = Presentation.age)) +
  geom_histogram(binwidth = 5, boundary = 0, closed = "left") +
  labs(x = "Age at presentation (years)", y = "Count") +
  theme_classic()

## read NZ population data
NZ_pop <- read.csv("D:/Sarcoma/Data/nz_population_age_sex_1971_2024.csv") 

####use linear regression to predict ethnicity population each year,by accepting the values for the census year
years_to_predict <- 1970:1990
years_with_actuals <- c(
  1971, 1976, 1981,1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000,2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 
  2009, 2010,2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024)

# Create an empty data frame to store results
predicted_pop <- data.frame()

age_groups <- unique(NZ_pop$age)

# Loop over each age group
for (age in age_groups) {
    # Filter data for the current age group, total population
  age_data <- NZ_pop %>% filter(Sex == "Total", age == !!age)
    
    # Fit the linear regression model for the current age group and ethnicity
  model <- lm(Population ~ year, data = age_data, na.action = na.exclude)
    
    # Create a data frame for all years to predict
  all_years <- data.frame(year = years_to_predict)
    
    # Predict population for the specified years
  predictions <- predict(model, newdata = all_years)
    
    # Round predictions to the nearest integer
  rounded_predictions <- round(predictions)
    
    # Combine predictions with year, age group, and ethnicity
  age_predictions <- data.frame(
      year = years_to_predict,
      age = age,
      predicted_pop = rounded_predictions)
    
    # Retrieve actual values for specified years
    actual_values <- age_data %>%
      filter(year %in% years_with_actuals) %>%
      select(year, age, Population) %>%
      rename(predicted_pop = Population)  # Rename Value to match the predictions data frame
    
    # Combine actual values and predictions
    combined_results <- rbind(
      actual_values,
      age_predictions %>% filter(!year %in% years_with_actuals)  # Filter out years with actual values
    )
    # Append the results to the cumulative data frame
    predicted_pop <- rbind(predicted_pop, combined_results)}

# add WHO pop
WHO_pop <- data.frame(
  age = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
          "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89"),
  WHO_pop = c(88569,85970,84670,82171,79272,76073,71475,65877,60379,
                        86870,53681,45484,37187,29590,22092,15195,9097,4398))
# left join
predicted_pop <- predicted_pop %>%
  left_join(WHO_pop, by = "age")

write.csv(predicted_pop, file="D:/Sarcoma/Result/predicted_pop.csv",row.names=FALSE)


########Total ASR from 1974-2024##############
ES_data_count <- ES_data %>%
  #filter(Presentation.age<30)%>%
  select(Presentation.year, age_group) %>%
  group_by(Presentation.year, age_group) %>%
  summarise(count = n()) %>%
  left_join(predicted_pop, by = c("Presentation.year" = "year","age_group" = "age"))


ES_data_year_total_ASR <-  ES_data_count %>%
  filter(age_group!="Unknown")%>%
  group_by(age_group)%>%
  summarise(count_total=sum(count), Population_total=sum(predicted_pop),pop=unique(WHO_pop))%>%
  summarise(age_adjust=list(ageadjust.direct(count=count_total, pop=Population_total/1e5, stdpop=pop, rate=NULL,conf.level = 0.95))) %>%
  mutate(age_adjust = map(age_adjust, ~as.data.frame.list(.))) %>%
  unnest(cols = c(age_adjust)) 

## Total people ASR each year, 1970-2024
ES_data_year_year <- ES_data_count %>%
  filter(age_group!="Unknown")%>%
  group_by(Presentation.year) %>%
  summarise(count_total=sum(count), Population_total=sum(predicted_pop),
            age_adjust=list(ageadjust.direct(count=count, pop=predicted_pop/1e5, stdpop=WHO_pop, rate=NULL,conf.level = 0.95))) %>%
  mutate(age_adjust = map(age_adjust, ~as.data.frame.list(.))) %>%
  unnest(cols = c(age_adjust)) %>%
  mutate(SE = adj.rate / sqrt(count_total))

write.csv(ES_data_year_year, file="D:/Sarcoma/Result/ES_data_year_year.csv",row.names=FALSE)


########age group ASR from 1974-2024##############
ES_data_count <- ES_data %>%
  filter(age_group=="0-4")%>%
  select(Presentation.year, age_group) %>%
  group_by(Presentation.year, age_group) %>%
  summarise(count = n()) %>%
  left_join(predicted_pop, by = c("Presentation.year" = "year","age_group" = "age"))


ES_data_year_total_ASR <-  ES_data_count %>%
  filter(age_group!="Unknown")%>%
  group_by(age_group)%>%
  summarise(count_total=sum(count), Population_total=sum(predicted_pop),pop=unique(WHO_pop))%>%
  summarise(age_adjust=list(ageadjust.direct(count=count_total, pop=Population_total/1e5, stdpop=pop, rate=NULL,conf.level = 0.95))) %>%
  mutate(age_adjust = map(age_adjust, ~as.data.frame.list(.))) %>%
  unnest(cols = c(age_adjust)) 

## age group ASR each period, 1970-2024
ES_data_count <- ES_data %>%
  filter(Presentation.age<30) %>%
  count(year_period, Presentation.year, age_group, name = "count") %>%
  left_join(predicted_pop, by = c("Presentation.year" = "year", "age_group" = "age"))

age_levels <- c("0-4","5-9","10-14","15-19","20-24","25-29")

ES_rate_by_period <- ES_data_count %>%
  mutate(age_group = factor(age_group, levels = age_levels))%>%
  group_by(year_period, age_group) %>%
  summarise(
    count_total = sum(count, na.rm = TRUE),
    Population_total = sum(predicted_pop, na.rm = TRUE),
    crude_rate_per100k = (count_total / Population_total) * 1e5,
    .groups = "drop")


ES_rate_by_period_age <- ggplot(
  ES_rate_by_period,
  aes(x = year_period, y = crude_rate_per100k, group = age_group, color = age_group)) +
  geom_line(linewidth=2) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, 0.2)) +
  labs(
    x = "Year period",
    y = "Crude incidence rate per 100,000",
    color = "Age group") +
  theme_classic()

##### ASR by EWSR1 testing from 1974-2024##############
ES_data_EWSR1_count <- ES_data %>%
  #filter(Presentation.age<30)%>%
  select(Presentation.year, age_group,FISH.for.EWSR1) %>%
  group_by(Presentation.year, age_group,FISH.for.EWSR1) %>%
  summarise(count = n()) %>%
  left_join(predicted_pop, by = c("Presentation.year" = "year","age_group" = "age"))


ES_data_EWSR1_total_ASR <-  ES_data_EWSR1_count %>%
  filter(FISH.for.EWSR1!="Missing")%>%
  group_by(FISH.for.EWSR1,age_group)%>%
  summarise(count_total=sum(count), Population_total=sum(predicted_pop),pop=unique(WHO_pop))%>%
  summarise(age_adjust=list(ageadjust.direct(count=count_total, pop=Population_total/1e5, stdpop=pop, rate=NULL,conf.level = 0.95))) %>%
  mutate(age_adjust = map(age_adjust, ~as.data.frame.list(.))) %>%
  unnest(cols = c(age_adjust)) 


## read NZ European population data
NZ_eth_pop <- read.csv("D:/Sarcoma/Data/Estimated Resident Population by Ethinicity.csv")%>%
  filter(Ethnic.group=="European or Other (including New Zealander)", Sex=="Total people",
         !age %in% c("Total people, age","80 Years and over","65 Years and over","95 Years and over","90  Years and over",
                     "0-14 Years","15-39 Years","40-64 Years"))%>%
  mutate(
    age = trimws(age) %>%
      gsub("\\s+Years$", "", .))

####use linear regression to predict ethnicity population each year,by accepting the values for the census year
years_with_actuals <- c(2018, 2013, 2006, 2001, 1996, 1991, 1986, 1981, 2023)
years_to_predict <- setdiff(1970:2024, years_with_actuals)

# Create an empty data frame to store results
predicted_pop <- data.frame()

age_groups <- unique(NZ_eth_pop$age)

# Loop over each age group
for (age in age_groups) {
  # Filter data for the current age group, total population
  age_data <- NZ_eth_pop %>% filter(age == !!age)
  
  # Fit the linear regression model for the current age group and ethnicity
  model <- lm(Population ~ Year, data = age_data, na.action = na.exclude)
  
  # Create a data frame for all years to predict
  all_years <- data.frame(Year = years_to_predict)
  
  # Predict population for the specified years
  predictions <- predict(model, newdata = all_years)
  
  # Round predictions to the nearest integer
  rounded_predictions <- round(predictions)
  
  # Combine predictions with year, age group, and ethnicity
  age_predictions <- data.frame(
    Year = years_to_predict,
    age = age,
    predicted_pop = rounded_predictions)
  
  # Retrieve actual values for specified years
  actual_values <- age_data %>%
    filter(Year %in% years_with_actuals) %>%
    select(Year, age, Population) %>%
    rename(predicted_pop = Population)  # Rename Value to match the predictions data frame
  
  # Combine actual values and predictions
  combined_results <- rbind(
    actual_values,
    age_predictions %>% filter(!Year %in% years_with_actuals)  # Filter out years with actual values
  )
  # Append the results to the cumulative data frame
  predicted_pop <- rbind(predicted_pop, combined_results)}

# add WHO pop
WHO_pop <- data.frame(
  age = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
          "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89"),
  WHO_pop = c(88569,85970,84670,82171,79272,76073,71475,65877,60379,
              86870,53681,45484,37187,29590,22092,15195,9097,4398))
# left join
predicted_EU_pop <- predicted_pop %>%
  left_join(WHO_pop, by = "age")

write.csv(predicted_EU_pop, file="D:/Sarcoma/Result/predicted_EU_pop.csv",row.names=FALSE)

########Total European ASR from 1970-2024##############
ES_data_count <- ES_data %>%
  filter(Ethnicity1=="European", FISH.for.EWSR1=="Positive")%>%
  select(Presentation.year, age_group) %>%
  group_by(Presentation.year, age_group) %>%
  summarise(count = n()) %>%
  left_join(predicted_EU_pop, by = c("Presentation.year" = "Year","age_group" = "age"))


ES_data_year_total_ASR <-  ES_data_count %>%
  filter(age_group!="Unknown")%>%
  group_by(age_group)%>%
  summarise(count_total=sum(count), Population_total=sum(predicted_pop),pop=unique(WHO_pop))%>%
  summarise(age_adjust=list(ageadjust.direct(count=count_total, pop=Population_total/1e5, stdpop=pop, rate=NULL,conf.level = 0.95))) %>%
  mutate(age_adjust = map(age_adjust, ~as.data.frame.list(.))) %>%
  unnest(cols = c(age_adjust)) 

##ASR each year, 1970-2024
ES_data_year_year <- ES_data_count %>%
  filter(age_group!="Unknown")%>%
  group_by(Presentation.year) %>%
  summarise(count_total=sum(count), Population_total=sum(predicted_pop),
            age_adjust=list(ageadjust.direct(count=count, pop=predicted_pop/1e5, stdpop=WHO_pop, rate=NULL,conf.level = 0.95))) %>%
  mutate(age_adjust = map(age_adjust, ~as.data.frame.list(.))) %>%
  unnest(cols = c(age_adjust)) %>%
  mutate(SE = adj.rate / sqrt(count_total))

write.csv(ES_data_year_year, file="D:/Sarcoma/Result/ES_data_year_year.csv",row.names=FALSE)


## read NZ Maori population data
NZ_eth_pop <- read.csv("D:/Sarcoma/Data/Estimated Resident Population by Ethinicity.csv")%>%
  filter(Ethnic.group=="Maori", Sex=="Total people",
         !age %in% c("Total people, age","80 Years and over","65 Years and over","95 Years and over","90  Years and over",
                     "0-14 Years","15-39 Years","40-64 Years"))%>%
  mutate(
    age = trimws(age) %>%
      gsub("\\s+Years$", "", .))

####use linear regression to predict ethnicity population each year,by accepting the values for the census year
years_with_actuals <- c(2018, 2013, 2006, 2001, 1996, 1991, 1986, 1973, 1970, 2023)
years_to_predict <- setdiff(1970:2024, years_with_actuals)

# Create an empty data frame to store results
predicted_pop <- data.frame()

age_groups <- unique(NZ_eth_pop$age)

# Loop over each age group
for (age in age_groups) {
  # Filter data for the current age group, total population
  age_data <- NZ_eth_pop %>% filter(age == !!age)
  
  # Fit the linear regression model for the current age group and ethnicity
  model <- lm(Population ~ Year, data = age_data, na.action = na.exclude)
  
  # Create a data frame for all years to predict
  all_years <- data.frame(Year = years_to_predict)
  
  # Predict population for the specified years
  predictions <- predict(model, newdata = all_years)
  
  # Round predictions to the nearest integer
  rounded_predictions <- round(predictions)
  
  # Combine predictions with year, age group, and ethnicity
  age_predictions <- data.frame(
    Year = years_to_predict,
    age = age,
    predicted_pop = rounded_predictions)
  
  # Retrieve actual values for specified years
  actual_values <- age_data %>%
    filter(Year %in% years_with_actuals) %>%
    select(Year, age, Population) %>%
    rename(predicted_pop = Population)  # Rename Value to match the predictions data frame
  
  # Combine actual values and predictions
  combined_results <- rbind(
    actual_values,
    age_predictions %>% filter(!Year %in% years_with_actuals)  # Filter out years with actual values
  )
  # Append the results to the cumulative data frame
  predicted_pop <- rbind(predicted_pop, combined_results)}

# add WHO pop
WHO_pop <- data.frame(
  age = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
          "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89"),
  WHO_pop = c(88569,85970,84670,82171,79272,76073,71475,65877,60379,
              86870,53681,45484,37187,29590,22092,15195,9097,4398))
# left join
predicted_maori_pop <- predicted_pop %>%
  left_join(WHO_pop, by = "age")

write.csv(predicted_maori_pop, file="D:/Sarcoma/Result/predicted_maori_pop.csv",row.names=FALSE)

########Total ASR from 1974-2024##############
ES_data_count <- ES_data %>%
  filter(Ethnicity1=="Maori",FISH.for.EWSR1=="Positive")%>%
  select(Presentation.year, age_group) %>%
  group_by(Presentation.year, age_group) %>%
  summarise(count = n()) %>%
  left_join(predicted_maori_pop, by = c("Presentation.year" = "Year","age_group" = "age"))


ES_data_year_total_ASR <-  ES_data_count %>%
  filter(age_group!="Unknown")%>%
  group_by(age_group)%>%
  summarise(count_total=sum(count), Population_total=sum(predicted_pop),pop=unique(WHO_pop))%>%
  summarise(age_adjust=list(ageadjust.direct(count=count_total, pop=Population_total/1e5, stdpop=pop, rate=NULL,conf.level = 0.95))) %>%
  mutate(age_adjust = map(age_adjust, ~as.data.frame.list(.))) %>%
  unnest(cols = c(age_adjust)) 

##ASR each year, 1970-2024
ES_data_year_year <- ES_data_count %>%
  filter(age_group!="Unknown")%>%
  group_by(Presentation.year) %>%
  summarise(count_total=sum(count), Population_total=sum(predicted_pop),
            age_adjust=list(ageadjust.direct(count=count, pop=predicted_pop/1e5, stdpop=WHO_pop, rate=NULL,conf.level = 0.95))) %>%
  mutate(age_adjust = map(age_adjust, ~as.data.frame.list(.))) %>%
  unnest(cols = c(age_adjust)) %>%
  mutate(SE = adj.rate / sqrt(count_total))

write.csv(ES_data_year_year, file="D:/Sarcoma/Result/ES_data_year_year.csv",row.names=FALSE)

##### read NZ Pacific population data #####
NZ_eth_pop <- read.csv("D:/Sarcoma/Data/Estimated Resident Population by Ethinicity.csv")%>%
  filter(Ethnic.group=="Pacific", Sex=="Total people",
         !age %in% c("Total people, age","80 Years and over","65 Years and over","95 Years and over","90  Years and over",
                     "0-14 Years","15-39 Years","40-64 Years"))%>%
  mutate(
    age = trimws(age) %>%
      gsub("\\s+Years$", "", .))

####use linear regression to predict ethnicity population each year,by accepting the values for the census year
years_with_actuals <- c(2018, 2013, 2006, 2001, 1996, 1991, 1986, 1981, 2023)
years_to_predict <- setdiff(1970:2024, years_with_actuals)

# Create an empty data frame to store results
predicted_pop <- data.frame()

age_groups <- unique(NZ_eth_pop$age)

# Loop over each age group
for (age in age_groups) {
  # Filter data for the current age group, total population
  age_data <- NZ_eth_pop %>% filter(age == !!age)
  
  # Fit the linear regression model for the current age group and ethnicity
  model <- lm(Population ~ Year, data = age_data, na.action = na.exclude)
  
  # Create a data frame for all years to predict
  all_years <- data.frame(Year = years_to_predict)
  
  # Predict population for the specified years
  predictions <- predict(model, newdata = all_years)
  
  # Round predictions to the nearest integer
  rounded_predictions <- round(predictions)
  
  # Combine predictions with year, age group, and ethnicity
  age_predictions <- data.frame(
    Year = years_to_predict,
    age = age,
    predicted_pop = rounded_predictions)
  
  # Retrieve actual values for specified years
  actual_values <- age_data %>%
    filter(Year %in% years_with_actuals) %>%
    select(Year, age, Population) %>%
    rename(predicted_pop = Population)  # Rename Value to match the predictions data frame
  
  # Combine actual values and predictions
  combined_results <- rbind(
    actual_values,
    age_predictions %>% filter(!Year %in% years_with_actuals)  # Filter out years with actual values
  )
  # Append the results to the cumulative data frame
  predicted_pop <- rbind(predicted_pop, combined_results)}

# add WHO pop
WHO_pop <- data.frame(
  age = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
          "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89"),
  WHO_pop = c(88569,85970,84670,82171,79272,76073,71475,65877,60379,
              86870,53681,45484,37187,29590,22092,15195,9097,4398))
# left join
predicted_Pacific_pop <- predicted_pop %>%
  left_join(WHO_pop, by = "age")

write.csv(predicted_Pacific_pop, file="D:/Sarcoma/Result/predicted_Pacific_pop.csv",row.names=FALSE)


### Total ASR from 1974-2024 ###
ES_data_Pacific_count <- ES_data %>%
  filter(Ethnicity1=="Pacific",FISH.for.EWSR1=="Positive")%>% #,FISH.for.EWSR1=="Positive"
  select(Presentation.year, age_group) %>%
  group_by(Presentation.year, age_group) %>%
  summarise(count = n()) %>%
  left_join(predicted_Pacific_pop, by = c("Presentation.year" = "Year","age_group" = "age"))


ES_data_Pacific_total_ASR <-  ES_data_Pacific_count %>%
  filter(age_group!="Unknown")%>%
  group_by(age_group)%>%
  summarise(count_total=sum(count), Population_total=sum(predicted_pop),pop=unique(WHO_pop))%>%
  summarise(age_adjust=list(ageadjust.direct(count=count_total, pop=Population_total/1e6, stdpop=pop, rate=NULL,conf.level = 0.95))) %>%
  mutate(age_adjust = map(age_adjust, ~as.data.frame.list(.))) %>%
  unnest(cols = c(age_adjust)) 

##ASR each year, 1970-2024
ES_data_year_year <- ES_data_count %>%
  filter(age_group!="Unknown")%>%
  group_by(Presentation.year) %>%
  summarise(count_total=sum(count), Population_total=sum(predicted_pop),
            age_adjust=list(ageadjust.direct(count=count, pop=predicted_pop/1e5, stdpop=WHO_pop, rate=NULL,conf.level = 0.95))) %>%
  mutate(age_adjust = map(age_adjust, ~as.data.frame.list(.))) %>%
  unnest(cols = c(age_adjust)) %>%
  mutate(SE = adj.rate / sqrt(count_total))

write.csv(ES_data_year_year, file="D:/Sarcoma/Result/ES_data_year_year.csv",row.names=FALSE)



##### survival #####
end_year <- 2025

ES_data_survival <- ES_data %>%
  filter(!is.na(Presentation.year)) %>%
  mutate(
    Decease.year = suppressWarnings(as.numeric(trimws(as.character(Decease.year)))),
    event = ifelse(Deceased == "Yes", 1, 0),
    exit_year = ifelse(event == 0, end_year, Decease.year),
    time = exit_year - Presentation.year)


ES_data_survival_rate <- survfit(Surv(time, event) ~ 1, data = ES_data_survival)
summary(ES_data_survival_rate)

# median follow-up
fit_fu <- survfit(Surv(time, 1 - event) ~ 1, data = ES_data_survival)
summary(fit_fu)$table[c("median","0.95LCL","0.95UCL")]

ggsurvplot(ES_data_survival_rate, data = ES_data_survival,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.height = 0.25,  
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE, 
           ggtheme = theme_classic(base_size = 12),
           xlab = "Years since presentation",
           ylab = "Overall survival probability")
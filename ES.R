packages <- c(
  "readxl", "dplyr", "tidyverse", "tableone", "epitools", "gtsummary", "dsr")

lapply(packages, require, character.only = TRUE)

ES_data <- read.csv("D:/Sarcoma/Data/EwingSarcoma_DATA_LABELS_2026-02-10_1703.csv") %>%
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
           Site %in% c("Femur","Tibia","Fibula","Foot","Thigh","Humerus") ~ "Appendicular",
           Site %in% c("Pelvis","Ilium","Sacrum",
                       "Cervical spine","Thoracic spine","Lumbar spine",
                       "Chest wall","Rib","Scapula","Clavicle") ~ "Axial",
          Site == "Missing" ~ "Missing",TRUE ~ "Other"),
        age_group = case_when(
          Presentation.age >= 0 & Presentation.age <= 90 ~ as.character(cut(
            Presentation.age,
            breaks = seq(0, 90, by = 5), 
            right = FALSE,
            labels = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
              "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89"))),
          TRUE ~ "Unknown"))

#Table 1
Table1 <- ES_data %>%
  select(Presentation.age,Gender,Ethnicity1,year_period,Rurality,Laterality,Location,Extraskeletal,Metastasis.at.diagnosis,
         Surgery,Chemotherapy,Radiotherapy,FISH.for.EWSR1) %>%
  mutate(Ethnicity1 = fct_relevel(Ethnicity1, c("Maori", "Pacific", "Asian", "European", "Other/Unknown"))) %>% 
  tbl_summary(
    percent = "column", 
    type = list(c(Ethnicity1,Gender,year_period,Rurality,Laterality,Location,Extraskeletal,Metastasis.at.diagnosis,
                  Surgery,Chemotherapy,Radiotherapy,FISH.for.EWSR1) ~ "categorical"),
    statistic = list(all_continuous() ~ "{median} ({min},{max})"),
    digits = list(all_categorical() ~ c(0, 1)))

Table1 %>% as_flex_table() %>% flextable::save_as_docx(path = "D:/Sarcoma/Result/Table1.docx")

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

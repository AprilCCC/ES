packages <- c(
  "readxl", "dplyr", "tidyverse", "tableone", "epitools", "gtsummary", "dsr")

lapply(packages, require, character.only = TRUE)

ES_data <- read.csv("D:/Sarcoma/Data/EwingSarcoma_DATA_LABELS_2026-02-10_1703.csv") %>%
  mutate(Presentation.age=as.numeric(Presentation.age),
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
          Site == "Missing" ~ "Missing",TRUE ~ "Other")) 

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


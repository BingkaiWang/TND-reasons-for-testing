library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(data.table)
library(ggplot2)
library(stringr)

# Load Immunization Data
end_date <- as.Date("2023-10-01")
if (end_date < as.Date("2022-08-01"))
{
  immunization_data <- read.csv("S:/DataDirect/HUM00164771 - Study associations between vaccines/LiyangYUAN/Flu_study/immunization/Immunizations_clean_0802.csv")
} else
{
  immunization_data <- read.csv("S:/DataDirect/HUM00164771 - Study associations between vaccines/Cleaned_Data_10032023/Immunizations_clean.csv")
}


# Load PCR Test Data

if (end_date < as.Date("2022-08-01"))
{
  pcr_data <- read.csv("S:/DataDirect/HUM00164771 - Study associations between vaccines/LiyangYUAN/Flu_study/Infections/PCR_Covid_0802.csv")
} else
{
  pcr_data <- read.csv("S:/DataDirect/HUM00164771 - Study associations between vaccines/Automated_Output_backup/HPI-6693_Data_2023-10-03/HPI-6693_COVID_Tests_2023-10-03.csv")
}


pcr_data <- pcr_data %>%
  mutate(COLLECTION_DATE = strptime(COLLECTION_DATE, format = "%m/%d/%Y %H:%M"))

sort_table = sort(table(pcr_data$ResultValue), decreasing =  TRUE)
ValidRes <- rownames(sort_table)[1:18]
labs <- pcr_data %>%
  filter(ResultValue %in% ValidRes)

## A lot of the responses of ResultValue in the PCR Data Set is not applicable, so we use the top 18, which yield 1044030 readings (99.5%) the amount of total data set.  

Rapid <- names(table(labs$OrderTestName))[grep(c("rapid"), tolower(names(table(labs$OrderTestName))))]
PCR <- names(table(labs$OrderTestName))[grep(c("pcr"), tolower(names(table(labs$OrderTestName))))]
all_Testname <- names(table(labs$OrderTestName))
PCR_Name = unique(c(
  #Regular PCR
  '2019 Novel Coronavirus Real-Time PCR, NP/OP', 
  'SARS Coronavirus with CoV-2 RNA, QL, Real-Time RT-PCR (UHS)', 
  'Novel Coronavirus (COVID-19), PCR',
  'ER SARS-CoV-2 (COVID-19) by PCR',
  'SARS-CoV-2 (COVID-19) by PCR, Health Care Worker/Employee',
  
  #Rapid PCR
  'Rapid Novel Coronavirus (COVID-19), PCR',
  'NC Rapid Novel Coronavirus(COVID-19)PCR',
  
  #Chen's new update
  "SARS-CoV-2, Flu A, Flu B and RSV, PCR, UHS",
  "SARS-CoV-2 (COVID-19), Flu A+B, and RSV by PCR",
  "POC SARS-CoV-2 (COVID-19)/ Influenza (A,B)/ RSV by PCR",
  
  Rapid,
  PCR,
  
  "POC COVID (SARS-CoV-2)"))
# We further filter COVID Tests by fitering tests with valid ResultTestName.
unk_dnr = c("UNK","DNR")

PCR_covid = labs %>%
  filter(!(grepl("First", ResultTestName, ignore.case= TRUE)))%>%
  # The question "First COVID Test?" May Interfere with Result/Symptom testing. 
  mutate(
    Result_PCR = ifelse(
      ResultValue %in% unk_dnr, "UNK_DNR", 
      (tolower(str_sub(ResultValue,1 , 8)) == 'detected') | (ResultValue == 'Presumptive Positive') | (OrderTestName %in% PCR_Name & tolower(ResultValue) == 'positive')
    )
  ) %>%
  mutate(Result_Sym = ifelse(
    ResultValue %in% unk_dnr, "UNK_DNR", 
    ((grepl("symptoms", ResultTestName, ignore.case= TRUE)&ResultValue =="Y") ) )
  ) %>%
  select(DEID_PatientID, COLLECTION_DATE, Result_PCR,Result_Sym)

Result_Aggregation <- function(res)
{
  if (any(res =="TRUE"))
    return ("TRUE")
  if (any(res =="FALSE"))
    return ("FALSE")
  return("UNK_DNR")
}

pcr_data_res <- PCR_covid %>%
  group_by(DEID_PatientID, COLLECTION_DATE) %>%
  summarise(
    Result_PCR = Result_Aggregation(Result_PCR),
    Result_Sym = Result_Aggregation(Result_Sym)
  ) %>%
  ungroup()

# The result of PCR tests are classified into 3 categories: "TRUE" if at least 1 result of the day is "Positive" or equivalent, "FALSE" if all results with reading of the day are "Negative" or equivalent, and "UNK/DNR" if all the results are without reading or with reading UNK/DNR. 

##We aggregate the results of tests, such that for all tests done in a day, the result is "TRUE" if at least one of the results in that day is related to "detected", "FALSE" if not all remaining results are UNK_DNR (Unknown or Did Not Report) and UNK_DNR if all remaining results are UNK_DNR. Next, we use the  immunization data by October 2023, where the variable Valid is 1 if both of the following cases are true,then value will be 1, otherwise = 0;

# concatenated vaccine name is reasonable, which means one of case in (Janssen, Pfizer, Pfizer_Pfizer, Moderna,  Moderna_Moderna);

# If date of the 1st vaccine is earlier than date of  the 2nd vaccine;

# Most Recent Vaccination level taken on or before test data (1st,2nd,1st booster, 2nd booster)- reflecting immunization levels. 
# Days since the most recent vaccination


combined_data_uniq <- pcr_data_res %>%
  left_join(immunization_data, by = c("DEID_PatientID" = "PatID"))

combined_data_uniq <- combined_data_uniq %>%
  mutate(across(starts_with("DateVax"),as.Date, format = "%Y-%m-%d"))

combined_data_long <- combined_data_uniq %>%
  mutate(
    Vax_Level = case_when(
      !is.na(DateVax4) & DateVax4 <= COLLECTION_DATE ~ 4,
      !is.na(DateVax3) & DateVax3 <= COLLECTION_DATE ~ 3,
      !is.na(DateVax2) & DateVax2 <= COLLECTION_DATE ~ 2,
      !is.na(DateVax1) & DateVax1 <= COLLECTION_DATE ~ 1,
      TRUE ~ 0
    ),
    Latest_Vax_Date = case_when(
      !is.na(DateVax4) & DateVax4 <= COLLECTION_DATE ~ DateVax4,
      !is.na(DateVax3) & DateVax3 <= COLLECTION_DATE ~ DateVax3,
      !is.na(DateVax2) & DateVax2 <= COLLECTION_DATE ~ DateVax2,
      !is.na(DateVax1) & DateVax1 <= COLLECTION_DATE ~ DateVax1,
      TRUE ~ NA_Date_
    ),
    Latest_Vax_Type = case_when(
      !is.na(DateVax4) & DateVax4 <= COLLECTION_DATE ~ ifelse(Dose4Type =="", NA, Dose4Type),
      !is.na(DateVax3) & DateVax3 <= COLLECTION_DATE ~ ifelse(Dose3Type =="", NA, Dose3Type),
      !is.na(DateVax2) & DateVax2 <= COLLECTION_DATE ~ ifelse(Dose2Type =="", NA, Dose2Type),
      !is.na(DateVax1) & DateVax1 <= COLLECTION_DATE ~ ifelse(Dose1Type =="", NA, Dose1Type),
      TRUE ~ NA
    ),
  )

combined_data_long <- combined_data_long %>%
  mutate(Days_Since_Vax = ifelse(is.na(Latest_Vax_Date ), NA, as.numeric( difftime(COLLECTION_DATE , as.POSIXct(Latest_Vax_Date), units = "days") )))

combined_data_long_reduced <- combined_data_long %>%
  select(DEID_PatientID,COLLECTION_DATE,Result_PCR,Result_Sym, Vaxname, DateVax1,DateVax2,DateVax3,DateVax4, Dose1Type,Dose2Type,Dose3Type,Dose4Type, Vax_Level,Latest_Vax_Date,Latest_Vax_Type, Days_Since_Vax)

## For each test, we collect the nearest vaccination date prior to test collection. We evaluate the time spent from most recent vaccination to testing, as well as time and type of previous doses. 
write.csv(combined_data_long_reduced , "S:/DataDirect/HUM00164771 - Study associations between vaccines/TomLiu_practice/pcr_vax.csv",row.names = FALSE)
## Now we work on ICD-10 Data. This csv file contains all ICD-10 records from Jan 2020 to Oct 2023, with each period of hospitalization identified by a deid_encounterid. 

file_paths <- c("S:/DataDirect/HUM00164771 - Study associations between vaccines/Automated_Output_backup/HPI-6693_Data_2023-10-03/HPI-6693_ICD10_Set1of4_2023-10-03.csv",
                
                "S:/DataDirect/HUM00164771 - Study associations between vaccines/Automated_Output_backup/HPI-6693_Data_2023-10-03/HPI-6693_ICD10_Set2of4_2023-10-03.csv",
                
                "S:/DataDirect/HUM00164771 - Study associations between vaccines/Automated_Output_backup/HPI-6693_Data_2023-10-03/HPI-6693_ICD10_Set3of4_2023-10-03.csv",
                
                "S:/DataDirect/HUM00164771 - Study associations between vaccines/Automated_Output_backup/HPI-6693_Data_2023-10-03/HPI-6693_ICD10_Set4of4_2023-10-03.csv"    )

data_list <- lapply(file_paths, fread)

filtered_data_list <- lapply(data_list, function(dt)
{
  dt[,admitdate:= as.POSIXct(admitdate, format = "%m/%d/%Y %H:%M")]
  dt[admitdate >= as.POSIXct("2020-01-01")]
}
)

### Filtering only Patients admitted after Jan.2020.
merged_data <- rbindlist(filtered_data_list)
head(merged_data)
encounter_data <- unique(merged_data[,.(DeID_PatientID, deid_EncounterID,admitdate, dischargedate,PatientClass, DxCode)])
demographics_data <- fread("S:/DataDirect/HUM00164771 - Study associations between vaccines/Automated_Output_backup/HPI-6693_Data_2023-10-03/HPI-6693_Demographics_2023-10-03.csv")
demographics_data[, DeceasedDate := as.POSIXct(DeceasedDate)]
(demographics_info<- colnames(demographics_data)[c(1:8,10:14)])
# Information of interest, subject to adjustment
setnames(demographics_data, "DEID_PatientID", "DeID_PatientID")
(demographics_subdata<- demographics_data[,c(1,14)])
joined_data_pcp<- merge(encounter_data,  demographics_subdata, by = "DeID_PatientID", all.x = TRUE )
joined_data_pcp<- joined_data_pcp[UM_PCP_PRECISION_HEALTH_DEF == "Y"]

## The following ICD-10 code will be of interest for our TND Design Study. 

COVID_ICD10_Code <- c("U07.1","U07.2")
symptom_ICD10_Code <- c("R50.9","R53.83")
case_tracing_ICD10_Code <- c("Z20.822")
non_symptomatic_ICD10_Code <- c("Z11.52")
ICD10_Code_ForTND <- c(symptom_ICD10_Code, case_tracing_ICD10_Code,  non_symptomatic_ICD10_Code,COVID_ICD10_Code)

joined_data_pcp[,DxCode:= ifelse(DxCode %in% ICD10_Code_ForTND, DxCode, "Other")]

joined_data_pcp <- unique(joined_data_pcp[,c(1,3:6)])
## Keep a row iff  (DxCode %in% ICD10_Code_ForTND). For a DxCode beyond our interest, we set it to "Other" for aggregation. 

fwrite(joined_data_pcp, "ICD10_data_with_pcp.csv")

## This csv file contains all ICD-10 records, and anyone in this csv has a primary care physician at UM, and thus considered as participant of our study. 

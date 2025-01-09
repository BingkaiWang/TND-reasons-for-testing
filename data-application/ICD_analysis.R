library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
library(lubridate)
library(data.table)
library(ggplot2)
library(stringr)
library(ggpubr)
library(tableone)

time_of_use = 3

pcr_vax_data <- fread("S:/DataDirect/HUM00164771 - Study associations between vaccines/TomLiu_practice/pcr_vax.csv")
joined_data_pcp <- fread("S:/DataDirect/HUM00164771 - Study associations between vaccines/TomLiu_practice/ICD10_data_with_pcp.csv")
pcr <- fread("S:/DataDirect/HUM00164771 - Study associations between vaccines/TomLiu_practice/pcr.csv")
immunization_data <- fread("S:/DataDirect/HUM00164771 - Study associations between vaccines/Cleaned_Data_10032023/Immunizations_clean.csv")
setnames(pcr_vax_data, "DEID_PatientID", "DeID_PatientID")
pcr_vax_data[,DDate := paste0(DeID_PatientID,COLLECTION_DATE)]
length(unique(joined_data_pcp[DeID_PatientID %in% pcr_vax_data$DeID_PatientID]$DeID_PatientID))
length(unique(pcr_vax_data$DeID_PatientID))


#Among 152044 Aggregated PCR Tests taken, 53224 of them can be associated with a patient possessing both a Primary Care Physicial at Precision Medicine, and an identifiable Hospitalization/DxCode of interest. 


joined_data_pcr <- merge( pcr_vax_data, joined_data_pcp , by = "DeID_PatientID", allow.cartesian = TRUE)


joined_data_pcr[, ':='(
  test_in_hp = (PatientClass !="Outpatient") & (COLLECTION_DATE >= admitdate &  COLLECTION_DATE <= dischargedate),
  test_3_prior = (PatientClass !="Outpatient") & (as.numeric(difftime(admitdate,COLLECTION_DATE,units="days"))>0 & as.numeric(difftime(admitdate,COLLECTION_DATE,units="days"))<=3)
)]

between01 <- function (x)
{
  return (x>0&x<=1)
}

joined_data_pcr[, ':='(
  R_td = ifelse( DxCode %in%c("R50.9","R53.83") & between01(abs(as.numeric(difftime(admitdate,COLLECTION_DATE,units="days")))),1,0),
  U_td = ifelse( DxCode %in%c("U07.1","U07.2") & between01(abs(as.numeric(difftime(admitdate,COLLECTION_DATE,units="days")))),1,0),
  Z20_td = ifelse(DxCode=="Z20.822" & between01(abs(as.numeric(difftime(admitdate,COLLECTION_DATE,units="days")))),1,0),
  Z11_td = ifelse(DxCode %in% c("Z11.52","Z02.89") & between01(abs(as.numeric(difftime(admitdate,COLLECTION_DATE,units="days")))),1,0)
)]

#View(joined_data_pcr)

demographics_data <- fread("S:/DataDirect/HUM00164771 - Study associations between vaccines/Automated_Output_backup/HPI-6693_Data_2023-10-03/HPI-6693_Demographics_2023-10-03.csv")

setnames(demographics_data, "DEID_PatientID","DeID_PatientID")
demographics_data[, Age.FirstDose:= ifelse(`Age@FirstDose`==">89",90,as.numeric(`Age@FirstDose`))]
demographics_data[, BMI:= ifelse(BMI >= 50,50,BMI)]
joined_data_final<- merge(joined_data_pcr,  demographics_data[,c(1,3:13)], by = "DeID_PatientID", all.x = TRUE )

#The aggregation data set is constructed such that each observation is uniquely identified by a combination of DeID_PatientID and COLLECTION_DATE of COVID tests. For each ID-Date combination, we know whether this test happens within at least 1 period of hospitalization, and that whether at least hospitalization happens within 3 days of the test. 

joined_data_final_agg <-  joined_data_final[,.(
  test_in_hp = ifelse(any(test_in_hp ==TRUE, na.rm=TRUE),1,0),
  test_3_prior = ifelse(any(test_3_prior == TRUE, na.rm = TRUE ),1,0),
  U_td = ifelse(any(U_td  ==1, na.rm = TRUE ),1,0),
  R_td = ifelse(any(R_td  ==1, na.rm = TRUE ),1,0),
  Z20_td = ifelse(any(Z20_td  ==1, na.rm = TRUE ),1,0),
  Z11_td = ifelse(any(Z11_td  ==1, na.rm = TRUE ),1,0)
),by = .(DeID_PatientID, COLLECTION_DATE)]

pcr_vax_data[,DDate := paste0(DeID_PatientID,COLLECTION_DATE)]
joined_data_final_agg[,DDate := paste0(DeID_PatientID,COLLECTION_DATE)]
pcr_vax_data[,c("DeID_PatientID","COLLECTION_DATE"):=NULL]
joined_data_final_agg<- merge(joined_data_final_agg, pcr_vax_data , by = "DDate", all.x=TRUE)
joined_data_final_agg<- merge(joined_data_final_agg,  demographics_data[,c(1,3:13)], by = "DeID_PatientID", all.x = TRUE )

joined_data_final_agg[, Age.FirstDose:= ifelse(`Age@FirstDose`==">89",90,as.numeric(`Age@FirstDose`))]
joined_data_final_agg$Omicron_Period_Flag <- with(joined_data_final_agg,
                                                  ifelse(COLLECTION_DATE >= as.POSIXct("2020-01-01") &  COLLECTION_DATE <= as.POSIXct("2021-06-30"),0,
                                                         ifelse(COLLECTION_DATE >= as.POSIXct("2021-07-01") &  COLLECTION_DATE <= as.POSIXct("2021-12-31"),1, 
                                                                ifelse(COLLECTION_DATE >= as.POSIXct("2022-01-01") &  COLLECTION_DATE <= as.POSIXct("2022-06-30"),2,
                                                                       ifelse(COLLECTION_DATE >= as.POSIXct("2022-07-01") &  COLLECTION_DATE <= as.POSIXct("2023-09-30"),3,  NA))))    
)                                              

joined_data_final_agg[joined_data_final_agg$BMI>60]$BMI = 60    
joined_data_final_agg[,TestingReason:= 9] 


# If a test happens outside Hospitalization and has a positive symptom result, we consider it a test due to Symptoms. 
joined_data_final_agg[,TestingReason:= ifelse(Result_Sym =="TRUE" & test_in_hp ==0,1,TestingReason)] 

# If a test happens within 1 day of the assignment for ICD-10 code "Z20.822" and has neither symptom nor hospitalization, we consider it a test due to case contact tracing. 
joined_data_final_agg[,TestingReason:= ifelse(Z20_td ==1,3,TestingReason)] 

joined_data_final_agg$TestingReason <- as.factor(joined_data_final_agg$TestingReason)
joined_data_final_agg_description <- joined_data_final_agg %>% 
  select(Result_PCR, Result_Sym, Vax_Level,Dose1Type, Dose2Type
         , Age.FirstDose, Gender,Race, BMI, TestingReason) %>%
  mutate(AgeCat = as.character(ifelse (Age.FirstDose <= 18, "<=18", 
                                       ifelse (Age.FirstDose >=60, ">=60", 
                                               "18~60") )))

descriptive_table <- CreateTableOne(
  vars = c("Gender", "Race", "AgeCat", "BMI", "Vax_Level", "Dose1Type","Dose2Type", "Result_PCR"),
  strata = "TestingReason",
  data = joined_data_final_agg_description,
  includeNA = TRUE,
  addOverall = FALSE
)

## Summary Statistics

common_xlim <- as.Date(c("2020-01-01","2023-10-01"))
plot_PC<- ggplot(data = subset(pcr, ( Result_PCR=="TRUE")),aes(x = floor_date(as.Date(COLLECTION_DATE), unit = "week") ))+
  geom_bar(stat = "count") +
  labs(title = "COVID-19 Tests (PCR Positive Cases)", x = "Date", y = "Count (Aggregated By Week)")+
  geom_vline(xintercept = as.Date("2021-07-01"), color = "darkgreen", linetype="dashed")+
  geom_vline(xintercept = as.Date("2021-12-31"), color = "red", linetype="dashed")+
  geom_vline(xintercept = as.Date("2022-06-30"), color = "magenta", linetype="dashed")+
  geom_vline(xintercept = as.Date("2023-09-30"), color = "orange", linetype="dashed")+
  # annotate("text", x = as.Date("2020-11-04"), y = 500, label = "First Vaccination\n Recorded",vjust = -0.5, color = "darkgreen" )+
  # annotate("text", x = as.Date("2021-09-26"), y = 500, label = "Omicron Variant\n Named",vjust = -0.5, color = "red" )+
  # annotate("text", x = as.Date("2022-05-29"), y = 400, label = "Omicron Subvariant\n Took Over",vjust = -0.5, color = "magenta" )+
  # annotate("text", x = as.Date("2023-07-31"), y = 500, label = "End of \n Data Collection",vjust = -0.5, color = "orange" )+
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y", limits = common_xlim) +
  theme_minimal()
plot_Vax<- ggplot(data = immunization_data,aes(x = floor_date(as.Date(DateVax1), unit = "week") ))+
  geom_bar(stat = "count") + 
  labs(title = "First Dose of COVID-19 Vaccination", x = "Date", y = "Count (Aggregated By Week)")+
  geom_vline(xintercept = as.Date("2021-07-01"), color = "darkgreen", linetype="dashed")+
  geom_vline(xintercept = as.Date("2021-12-31"), color = "red", linetype="dashed")+
  geom_vline(xintercept = as.Date("2022-06-30"), color = "magenta", linetype="dashed")+
  geom_vline(xintercept = as.Date("2023-09-30"), color = "orange", linetype="dashed")+
  annotate("text", x = as.Date("2020-11-04"), y = 30000, label = "Delta Variant\n Named",vjust = -0.5, color = "darkgreen" )+
  annotate("text", x = as.Date("2021-09-26"), y = 30000, label = "Omicron Variant\n Named",vjust = -0.5, color = "red" )+
  annotate("text", x = as.Date("2022-05-29"), y = 20000, label = "Omicron Subvariant\n Took Over",vjust = -0.5, color = "magenta" )+
  annotate("text", x = as.Date("2023-07-31"), y = 20000, label = "End of \n Data Collection",vjust = -0.5, color = "orange" )+
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y", limits = common_xlim) +
  theme_minimal()


demographics_data_plt <- demographics_data %>%
  filter(!is.na(Age.FirstDose),!is.na(BMI))

plot_age <- demographics_data_plt%>%
  ggplot(aes(x = Age.FirstDose))+
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black")+
  labs(title ="Age at first Dose of COVID Vaccination", x = "Age at 1st Dose")+
  theme_minimal()

demographics_data_plt <- demographics_data_plt%>%
  mutate(BMI_Category = case_when(
    BMI < 18.5 ~ "Underweight",
    BMI >= 18.5 & BMI < 25 ~ "Healthy Weight",
    BMI >= 25 & BMI < 30 ~ "Overweight",
    BMI >= 30 & BMI < 40  ~ "Obesity Class 1&2",
    BMI >= 40  ~ "Obesity Class 3"
  ))

bmi_cnt  <- demographics_data_plt %>% count(BMI_Category)

plot_BMI <- bmi_cnt %>%
  ggplot(aes(x="",y=n, fill = BMI_Category)) + 
  geom_bar(stat = "identity", width = 1) + 
  coord_polar("y", start = 0) + 
  labs(title = "BMI Categories") + 
  theme_void() + 
  theme(legend.title = element_blank())


ICD10_data <-fread("ICD10_data_with_pcp.csv")
plot_EHR<- ggplot(data = ICD10_data,aes(x = floor_date(as.Date(admitdate), unit = "week") ))+
  geom_bar(stat = "count") +
  labs(title = "Number of Admissions", x = "Date", y = "Count (Aggregated By Week)")+
  geom_vline(xintercept = as.Date("2020-11-04"), color = "darkgreen", linetype="dashed")+
  geom_vline(xintercept = as.Date("2021-11-26"), color = "red", linetype="dashed")+
  geom_vline(xintercept = as.Date("2022-03-29"), color = "magenta", linetype="dashed")+
  geom_vline(xintercept = as.Date("2023-09-30"), color = "orange", linetype="dashed")+
  # annotate("text", x = as.Date("2020-11-04"), y = 500, label = "First Vaccination\n Recorded",vjust = -0.5, color = "darkgreen" )+
  # annotate("text", x = as.Date("2021-09-26"), y = 500, label = "Omicron Variant\n Named",vjust = -0.5, color = "red" )+
  # annotate("text", x = as.Date("2022-05-29"), y = 400, label = "Omicron Subvariant\n Took Over",vjust = -0.5, color = "magenta" )+
  # annotate("text", x = as.Date("2023-07-31"), y = 500, label = "End of \n Data Collection",vjust = -0.5, color = "orange" )+
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y", limits = common_xlim) +
  theme_minimal()

PatientClass_cnt  <- ICD10_data %>% count(PatientClass)

plot_PatientClass <- PatientClass_cnt%>%
  ggplot(aes(x="",y=n, fill = PatientClass)) + 
  geom_bar(stat = "identity", width = 1) + 
  coord_polar("y", start = 0) + 
  labs(title = "PatientClass Categories") + 
  theme_void() + 
  theme(legend.title = element_blank())

## Data Analysis for Delta Period

df_deltap <- joined_data_final_agg %>%
  filter(Omicron_Period_Flag ==1)

df_delta_naive <- df_deltap%>%
  filter(TestingReason %in% c(1,3)) %>%
  mutate(d.Test = as.factor(ifelse(Result_PCR=="TRUE" & Result_Sym == "TRUE",1,0))) %>%
  mutate(V.Test = ifelse(Vax_Level>0, 1, 0))

table(df_delta_naive$d.Test)
dim(df_delta_naive)

r0 <-  glm(d.Test ~ V.Test + Age.FirstDose + Gender + BMI,family = "binomial", data = df_delta_naive)
summary(r0)
beta_v0 <- summary(r0)$coefficient[2,1]
sd_beta_v0 <- summary(r0)$coefficient[2,2]
ve0 <- 1 - exp(beta_v0)
sd_ve0<- exp(beta_v0)*sd_beta_v0
ve0
sd_ve0


df_delta_1p <- df_deltap %>%
  filter(TestingReason== 1) %>%
  mutate(d.Test = as.factor(ifelse(Result_PCR=="TRUE" & Result_Sym == "TRUE",1,0))) %>%
  mutate(V.Test = ifelse(Vax_Level>0, 1, 0))
dim(df_delta_1p)  
r1 <- glm(d.Test ~ V.Test + Age.FirstDose + Gender + BMI,family = "binomial", data = df_delta_1p) 
beta_v1 <- summary(r1)$coefficient[2,1]
sd_beta_v1 <- summary(r1)$coefficient[2,2]
ve1 <- 1 - exp(beta_v1)
sd_ve1<- exp(beta_v1)*sd_beta_v1

summary(r1)

df_delta_3p <- df_deltap %>%
  filter(TestingReason== 3) %>%
  mutate(d.Test = as.factor(ifelse(Result_PCR=="TRUE" & Result_Sym == "TRUE",1,0))) %>%
  mutate(V.Test = ifelse(Vax_Level>0, 1, 0))
dim(df_delta_3p)  
r3 <- glm(d.Test ~ V.Test + Age.FirstDose + Gender + BMI,family = "binomial", data = df_delta_3p) 
summary(r3)

beta_v3 <- summary(r3)$coefficient[2,1]
sd_beta_v3 <- summary(r3)$coefficient[2,2]
ve3 <- 1 - exp(beta_v3)
sd_ve3<- exp(beta_v3)*sd_beta_v3

## Bootstrap visualizations of mean VE vs. Age in the Delta Period

n_boot = 500
age_seq = seq(18,60,1)
formula_main = d.Test ~ V.Test * Age.FirstDose
xlim_vec = c(18,60)
ylim_vec = c(-1/2,1)
library(tidyverse)

fit_symptom <- glm(d.Test ~ V.Test * Age.FirstDose, data = df_delta_1p, family = binomial(link = "logit"))
fit_casetracing <- glm(d.Test ~ V.Test * Age.FirstDose, data = df_delta_3p, family = binomial(link = "logit"))
get_ve <- function(fitnow,age){
  
  df_vax <- data.frame( V.Test = 1,Age.FirstDose = age)
  df_unv <- data.frame( V.Test = 0,Age.FirstDose = age)
  rhs = ~V.Test * Age.FirstDose
  x_vax <- model.matrix(terms(rhs), df_vax)
  x_unv <- model.matrix(terms(rhs),df_unv)
  coefs <- coef(fitnow)
  delta <- sum(x_vax * coefs) - sum(x_unv * coefs)
  ve = 1-exp(delta)
  return (ve)
}
ve_symptom <- map_dbl(age_seq, ~get_ve(fit_symptom,.x))
ve_casetracing <- map_dbl(age_seq, ~get_ve(fit_casetracing,.x))
ve_mat_symptom <- matrix(NA, nrow = n_boot, ncol = length(age_seq))
ve_mat_casetracing <- ve_mat_symptom
set.seed(619)
N1 <- nrow(df_delta_1p)
N3 <- nrow(df_delta_3p)
for (b in 1:n_boot)
{
  idx1 <- sample(seq_len(N1), size = N1, replace = TRUE)
  idx3 <- sample(seq_len(N3), size = N3, replace = TRUE)
  data_b_symptom <- df_delta_1p[idx1,]
  data_b_casetracing <- df_delta_3p[idx3,]
  fit_b_symptom <- glm(formula_main, data = data_b_symptom, family = binomial(link = "logit"))
  fit_b_casetracing <- glm(formula_main, data = data_b_casetracing, family = binomial(link = "logit"))
  ve_b_symptom <- map_dbl(age_seq, ~get_ve(fit_b_symptom, .x))
  ve_mat_symptom[b,] <- ve_b_symptom
  ve_b_casetracing <- map_dbl(age_seq, ~get_ve(fit_b_casetracing, .x))
  ve_mat_casetracing[b,] <- ve_b_casetracing
  
  
}
ve_mean_symptom <- colMeans(ve_mat_symptom, na.rm = TRUE)
ve_lower_symptom <- apply(ve_mat_symptom, 2, quantile, probs = 0.025, na.rm = TRUE)
ve_upper_symptom <- apply(ve_mat_symptom, 2, quantile, probs = 0.975, na.rm = TRUE)

ve_mean_casetracing <- colMeans(ve_mat_casetracing, na.rm = TRUE)
ve_lower_casetracing <- apply(ve_mat_casetracing, 2, quantile, probs = 0.025, na.rm = TRUE)
ve_upper_casetracing <- apply(ve_mat_casetracing, 2, quantile, probs = 0.975, na.rm = TRUE)


plot_data <- tibble(
  AgeVal = age_seq,
  ve_mean_casetracing = ve_mean_casetracing,
  ve_lower_casetracing = ve_lower_casetracing,
  ve_upper_casetracing = ve_upper_casetracing,
  ve_mean_symptom = ve_mean_symptom,
  ve_lower_symptom = ve_lower_symptom,
  ve_upper_symptom = ve_upper_symptom
)

plot_data_symptoms <- plot_data %>%
  select(AgeVal, ve_mean_symptom, ve_lower_symptom, ve_upper_symptom) %>%
  rename(ve_lower = ve_lower_symptom, 
         ve_upper = ve_upper_symptom, 
         ve_mean = ve_mean_symptom, 
  ) %>%
  mutate(Group = "Symptoms")

plot_data_case <- plot_data %>%
  select(AgeVal, ve_mean_casetracing, ve_lower_casetracing,
         ve_upper_casetracing) %>%
  rename(ve_lower = ve_lower_casetracing, 
         ve_upper = ve_upper_casetracing, 
         ve_mean = ve_mean_casetracing, 
  ) %>%
  mutate(Group = "Case-Contact")

plot_data_long <- bind_rows(plot_data_case,plot_data_symptoms)

plot_data_long <- plot_data_long%>%mutate(Gp = factor(Group))

p <- ggplot(plot_data_long, aes(x = AgeVal, color = Group)) + 
  geom_ribbon(aes(ymin = ve_lower, ymax = ve_upper, fill = Group), 
              alpha = 0.4) +
  geom_line(aes(y=ve_mean, color = Group), size = 1) +
  coord_cartesian(xlim = xlim_vec, ylim = ylim_vec) +
  theme_bw(base_size = 14) + 
  labs(
    x="Age",
    y="VE",
    fill = "Group",
    color = "Group",
    title = "VE vs. Age (Delta Period)"
  )+    scale_fill_manual(values = c("Symptoms" = "red", "Case-Contact" = "blue"))+
  scale_color_manual(values = c("Symptoms" = "red", "Case-Contact" = "blue"))

## Data Analysis for Omicron Period

df_omicronp <- joined_data_final_agg %>%
  filter(Omicron_Period_Flag ==2)

df_omicron_naive <- df_omicronp%>%
  filter(TestingReason %in% c(1,3)) %>%
  mutate(d.Test = as.factor(ifelse(Result_PCR=="TRUE" & Result_Sym == "TRUE",1,0))) %>%
  mutate(V.Test = ifelse(Vax_Level>2, 1, 0))

r0 <-  glm(d.Test ~ V.Test + Age.FirstDose + Gender + BMI,family = "binomial", data = df_omicron_naive)
dim(df_omicron_naive)
summary(r0)
beta_v0 <- summary(r0)$coefficient[2,1]
sd_beta_v0 <- summary(r0)$coefficient[2,2]
ve0 <- 1 - exp(beta_v0)
sd_ve0<- exp(beta_v0)*sd_beta_v0
cat("ve_0 = ", ve0)
cat("\n sd_ve_0 =", sd_ve0)

table(df_omicron_naive$d.Test)

df_omicron_1p <- df_omicronp %>%
  filter(TestingReason== 1) %>%
  mutate(d.Test = as.factor(ifelse(Result_PCR=="TRUE" & Result_Sym == "TRUE",1,0))) %>%
  mutate(V.Test = ifelse(Vax_Level>2  , 1, 0) )
dim(df_omicron_1p)  
r1 <- glm(d.Test ~ V.Test + Age.FirstDose + Gender + BMI,family = "binomial", data = df_omicron_1p) 
beta_v1 <- summary(r1)$coefficient[2,1]
sd_beta_v1 <- summary(r1)$coefficient[2,2]
ve1 <- 1 - exp(beta_v1)
sd_ve1<- exp(beta_v1)*sd_beta_v1

summary(r1)
cat("ve_1 = ", ve1)
cat("\n sd_ve_1 =", sd_ve1)

df_omicron_3p <- df_omicronp %>%
  filter(TestingReason== 3) %>%
  mutate(d.Test = as.factor(ifelse(Result_PCR=="TRUE" & Result_Sym == "TRUE",1,0))) %>%
  mutate(V.Test = ifelse(Vax_Level>2 , 1, 0))
dim(df_omicron_3p)  
r3 <- glm(d.Test ~ V.Test + Age.FirstDose + Gender + BMI,family = "binomial", data = df_omicron_3p) 
summary(r3)


beta_v3 <- summary(r3)$coefficient[2,1]
sd_beta_v3 <- summary(r3)$coefficient[2,2]
ve3 <- 1 - exp(beta_v3)
sd_ve3<- exp(beta_v3)*sd_beta_v3
cat("ve_3 = ", ve3)
cat("\n sd_ve_3 =", sd_ve3)

## Bootstrap visualizations of mean VE vs. Age in the Omicron Period

n_boot = 500
age_seq = seq(18,60,1)
formula_main = d.Test ~ V.Test * Age.FirstDose
xlim_vec = c(18,60)
ylim_vec = c(-1/2,1)
library(tidyverse)

fit_symptom <- glm(d.Test ~ V.Test * Age.FirstDose, data = df_omicron_1p, family = binomial(link = "logit"))
fit_casetracing <- glm(d.Test ~ V.Test * Age.FirstDose, data = df_omicron_3p, family = binomial(link = "logit"))
get_ve <- function(fitnow,age){
  
  df_vax <- data.frame( V.Test = 1,Age.FirstDose = age)
  df_unv <- data.frame( V.Test = 0,Age.FirstDose = age)
  rhs = ~V.Test * Age.FirstDose
  x_vax <- model.matrix(terms(rhs), df_vax)
  x_unv <- model.matrix(terms(rhs),df_unv)
  coefs <- coef(fitnow)
  delta <- sum(x_vax * coefs) - sum(x_unv * coefs)
  ve = 1-exp(delta)
  return (ve)
}
ve_symptom <- map_dbl(age_seq, ~get_ve(fit_symptom,.x))
ve_casetracing <- map_dbl(age_seq, ~get_ve(fit_casetracing,.x))
ve_mat_symptom <- matrix(NA, nrow = n_boot, ncol = length(age_seq))
ve_mat_casetracing <- ve_mat_symptom
set.seed(619)
N1 <- nrow(df_omicron_1p)
N3 <- nrow(df_omicron_3p)
for (b in 1:n_boot)
{
  idx1 <- sample(seq_len(N1), size = N1, replace = TRUE)
  idx3 <- sample(seq_len(N3), size = N3, replace = TRUE)
  data_b_symptom <- df_omicron_1p[idx1,]
  data_b_casetracing <- df_omicron_3p[idx3,]
  fit_b_symptom <- glm(formula_main, data = data_b_symptom, family = binomial(link = "logit"))
  fit_b_casetracing <- glm(formula_main, data = data_b_casetracing, family = binomial(link = "logit"))
  ve_b_symptom <- map_dbl(age_seq, ~get_ve(fit_b_symptom, .x))
  ve_mat_symptom[b,] <- ve_b_symptom
  ve_b_casetracing <- map_dbl(age_seq, ~get_ve(fit_b_casetracing, .x))
  ve_mat_casetracing[b,] <- ve_b_casetracing
  
  
}
ve_mean_symptom <- colMeans(ve_mat_symptom, na.rm = TRUE)
ve_lower_symptom <- apply(ve_mat_symptom, 2, quantile, probs = 0.025, na.rm = TRUE)
ve_upper_symptom <- apply(ve_mat_symptom, 2, quantile, probs = 0.975, na.rm = TRUE)

ve_mean_casetracing <- colMeans(ve_mat_casetracing, na.rm = TRUE)
ve_lower_casetracing <- apply(ve_mat_casetracing, 2, quantile, probs = 0.025, na.rm = TRUE)
ve_upper_casetracing <- apply(ve_mat_casetracing, 2, quantile, probs = 0.975, na.rm = TRUE)

plot_data <- tibble(
  AgeVal = age_seq,
  ve_mean_casetracing = ve_mean_casetracing,
  ve_lower_casetracing = ve_lower_casetracing,
  ve_upper_casetracing = ve_upper_casetracing,
  ve_mean_symptom = ve_mean_symptom,
  ve_lower_symptom = ve_lower_symptom,
  ve_upper_symptom = ve_upper_symptom
)

plot_data_symptoms <- plot_data %>%
  select(AgeVal, ve_mean_symptom, ve_lower_symptom, ve_upper_symptom) %>%
  rename(ve_lower = ve_lower_symptom, 
         ve_upper = ve_upper_symptom, 
         ve_mean = ve_mean_symptom, 
  ) %>%
  mutate(Group = "Symptoms")

plot_data_case <- plot_data %>%
  select(AgeVal, ve_mean_casetracing, ve_lower_casetracing,
         ve_upper_casetracing) %>%
  rename(ve_lower = ve_lower_casetracing, 
         ve_upper = ve_upper_casetracing, 
         ve_mean = ve_mean_casetracing, 
  ) %>%
  mutate(Group = "Case-Contact")

plot_data_long <- bind_rows(plot_data_case,plot_data_symptoms)
plot_data_long <- plot_data_long%>%mutate(Gp = factor(Group))

p2 <- ggplot(plot_data_long, aes(x = AgeVal, color = Group)) + 
  
  geom_ribbon(aes(ymin = ve_lower, ymax = ve_upper, fill = Group), 
              alpha = 0.4) +
  geom_line(aes(y=ve_mean, color = Group), size = 1) +
  coord_cartesian(xlim = xlim_vec, ylim = ylim_vec) +
  theme_bw(base_size = 14) + 
  labs(
    x="Age",
    y="VE",
    fill = "Group",
    color = "Group",
    title = "VE vs. Age (Omicron Period)"
  )+    scale_fill_manual(values = c("Symptoms" = "red", "Case-Contact" = "blue"))+
  scale_color_manual(values = c("Symptoms" = "red", "Case-Contact" = "blue"))

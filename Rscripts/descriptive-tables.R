# this script produces summary tables for baseline 
# and for dates of entry/withdrawal from the study
# the datasets used here are produced 
# in the script "raw-to-td" 



#############################################
# load libraries 
library(dplyr)
library(tidyr)



# time from randomization to treatment discountinuation
# time from randomization to withdrawal/lost to follow up

# load data
setwd('C:/Users/ericapo/Desktop/CTU/Dividint/DividintStat/data/')
dates.db <- read.csv("td/dates.csv", sep = '')
#head(dates.db)

# from character to dates 
db <- dates.db %>% mutate(rand = as.Date(dates.db$date_rand, format = "%d-%m-%Y"))
db <- db %>% mutate(out = as.Date(dates.db$date_end, format = "%d-%m-%Y"))
db <- db %>% mutate(timetoend = out-rand)


# load treatment information
### FOR UNBLINDED ANALYSIS: 
### CHANGE NAME OF INPUT FILE:
### UNCOMMENT FOLLOWING LINES
#treatment_alloc = read.table(file="td/TreatmentAllocFake.csv",sep="")
#names(treatment_alloc) = c("study_id","treatment")
#treatment_alloc$trt[treatment_alloc$treatment == "Placebo"] = 0
#treatment_alloc$trt[treatment_alloc$treatment == "Pleconaril/Ribavirin"] = 1
#treatment_alloc$trt = as.numeric(treatment_alloc$trt)

# exclude those who withdrew before randomization
db <- inner_join(db, treatment_alloc)

table_w <- db %>% group_by(treatment) %>% summarise(withdrew = sum(!is.na(date_withdraw)))
table_in <- db %>% group_by(treatment) %>% summarise(included = sum(!is.na(date_baseline)))
table_lf <- db %>% group_by(treatment) %>% filter(!is.na(date_firstyear)) %>% summarise(lf = sum(!is.na(date_withdraw)))
table_c <- db %>% group_by(treatment) %>% summarise(completed = sum(!is.na(date_last)))
table_c1 <- db %>% group_by(treatment) %>% summarise(completedfirstyear = sum(!is.na(date_firstyear)))

table1 <- left_join(table_in, table_c1)
table1 <- left_join(table1, table_c)
table1 <- left_join(table1, table_w)
table1 <- left_join(table1, table_lf)


# report median and IQR of time to end 
# summarized by group
Median <- function(x) {median(x, na.rm = TRUE)}
IQR <- function(x) {quantile(x, probs = 0.75, na.rm = TRUE) - quantile(x, probs = 0.25, na.rm = TRUE)}
IQR75 <- function(x) {quantile(x, probs = 0.75, na.rm = TRUE)}
IQR25 <- function(x) {quantile(x, probs = 0.25, na.rm = TRUE)}

IQf <- function(x) {
  IQR <- IQR(x)
  IQR75 <- IQR75(x)
  IQR25 <- IQR25(x)
  string <- paste(IQR, ' days (', IQR25, '-', IQR75, ')', sep ='')
  return(string)
}

stats = c('Median', 'IQf')
table_dates <- db %>% group_by(treatment) %>% summarise_at("timetoend", stats) 
#setwd('..')
#save(table_dates, file = 'results/table_dates.RData')

# patients demographics and baseline characteristics
# age (years), gender, BMI, 
# medical conditions
# family member medical condition
# menarchal status
# breast development girls
# pubic hair girls
# pubic hair boys 
# tanner stage

# first (baseline) characteristics are saved in ADSL

baseline_db <- read.csv('td/ADSL.csv', sep = '')
# report n, mean and SD of age and BMI
# report n and % for gender
# summarized by group and overall
No = function(x) {sum(!is.na(x))}
Mean = function(x) {mean(x, na.rm = TRUE)}
Sd = function(x) {sd(x, na.rm = TRUE)}
stats = c('Mean', 'Sd','No')
stats_cat = c('n', 'freq')

table_baseline <- baseline_db %>% 
  summarise_at(c('age', 'bmi', 'weight' ), stats) 

table_baseline_bt <- baseline_db %>% group_by(treatment) %>% 
                  summarise_at(c('age', 'bmi', 'weight' ), stats) 

table_gender <- baseline_db %>% group_by(gender) %>%  summarise(n_gender = n()) %>%
  mutate(freq_gender = n_gender / sum(n_gender)*100) 

table_gender_bt <- baseline_db %>% group_by(treatment, gender) %>%  summarise(n_gender = n()) %>%
  mutate(freq_gender = n_gender / sum(n_gender)*100) 

# other info is in raw data

# load raw data in list

file_list <- list.files(path="locked/")

data <- lapply(file_list, function(x) read.csv(file=paste("locked/", x, sep = ''),
                                               sep = ";", header = TRUE))
names(data) <- file_list

# medical condition
# check all file names!
medicalcond <- data$DiViDint_MEDICAL_HISTORY_2022.08.31.csv
# filter on visit 2 only (baseline)
medicalcond <- medicalcond %>% filter(visit == 'DIVIDInt Visit 2: Baseline')
medicalcond <- medicalcond  %>% filter(study_id %in% id_list)

# filter on diagnosed = yes
medicalcond <- medicalcond %>% filter(diagnosed == 'YES')
# exclude those who withdrew before randomization
medicalcond <- inner_join(medicalcond, treatment_alloc)


Ntot <- nrow(ADSL)
table_medcond <- medicalcond %>% group_by(condition)  %>% summarise(n_cond = n()) %>%
  mutate(freq_cond = n_cond / Ntot*100) 



table_medcond_bt <- medicalcond %>% group_by(treatment, condition) %>% summarise(n_cond = n()) 

Nactive <- nrow(ADSL[ADSL$treatment == 'Pleconaril/Ribavirin', ])
Nplacebo <- nrow(ADSL[ADSL$treatment == 'Placebo', ])

table_medcond_bt$freq_cond <- 0
table_medcond_bt[table_medcond_bt$treatment == 'Pleconaril/Ribavirin', ]$freq_cond <- table_medcond_bt[table_medcond_bt$treatment == 'Pleconaril/Ribavirin', ]$n_cond/Nactive
table_medcond_bt[table_medcond_bt$treatment == 'Placebo', ]$freq_cond <- table_medcond_bt[table_medcond_bt$treatment == 'Placebo', ]$n_cond/Nplacebo

# family member medical condition
fammedicalcond <- data$DiViDint_FAMILY_MEDICAL_HISTORY_2022.08.31.csv

# filter on diagnosed = yes
fammedicalcond <- fammedicalcond %>% filter(fam_diagnosed == 'YES')
# exclude those who withdrew before randomization
fammedicalcond <- inner_join(fammedicalcond, treatment_alloc)

fammedicalcond <- fammedicalcond  %>% filter(study_id %in% id_list)

table_fammedcond <- fammedicalcond %>% group_by(fam_condition)  %>% summarise(n_cond = n()) %>%
  mutate(freq_cond = n_cond / Ntot*100) 

table_fammedcond_bt <- fammedicalcond %>% group_by(treatment, fam_condition) %>% summarise(n_cond = n()) 

table_fammedcond_bt$freq_cond <- 0
table_fammedcond_bt[table_fammedcond_bt$treatment == 'Pleconaril/Ribavirin', ]$freq_cond <- table_fammedcond_bt[table_fammedcond_bt$treatment == 'Pleconaril/Ribavirin', ]$n_cond/Nactive
table_fammedcond_bt[table_fammedcond_bt$treatment == 'Placebo', ]$freq_cond <- table_fammedcond_bt[table_fammedcond_bt$treatment == 'Placebo', ]$n_cond/Nplacebo



# menarche and pubic hair in clinical status
clinstat <- data$DiViDint_CLINICAL_STATUS_2022.08.31.csv

clinstat <- clinstat  %>% filter(study_id %in% id_list)
# code all in the same way!
clinstat[clinstat$pubic_hair_boys == '1', ]$pubic_hair_boys <-'P1'
clinstat[clinstat$pubic_hair_boys == '4', ]$pubic_hair_boys <-'P4'
clinstat[clinstat$pubic_hair_boys == 'p1', ]$pubic_hair_boys <-'P1'
clinstat[clinstat$pubic_hair_boys == 'p2', ]$pubic_hair_boys <-'P2'



# filter on visit 2 only (baseline)
clinstat <- clinstat %>% filter(visit == 'DIVIDInt Visit 2: Baseline')
# exclude those who withdrew before randomization
clinstat <- inner_join(clinstat, treatment_alloc)

# add information about gender
clinstat <- inner_join(clinstat, baseline_db)

table_menarche <- clinstat  %>% group_by(menarche) %>% summarise(n = n()) %>%
  mutate(freq = n / Ntot*100) 

table_breast <- clinstat  %>% filter(gender == 'FEMALE') %>% 
  group_by(breast_development_girls) %>% summarise(n = n()) %>%
  mutate(freq = n / Ntot*100) 

table_phg <- clinstat %>% filter(gender == 'FEMALE')%>% 
             group_by(pubic_hair_girls) %>% summarise(n = n()) %>%
  mutate(freq = n / Ntot*100) 

table_phb <- clinstat  %>% filter(gender == 'MALE')%>% 
  group_by(pubic_hair_boys) %>% summarise(n = n()) %>%
  mutate(freq = n / Ntot*100) 

table_menarche_bt <- clinstat  %>% group_by(treatment,menarche) %>% summarise(n = n()) 

table_menarche_bt$freq <- 0
table_menarche_bt[table_menarche_bt$treatment == 'Pleconaril/Ribavirin', ]$freq <- table_menarche_bt[table_menarche_bt$treatment == 'Pleconaril/Ribavirin', ]$n/Nactive
table_menarche_bt[table_menarche_bt$treatment == 'Placebo', ]$freq <- table_menarche_bt[table_menarche_bt$treatment == 'Placebo', ]$n/Nplacebo


table_breast_bt <- clinstat  %>% filter(gender == 'FEMALE') %>% 
  group_by(treatment,breast_development_girls) %>% summarise(n = n()) 



table_breast_bt$freq <- 0
table_breast_bt[table_breast_bt$treatment == 'Pleconaril/Ribavirin', ]$freq <- table_breast_bt[table_breast_bt$treatment == 'Pleconaril/Ribavirin', ]$n/Nactive
table_breast_bt[table_breast_bt$treatment == 'Placebo', ]$freq <- table_breast_bt[table_breast_bt$treatment == 'Placebo', ]$n/Nplacebo


table_phg_bt <- clinstat %>% filter(gender == 'FEMALE')%>% 
  group_by(treatment,pubic_hair_girls) %>% summarise(n = n()) 

table_phg_bt$freq <- 0
table_phg_bt[table_phg_bt$treatment == 'Pleconaril/Ribavirin', ]$freq <- table_phg_bt[table_phg_bt$treatment == 'Pleconaril/Ribavirin', ]$n/Nactive
table_phg_bt[table_phg_bt$treatment == 'Placebo', ]$freq <- table_phg_bt[table_phg_bt$treatment == 'Placebo', ]$n/Nplacebo


table_phb_bt <- clinstat  %>% filter(gender == 'MALE')%>% 
  group_by(treatment, pubic_hair_boys) %>% summarise(n = n()) 


table_phb_bt$freq <- 0
table_phb_bt[table_phb_bt$treatment == 'Pleconaril/Ribavirin', ]$freq <- table_phb_bt[table_phb_bt$treatment == 'Pleconaril/Ribavirin', ]$n/Nactive
table_phb_bt[table_phb_bt$treatment == 'Placebo', ]$freq <- table_phb_bt[table_phb_bt$treatment == 'Placebo', ]$n/Nplacebo


# add height
height <- data$DiViDint_BMI_2022.08.31.csv
# filter on visit 2 only (baseline)
height <- height %>% filter(visit == 'DIVIDInt Visit 2: Baseline')
# exclude those who withdrew before randomization
height <- inner_join(height, treatment_alloc)
height <- height  %>% filter(study_id %in% id_list)



table_height  <- height  %>% 
  summarise_at(c('height' ), stats) 

names(table_height) <- c('height_Mean', 'height_Sd', 'height_No')
table_baseline <- cbind(table_baseline, table_height)

table_height_bt <- height  %>% group_by(treatment) %>% 
  summarise_at(c('height'), stats) 

names(table_height_bt) <- c('treatment', 'height_Mean', 'height_Sd', 'height_No')
table_baseline_bt <- inner_join(table_baseline_bt, table_height_bt)

# merge all tables to obtain table 7 and 8 from masterSAP
# save two tables one overall and one by treatment

# continuous variables to be reformatted
# as mean (SD)
cont_var <- c('bmi', 'height', 'weight', 'age')
table_baseline_t <- table_baseline_bt %>% filter(treatment != 'Placebo')
table_baseline_c <- table_baseline_bt %>% filter(treatment == 'Placebo')
summary_cont <- c()
summary_cont_treat <- c()
summary_cont_ctrl <- c()

N_cont <- c()
N_cont_treat <- c()
N_cont_ctrl <- c()


for (var in cont_var){
  # look for correct column in the table
  funpaste <- function(x) {paste(var, x, sep = '_')}
  names <- sapply(stats, funpaste)
  names <- as.character(names)
  # rewrite with correct format for summary
  summary_cont[var] <-  paste(round(table_baseline %>% pull(names[1]), 2), 
                                    '(', round(table_baseline %>% pull(names[2]), 2),
                                    ')', sep ='')
  
  N_cont[var] <- table_baseline %>% pull(names[3])
  
  
  summary_cont_treat[var] <-  paste(round(table_baseline_t %>% pull(names[1]), 2), 
                              '(', round(table_baseline_t %>% pull(names[2]), 2),
                              ')', sep ='')
  
  summary_cont_ctrl[var] <-  paste(round(table_baseline_c %>% pull(names[1]), 2), 
                              '(', round(table_baseline_c %>% pull(names[2]), 2),
                              ')', sep ='')
  
  N_cont_treat[var] <- table_baseline_t %>% pull(names[3])
  
  N_cont_ctrl[var] <- table_baseline_c %>% pull(names[3])
 
}


# current insulin regimen by group and visit
insulin <- data$DiViDint_DIABETES_CARE_2022.08.31.csv
insulin <- insulin  %>% filter(study_id %in% id_list)


# exclude those who withdrew before randomization
insulin<- inner_join(insulin, treatment_alloc)
insulin_baseline <- insulin %>% filter(visit == 'DIVIDInt Visit 2: Baseline')

table_insulin <- insulin_baseline %>% group_by(current_insulin_regimen)  %>% summarise(n = n()) %>%
  mutate(freq = n / Ntot*100) 

table_insulin_bt <- insulin_baseline %>% group_by(treatment, current_insulin_regimen)  %>% summarise(n = n()) 


table_insulin_bt$freq <- 0
table_insulin_bt[table_insulin_bt$treatment == 'Pleconaril/Ribavirin', ]$freq <- table_insulin_bt[table_insulin_bt$treatment == 'Pleconaril/Ribavirin', ]$n/Nactive
table_insulin_bt[table_insulin_bt$treatment == 'Placebo', ]$freq <- table_insulin_bt[table_insulin_bt$treatment == 'Placebo', ]$n/Nplacebo

# insulin at each visit

table_insulin_byvisit <- insulin %>% filter(visit != 'DIVIDInt Visit 1: Inclusion') %>%  group_by(visit, current_insulin_regimen) %>% summarise(n = n()) %>%
  mutate(freq = n / Ntot*100) 

table_insulin_byvisit_bt <- insulin %>% filter(visit != 'DIVIDInt Visit 1: Inclusion') %>%  group_by(visit, treatment, current_insulin_regimen) %>% summarise(n = n()) 


table_insulin_byvisit_bt$freq <- 0
table_insulin_byvisit_bt[table_insulin_byvisit_bt$treatment == 'Pleconaril/Ribavirin', ]$freq <- table_insulin_byvisit_bt[table_insulin_byvisit_bt$treatment == 'Pleconaril/Ribavirin', ]$n/Nactive
table_insulin_byvisit_bt[table_insulin_byvisit_bt$treatment == 'Placebo', ]$freq <- table_insulin_byvisit_bt[table_insulin_byvisit_bt$treatment == 'Placebo', ]$n/Nplacebo



summary_baseline_overall <- data.frame(cont_var, N_cont, summary_cont)

summary_baseline_bytreatment <- data.frame(cont_var,  
                                           N_cont_treat,
                                           summary_cont_treat,
                                           N_cont_ctrl,
                                           summary_cont_ctrl)


summary_overall <- list(summary_baseline_overall, 
                        table_breast, 
                        table_fammedcond, 
                        table_medcond, 
                        table_gender, 
                        table_menarche, 
                       # table_phb, 
                       # table_phg, 
                        table_insulin)

summary_bytreatment <- list(summary_baseline_bytreatment, 
                        table_breast_bt, 
                        table_fammedcond_bt, 
                        table_medcond_bt, 
                        table_gender_bt, 
                        table_menarche_bt, 
                       # table_phb_bt, 
                       # table_phg_bt, 
                        table_insulin_bt)

names(summary_overall) <- names(summary_bytreatment) <- c('baseline', 
                                                          'breast', 
                                                          'fammedcond', 
                                                          'medcond', 
                                                          'gender', 
                                                          'menarche', 
                                                        #  'phb', 
                                                        #  'phg', 
                                                          'insulin')

#setwd('..')
#save(summary_overall, file = 'results/table_summaryoverall.RData')
#save(summary_bytreatment, file = 'results/table_summarybytreatment.RData')

visits <- levels(table_insulin_byvisit_bt$visit)[-1]
tab_insreg <- list()
for (vv in 1:length(visits)){
  tab_tr <- table_insulin_byvisit_bt %>% filter(treatment == 'Pleconaril/Ribavirin') %>%
    filter(visit == visits[vv])
  
  tab_pl <- table_insulin_byvisit_bt %>% filter(treatment == 'Placebo') %>%
    filter(visit == visits[vv])
  
  db <- merge(tab_tr, tab_pl, by = "current_insulin_regimen", all =TRUE) 
  db <- db %>%  mutate(N.tr = replace_na(n.x, 0), 
                       Freq.tr = replace_na(freq.x, 0), 
                       N.pl = replace_na(n.y, 0), 
                       Freq.pl = replace_na(freq.y, 0))  %>%
    dplyr::select(current_insulin_regimen, N.tr, Freq.tr, N.pl, Freq.pl)
  colnames(db) <- c("current_insulin_regimen", 
                    paste("N.tr", vv, sep =''), 
                    paste("Freq.tr", vv, sep =''), 
                    paste("N.pl", vv, sep =''), 
                    paste("Freq.pl", vv, sep =''))
  
  if (vv == 1) tab_insreg  <- db
  if (vv > 1) tab_insreg  <- merge(tab_insreg, db, by = "current_insulin_regimen", all =TRUE) 
  
}

tab_insreg[is.na(tab_insreg)] = 0

save(tab_insreg , file ='../results/tab_insreg.RData')








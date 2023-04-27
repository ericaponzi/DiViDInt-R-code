# this script produces the datasets in "td" folder
# from the raw data to ADSL data and dates dataset

# structure of ADSL data: 
# patient ID, treatment, stratification factors, FAS, PP, data rand, data left
# stratification factors here are: 
# site, age at baseline, gender, bmi at baseline
# FAS and PP: 
# exclude the 10 who withdrew
# add compliance information



#############################################
# load libraries 
library(dplyr)
library(tidyr)

# extract data

### FOR UNBLINDED ANALYSIS: 
### change this two lines to local directory
setwd("../data/")
file_list <- list.files(path="locked/")

data <- lapply(file_list, function(x) read.csv(file=paste("locked/", x, sep = ''),
                                               sep = ";", header = TRUE))
names(data) <- file_list
#file_list



# allocation dataset with only two columns: id and treatment
# !!!!!
# create fake allocation variable for blinded analysis
### FOR UNBLINDED ANALYSIS: 
### COMMENT THE FOUR FOLLOWING LINES 
id_list <- unique(data[[4]]$study_id) #these should be 98 # check that the file is the BMI one
id_list <- as.character(id_list)

# the following IDs have to be excluded (total 96)
which(id_list == 'DVD-06001-00001')
which(id_list == 'DVD-19000-00038')
id_list <- id_list[-c(83, 92)]

# set.seed(4224)
# alloc <- sample(x=c(0,1), size=length(id_list), replace = T, prob = NULL)
# treatment_alloc <- data.frame(study_id = id_list, trt = alloc)
# create character variable for labels in plots
# treatment_alloc$treatment[treatment_alloc$trt == 0] = "Placebo"
# treatment_alloc$treatment[treatment_alloc$trt == 1] = "Active"
# head(treatment_alloc)
# table(treatment_alloc$trt)
# write.table(treatment_alloc, 'td/TreatmentAllocFake.csv')





# Demographics contains information on age and sex, and site
Demo <- inner_join(data$DiViDint_DEMOGRAPHICS_2022.08.31.csv, 
                   data$DiViDint_VISIT_METADATA_2022.08.31.csv)

# this contains two patients extra than the 96
# they withdrew
# exclude them
Demo <- Demo %>% filter(as.character(study_id) %in% id_list)

# extract date of birth and date of visit
Demo <- Demo %>% separate(date_of_birth,
                          into = c('birthmonth', 'birthyear'), sep = '-')
# some checks 
#table(Demo$birthmonth)
#table(Demo$birthyear)

Demo <- Demo %>% separate(date_of_visit,
                          into = c('visitday', 'visitmonth', 'visityear'), sep = '-')
Demo[Demo$visitday == '', ]$visitday <- NA 

# some checks 
#table(Demo$visitday)
#table(Demo$visitmonth)
#table(Demo$visityear)

#sum(is.na(Demo$visitday))
#sum(is.na(Demo$visitmonth))
#sum(is.na(Demo$visityear))

#which(is.na(Demo$visitday)) == which(is.na(Demo$visitmonth))
#which(is.na(Demo$visitday)) == which(is.na(Demo$visityear))


# calculate age as difference between dates
Demo$age <- as.numeric(Demo$visityear)-
  as.numeric(Demo$birthyear) +
  (as.numeric(Demo$visitmonth)-as.numeric(Demo$birthmonth))/12


#summary(Demo$age)


# only age we will use as stratification factor in models will be the 
# age at baseline (pre treatment)
# extract age at baseline, gender and site information
AgeBaseline <- Demo %>% filter(visit == 'DIVIDInt Visit 2: Baseline') %>%
  dplyr::select(c('study_id', 'site', 'gender', 'age'))




# BMI information: 
BMI <- data$DiViDint_BMI_2022.08.31.csv


# this contains two patients extra than the 96
# they withdrew
# exclude them
BMI <- BMI %>% filter(as.character(study_id) %in% id_list)


# only bmi we will use as stratification factor in models will be the 
# bmi at baseline (pre treatment)
BMIBaseline <- BMI %>% filter(visit == 'DIVIDInt Visit 2: Baseline') %>%
  dplyr::select(c('study_id', 'bmi', 'weight'))



# information about dates of randomization (inclusion visit) and dates they left 

dates <- data$DiViDint_VISIT_METADATA_2022.08.31.csv
dates <- dates %>% filter(as.character(study_id) %in% id_list)

date_rand <- dates %>% filter(visit =='DIVIDInt Visit 1: Inclusion') %>% 
  dplyr::select(study_id, date_of_visit)
names(date_rand) <- c('study_id', 'date_rand')

date_firstyear <- dates %>% filter(visit =='DIVIDInt Visit 7: 1 year from first IMP') %>% 
  dplyr::select(study_id, date_of_visit)
names(date_firstyear) <- c('study_id', 'date_firstyear')

date_firstyear_extra <- dates %>% filter(visit =='DIVIDInt Visit 7: 1 year from first IMP EXTRA VISIT') %>% 
  dplyr::select(study_id, date_of_visit)
names(date_firstyear_extra) <- c('study_id', 'date_firstyear')

date_last <- dates %>% filter(visit =='DIVIDInt Visit 9: 3 years from first IMP') %>% 
  dplyr::select(study_id, date_of_visit) 

names(date_last) <- c('study_id', 'date_last')


# date of withdraw
date_withdrawal <- dates %>% filter(visit =='DVD Withdrawal Visit') %>% 
  dplyr::select(study_id, date_of_visit) 
names(date_withdrawal) <- c('study_id', 'date_withdraw')





date <- left_join(date_rand, date_last)
date <- left_join(date, date_firstyear)
date <- left_join(date, date_withdrawal)
# end date is either last visit or withdrawal visit
date$date_left <- date$date_last
date[!is.na(date$date_withdraw),
]$date_left <- date[!is.na(date$date_withdraw), 
]$date_withdraw

date <- date %>%  dplyr::select(study_id, date_rand, date_firstyear, date_left)


# pool data together

ADSL <- inner_join(BMIBaseline, AgeBaseline)
ADSL <- inner_join(ADSL, treatment_alloc)
# this excludes the two extra from above 
# (they withdrew before treatment was assigned)
ADSL <- inner_join(ADSL, date)

# FAS: all individuals
ADSL$FAS <- 1 

withdraw <- data$DiViDint_WITHDRAWAL_2022.08.31.csv
withdraw <- withdraw[-which(withdraw$study_id == "DVD-06001-00003"),]

withdraw <- withdraw %>% filter(study_id %in% id_list)
# compliance
# extract compliance data from IMP file (extended version!)



IMP <- data$DiViDint_IMP_EXTENDED_2022.08.31.csv

IMP <- IMP %>% filter(as.character(study_id) %in% id_list)


#head(IMP)
#str(IMP)


# according to meeting w Ida on 29.11
# compliance is obtained as ratio between:
# "IMP$v6_first_disp_used_amount_pleconaril_placebo" 
# and 
# "v6_first_disp_estimated_dose_pleconaril_placebo"
# these values should be only available at 6 months


# filter data to 6months visit only
Comp <- IMP %>% filter(visit == 'DIVIDInt Visit 6: 6 months from first IMP') %>% 
  dplyr::select('study_id', 'v6_first_disp_used_amount_pleconaril_placebo', 'v6_first_disp_used_amount_ribavirin_placebo', 
         'v6_first_disp_estimated_dose_pleconaril_placebo', 'v6_first_disp_estimated_dose_ribavirin_placebo', 
         'v6_total_used_amount_pleconaril_placebo', 'v6_total_estimated_dose_pleconaril_placebo', 
         'v6_total_used_amount_ribavirin_placebo', 'v6_total_estimated_dose_ribavirin_placebo')


Comp <- Comp %>% mutate(Compliance_Plec = v6_total_used_amount_pleconaril_placebo/v6_total_estimated_dose_pleconaril_placebo) %>%
  mutate(Compliance_Rib = v6_total_used_amount_ribavirin_placebo/v6_total_estimated_dose_ribavirin_placebo)%>% 
  mutate(Compliance_Plec3 = v6_first_disp_used_amount_pleconaril_placebo/v6_first_disp_estimated_dose_pleconaril_placebo) %>%
  mutate(Compliance_Rib3 = v6_first_disp_used_amount_ribavirin_placebo/v6_first_disp_estimated_dose_ribavirin_placebo)


# add this to ADSL data
Compliance <- Comp %>%  dplyr::select(study_id, Compliance_Plec, Compliance_Rib, Compliance_Plec3, Compliance_Rib3)
ADSL <- left_join(ADSL, Compliance)












# structure of dates dataset: 
# patient id, inclusion date, randomization date, 1treatment, last treatment, 
# end of study, withdrawal yes/no, reason for withdrawing
# this will be used for flowchart

# date of inclusion
date_incl <- dates %>% filter(visit == 'DIVIDInt Visit 1: Inclusion')%>% 
  dplyr::select(study_id, date_of_visit)
names(date_incl) <- c('study_id', 'date_inclusion')

# date of baseline visit (also 1st treatment)
date_first <- dates %>% filter(visit =='DIVIDInt Visit 2: Baseline') %>% 
  dplyr::select(study_id, date_of_visit)
names(date_first) <- c('study_id', 'date_baseline')

# date of one year visit 
date_firstyear <- dates %>% filter(visit =='DIVIDInt Visit 7: 1 year from first IMP') %>% 
  dplyr::select(study_id, date_of_visit)
names(date_firstyear) <- c('study_id', 'date_firstyear')

# date of last visit
date_last <- dates %>% filter(visit =='DIVIDInt Visit 9: 3 years from first IMP') %>% 
  dplyr::select(study_id, date_of_visit) 
names(date_last) <- c('study_id', 'date_last')

# date of withdrawal
date_withdrawal <- dates %>% filter(visit =='DVD Withdrawal Visit') %>% 
  dplyr::select(study_id, date_of_visit) 
names(date_withdrawal) <- c('study_id', 'date_withdraw')

# pool data together
dates_all <- left_join(date_incl, date_rand)
dates_all <- left_join(dates_all, date_first)
dates_all <- left_join(dates_all, date_firstyear)
dates_all <- left_join(dates_all, date_last)
dates_all <- left_join(dates_all, date_withdrawal)

dates_all <- dates_all[dates_all$study_id %in% id_list, ]

# end date is either last visit or withdrawal visit
dates_all$date_end <- dates_all$date_last
dates_all[!is.na(dates_all$date_withdraw),
          ]$date_end <- dates_all[!is.na(dates_all$date_withdraw), 
                                  ]$date_withdraw

# add reason for withdraw
reason <- withdraw %>%  dplyr::select(study_id, withdrawal_reason)

dates_all <- left_join(dates_all, reason)


# save dates data
write.table(dates_all, 'td/dates.csv')


# add compliance information to the ADSL data

withdrawn_ever <- withdraw$study_id
withdrawn_first_db <- dates_all %>% filter(!is.na(date_firstyear)) %>% filter(!is.na(date_withdraw)) 
withdrawn_first <- withdrawn_first_db$study_id #withdrew after first year visit (so they have primary)


ADSL$withdrawn_ever <- 0
ADSL[ADSL$study_id %in% withdrawn_ever, ]$withdrawn_ever <- 1

ADSL$withdrawn_first <- 0
ADSL[ADSL$study_id %in% withdrawn_first, ]$withdrawn_first <- 1




ADSL$PP70 <- ADSL$PP60 <- ADSL$PP80 <- 1
ADSL$avg.compl <- (ADSL$Compliance_Plec + ADSL$Compliance_Rib)/2
noncomp70 <- which(ADSL$avg.compl < 0.7)
noncomp60 <- which(ADSL$avg.compl < 0.6)
noncomp80 <- which(ADSL$avg.compl < 0.8)


  
ADSL[noncomp70, ]$PP70 <- 0
ADSL[noncomp60, ]$PP60 <- 0
ADSL[noncomp80, ]$PP80 <- 0

ADSL[is.na(ADSL$Compliance_Plec), ]$PP70 <- 0
ADSL[is.na(ADSL$Compliance_Plec), ]$PP60 <- 0
ADSL[is.na(ADSL$Compliance_Plec), ]$PP80 <- 0

# save ADSL data
write.table(ADSL, 'td/ADSL.csv')

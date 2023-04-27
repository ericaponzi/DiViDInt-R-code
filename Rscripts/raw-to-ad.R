# this script produces the datasets in "ad" folder
# from the raw data and ADSL data to datasets for 
# statistical analysis



# efficacy analysis
# this dataset contains all information for the analysis of
# all primary and secondary outcomes
# structure of the dataset: 
# study_id, stratification variables (from ADLS data), treatment, 
# FAS, PP, all outcomes

# the following is used to extract all outcomes and then merge with ADSL dataset
# this script produces data in longitudinal format


# load libraries 
library(dplyr)
library(tidyr)
library(readxl)

setwd("../data/")
file_list <- list.files(path="locked/")

data <- lapply(file_list, function(x) read.csv(file=paste("locked/", x, sep = ''),
                                               sep = ";", header = TRUE))
names(data) <- file_list
#file_list


# PRIMARY OUTCOME: change from baseline in AUC after 12 months
# Secondary outcomes: change from baseline in AUC after 3, 6, 24 and 36 months
# calculate AUC for each visit

#MMTTold <- data$erica.ponzi_2022_6_29_12_9_10_DVD_COMBINED_C_PEPTIDE.csv
#MMTT <- read_excel('raw/MMTT.xlsx')
MMTT <- read.csv('locked/DiViDint_MMTT_C-peptide 2022.06.30.csv', sep = ';')
MMTT$result_c_pep_mmtt <- as.numeric(as.character(MMTT$result_c_pep_mmtt)) #!!!

MMTT <- MMTT %>% filter(study_id %in% id_list)

# formula for AUC (see mSAP)
# the -10 value is not used 
# (MMTTpep0+MMTTpep15) /3*15+(MMTTpep15+MMTTpep30)/2*15+(MMTTpep60+MMTTpep30)/2*30+(MMTTpep60+MMTTpep90)/2*30+(MMTTpep120+MMTTpep90)/2*30)/120000
# then transformed on the log scale

# calculate AUC for baseline visit

Baseline0 <- MMTT %>% filter(visit == 'DIVIDInt Visit 2: Baseline') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide 0 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Baseline0) <- c('study_id', 'AUC0')

Baseline15 <- MMTT %>% filter(visit == 'DIVIDInt Visit 2: Baseline') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +15 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Baseline15) <- c('study_id', 'AUC15')

Baseline30 <- MMTT %>% filter(visit == 'DIVIDInt Visit 2: Baseline') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +30 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Baseline30) <- c('study_id', 'AUC30')

Baseline60 <- MMTT %>% filter(visit == 'DIVIDInt Visit 2: Baseline') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +60 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Baseline60) <- c('study_id', 'AUC60')

Baseline90 <- MMTT %>% filter(visit == 'DIVIDInt Visit 2: Baseline') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +90 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Baseline90) <- c('study_id', 'AUC90')

Baseline120 <- MMTT %>% filter(visit == 'DIVIDInt Visit 2: Baseline') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +120 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Baseline120) <- c('study_id', 'AUC120')

Baseline <- full_join(Baseline0, Baseline15)
Baseline <- full_join(Baseline, Baseline30)
Baseline <- full_join(Baseline, Baseline60)
Baseline <- full_join(Baseline, Baseline90)
Baseline <- full_join(Baseline, Baseline120)

Baseline <- Baseline %>% mutate(AUCTrap = ((AUC0 + AUC15)/2*15 +
                     (AUC15 + AUC30)/2*15 +
                     (AUC30 + AUC60)/2*30 + 
                     (AUC60+ AUC90)/2*30 + 
                     (AUC90 + AUC120)/2*30)/120000)
Baseline$AUC <- log(Baseline$AUCTrap + 1)


# some have missing timepoints: 
# calculate AUC with available timepoints (?)


weights <- c(7.5, 15, 22.5, 30, 30, 15)
Baseline$AUCimputed = apply(Baseline[, 2:7], 
                            MARGIN=1, FUN=weighted.mean,
                            w=weights, na.rm =TRUE)/1000

# Check they are the same for non missing
#Baseline %>% filter(!is.na(AUCTrap)) %>% mutate(Diff = AUCTrap-AUCimputed) 
# check if the missing ones are correct
#Baseline %>% filter(is.na(AUCTrap))


Baseline$AUCimp <- log(Baseline$AUCimputed +1)


# calculate AUC for one year visit

OneYear0 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide 0 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(OneYear0) <- c('study_id', 'AUC0')

OneYear15 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +15 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(OneYear15) <- c('study_id', 'AUC15')

OneYear30 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +30 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(OneYear30) <- c('study_id', 'AUC30')

OneYear60 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +60 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(OneYear60) <- c('study_id', 'AUC60')

OneYear90 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +90 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(OneYear90) <- c('study_id', 'AUC90')

OneYear120 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +120 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(OneYear120) <- c('study_id', 'AUC120')

OneYear <- full_join(OneYear0, OneYear15)
OneYear <- full_join(OneYear, OneYear30)
OneYear <- full_join(OneYear, OneYear60)
OneYear <- full_join(OneYear, OneYear90)
OneYear <- full_join(OneYear, OneYear120)

OneYear <- OneYear %>% mutate(AUCTrap = ((AUC0 + AUC15)/2*15 +
                                           (AUC15 + AUC30)/2*15 +
                                           (AUC30 + AUC60)/2*30 + 
                                           (AUC60+ AUC90)/2*30 + 
                                           (AUC90 + AUC120)/2*30)/120000)

OneYear$AUC <- log(OneYear$AUCTrap + 1)



# add those with extra visit?
Extra <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP EXTRA VISIT') 
#Extra$study_id
# two participants
# missing in OneYear: need to add


Extra0 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP EXTRA VISIT') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide 0 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Extra0) <- c('study_id', 'AUC0')

Extra15 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP EXTRA VISIT') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +15 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Extra15) <- c('study_id', 'AUC15')

Extra30 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP EXTRA VISIT') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +30 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Extra30) <- c('study_id', 'AUC30')

Extra60 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP EXTRA VISIT') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +60 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Extra60) <- c('study_id', 'AUC60')

Extra90 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP EXTRA VISIT') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +90 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Extra90) <- c('study_id', 'AUC90')

Extra120 <- MMTT %>% filter(visit == 'DIVIDInt Visit 7: 1 year from first IMP EXTRA VISIT') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +120 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(Extra120) <- c('study_id', 'AUC120')

Extra <- full_join(Extra0, Extra15)
Extra <- full_join(Extra, Extra30)
Extra <- full_join(Extra, Extra60)
Extra <- full_join(Extra, Extra90)
Extra <- full_join(Extra, Extra120)

Extra <- Extra %>% mutate(AUCTrap = ((AUC0 + AUC15)/2*15 +
                                           (AUC15 + AUC30)/2*15 +
                                           (AUC30 + AUC60)/2*30 + 
                                           (AUC60+ AUC90)/2*30 + 
                                           (AUC90 + AUC120)/2*30)/120000)

Extra$AUC <- log(Extra$AUCTrap + 1)

OneYear <- rbind(OneYear, Extra)

# some have missing timepoints: 
# calculate AUC with available timepoints (?)


weights <- c(7.5, 15, 22.5, 30, 30, 15)
OneYear$AUCimputed = apply(OneYear[, 2:7], 
                           MARGIN=1, FUN=weighted.mean,
                           w=weights, na.rm =TRUE)/1000

# Check they are the same for non missing
#OneYear %>% filter(!is.na(AUCTrap)) %>% mutate(Diff = AUCTrap-AUCimputed) 
OneYear$AUCimp <- log(OneYear$AUCimputed +1)


# add AUC after 3, 6 months and 2 and 3 years

# calculate AUC for three month visit
ThreeM0 <- MMTT %>% filter(visit == 'DIVIDInt Visit 5: 3 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide 0 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeM0) <- c('study_id', 'AUC0')

ThreeM15 <- MMTT %>% filter(visit == 'DIVIDInt Visit 5: 3 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +15 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeM15) <- c('study_id', 'AUC15')

ThreeM30 <- MMTT %>% filter(visit == 'DIVIDInt Visit 5: 3 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +30 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeM30) <- c('study_id', 'AUC30')

ThreeM60 <- MMTT %>% filter(visit == 'DIVIDInt Visit 5: 3 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +60 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeM60) <- c('study_id', 'AUC60')

ThreeM90 <- MMTT %>% filter(visit == 'DIVIDInt Visit 5: 3 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +90 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeM90) <- c('study_id', 'AUC90')

ThreeM120 <- MMTT %>% filter(visit == 'DIVIDInt Visit 5: 3 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +120 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeM120) <- c('study_id', 'AUC120')


ThreeM <- full_join(ThreeM0, ThreeM15)
ThreeM <- full_join(ThreeM, ThreeM30)
ThreeM <- full_join(ThreeM, ThreeM60)
ThreeM <- full_join(ThreeM, ThreeM90)
ThreeM <- full_join(ThreeM, ThreeM120)

ThreeM <- ThreeM %>% mutate(AUCTrap = ((AUC0 + AUC15)/2*15 +
                                         (AUC15 + AUC30)/2*15 +
                                         (AUC30 + AUC60)/2*30 + 
                                         (AUC60+ AUC90)/2*30 + 
                                         (AUC90 + AUC120)/2*30)/120000)

ThreeM$AUC <- log(ThreeM$AUCTrap + 1)

# some have missing timepoints: 
# calculate AUC with available timepoints (?)


weights <- c(7.5, 15, 22.5, 30, 30, 15)
ThreeM$AUCimputed = apply(ThreeM[, 2:7], 
                           MARGIN=1, FUN=weighted.mean,
                           w=weights, na.rm =TRUE)/1000

# Check they are the same for non missing
#ThreeM %>% filter(!is.na(AUCTrap)) %>% mutate(Diff = AUCTrap-AUCimputed) 
ThreeM$AUCimp <- log(ThreeM$AUCimputed +1)


# 6 months

SixM0 <- MMTT %>% filter(visit == 'DIVIDInt Visit 6: 6 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide 0 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(SixM0) <- c('study_id', 'AUC0')

SixM15 <- MMTT %>% filter(visit == 'DIVIDInt Visit 6: 6 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +15 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(SixM15) <- c('study_id', 'AUC15')

SixM30 <- MMTT %>% filter(visit == 'DIVIDInt Visit 6: 6 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +30 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(SixM30) <- c('study_id', 'AUC30')

SixM60 <- MMTT %>% filter(visit == 'DIVIDInt Visit 6: 6 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +60 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(SixM60) <- c('study_id', 'AUC60')

SixM90 <- MMTT %>% filter(visit == 'DIVIDInt Visit 6: 6 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +90 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(SixM90) <- c('study_id', 'AUC90')

SixM120 <- MMTT %>% filter(visit == 'DIVIDInt Visit 6: 6 months from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +120 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(SixM120) <- c('study_id', 'AUC120')

SixM <- full_join(SixM0, SixM15)
SixM <- full_join(SixM, SixM30)
SixM <- full_join(SixM, SixM60)
SixM <- full_join(SixM, SixM90)
SixM <- full_join(SixM, SixM120)

SixM <- SixM %>% mutate(AUCTrap = ((AUC0 + AUC15)/2*15 +
                                     (AUC15 + AUC30)/2*15 +
                                     (AUC30 + AUC60)/2*30 + 
                                     (AUC60+ AUC90)/2*30 + 
                                     (AUC90 + AUC120)/2*30)/120000)

SixM$AUC <- log(SixM$AUCTrap + 1)

# some have missing timepoints: 
# calculate AUC with available timepoints (?)


weights <- c(7.5, 15, 22.5, 30, 30, 15)
SixM$AUCimputed = apply(SixM[, 2:7], 
                          MARGIN=1, FUN=weighted.mean,
                          w=weights, na.rm =TRUE)/1000

# Check they are the same for non missing
#SixM %>% filter(!is.na(AUCTrap)) %>% mutate(Diff = AUCTrap-AUCimputed) 
SixM$AUCimp <- log(SixM$AUCimputed +1)

# 2 years

TwoY0 <- MMTT %>% filter(visit == 'DIVIDInt Visit 8: 2 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide 0 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(TwoY0) <- c('study_id', 'AUC0')

TwoY15 <- MMTT %>% filter(visit == 'DIVIDInt Visit 8: 2 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +15 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(TwoY15) <- c('study_id', 'AUC15')

TwoY30 <- MMTT %>% filter(visit == 'DIVIDInt Visit 8: 2 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +30 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(TwoY30) <- c('study_id', 'AUC30')

TwoY60 <- MMTT %>% filter(visit == 'DIVIDInt Visit 8: 2 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +60 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(TwoY60) <- c('study_id', 'AUC60')

TwoY90 <- MMTT %>% filter(visit == 'DIVIDInt Visit 8: 2 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +90 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(TwoY90) <- c('study_id', 'AUC90')

TwoY120 <- MMTT %>% filter(visit == 'DIVIDInt Visit 8: 2 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +120 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(TwoY120) <- c('study_id', 'AUC120')

TwoY <- full_join(TwoY0, TwoY15)
TwoY <- full_join(TwoY, TwoY30)
TwoY <- full_join(TwoY, TwoY60)
TwoY <- full_join(TwoY, TwoY90)
TwoY <- full_join(TwoY, TwoY120)

TwoY <- TwoY %>% mutate(AUCTrap = ((AUC0 + AUC15)/2*15 +
                                     (AUC15 + AUC30)/2*15 +
                                     (AUC30 + AUC60)/2*30 + 
                                     (AUC60+ AUC90)/2*30 + 
                                     (AUC90 + AUC120)/2*30)/120000)

TwoY$AUC <- log(TwoY$AUCTrap + 1)



# some have missing timepoints: 
# calculate AUC with available timepoints (?)


weights <- c(7.5, 15, 22.5, 30, 30, 15)
TwoY$AUCimputed = apply(TwoY[, 2:7], 
                        MARGIN=1, FUN=weighted.mean,
                        w=weights, na.rm =TRUE)/1000

# Check they are the same for non missing
#TwoY %>% filter(!is.na(AUCTrap)) %>% mutate(Diff = AUCTrap-AUCimputed) 
TwoY$AUCimp <- log(TwoY$AUCimputed +1)

# 3 years
# still empty

ThreeY0 <- MMTT %>% filter(visit == 'DIVIDInt Visit 9: 3 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide 0 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeY0) <- c('study_id', 'AUC0')

ThreeY15 <- MMTT %>% filter(visit == 'DIVIDInt Visit 9: 3 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +15 sample Aliquot')%>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeY15) <- c('study_id', 'AUC15')

ThreeY30 <- MMTT %>% filter(visit == 'DIVIDInt Visit 9: 3 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +30 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeY30) <- c('study_id', 'AUC30')

ThreeY60 <- MMTT %>% filter(visit == 'DIVIDInt Visit 9: 3 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +60 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeY60) <- c('study_id', 'AUC60')

ThreeY90 <- MMTT %>% filter(visit == 'DIVIDInt Visit 9: 3 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +90 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeY90) <- c('study_id', 'AUC90')

ThreeY120 <- MMTT %>% filter(visit == 'DIVIDInt Visit 9: 3 years from first IMP') %>% 
  filter(sample_type_c_pep == 'MMTT C-peptide +120 sample Aliquot') %>%  dplyr::select(study_id, result_c_pep_mmtt)
names(ThreeY120) <- c('study_id', 'AUC120')

ThreeY <- full_join(ThreeY0, ThreeY15)
ThreeY <- full_join(ThreeY, ThreeY30)
ThreeY <- full_join(ThreeY, ThreeY60)
ThreeY <- full_join(ThreeY, ThreeY90)
ThreeY <- full_join(ThreeY, ThreeY120)

ThreeY <- ThreeY %>% mutate(AUCTrap = ((AUC0 + AUC15)/2*15 +
                                         (AUC15 + AUC30)/2*15 +
                                         (AUC30 + AUC60)/2*30 + 
                                         (AUC60+ AUC90)/2*30 + 
                                         (AUC90 + AUC120)/2*30)/120000)

ThreeY$AUC <- log(ThreeY$AUCTrap + 1)

# some have missing timepoints: 
# calculate AUC with available timepoints (?)


weights <- c(7.5, 15, 22.5, 30, 30, 15)
ThreeY$AUCimputed = apply(ThreeY[, 2:7], 
                        MARGIN=1, FUN=weighted.mean,
                        w=weights, na.rm =TRUE)/1000

# Check they are the same for non missing
#ThreeY %>% filter(!is.na(AUCTrap)) %>% mutate(Diff = AUCTrap-AUCimputed) 
ThreeY$AUCimp <- log(ThreeY$AUCimputed +1)





B <- Baseline %>%  dplyr::select(study_id, AUC, AUCimp) 
Y <- OneYear %>%  dplyr::select(study_id, AUC, AUCimp) 
M3 <- ThreeM %>%  dplyr::select(study_id, AUC, AUCimp)
M6 <- SixM %>%  dplyr::select(study_id, AUC, AUCimp)
Y2 <- TwoY %>%  dplyr::select(study_id, AUC, AUCimp)
Y3 <- ThreeY %>%  dplyr::select(study_id, AUC,  AUCimp)

ChangeAUC <- rbind(B, Y, M3, M6, Y2, Y3)

# check if we use time like this (in months) or more precise 
# (eg days exactly from visit db)
ChangeAUC$time <- c(rep(0, nrow(B)), 
                    rep(12, nrow(Y)), 
                    rep(3, nrow(M3)), 
                    rep(6, nrow(M6)), 
                    rep(24, nrow(Y2)), 
                    rep(36, nrow(Y3)))

ChangeAUC$timefactor <- c(rep('0', nrow(B)), 
                    rep('3', nrow(Y)), 
                    rep('1', nrow(M3)), 
                    rep('2', nrow(M6)), 
                    rep('4', nrow(Y2)), 
                    rep('5', nrow(Y3)))
#names(B) <- c('study_id', 'baseline')
#ChangeAUC <- left_join(ChangeAUC, B)
  
  
# add covariates info
ADSL <- read.table('td/ADSL.csv')

# Merge all into efficacy analysis set
EAB <- left_join(ADSL, ChangeAUC)
write.table(EAB, 'ad/EAB.csv')


# add secondary outcomes: 


# Peak AUC > 200mmol

Baseline$max <- apply(Baseline[, 2:7], 
                  MARGIN=1, FUN=max,
                  na.rm =TRUE)

Baseline$peak <- as.numeric(Baseline$max > 200)


ThreeM$max <- apply(ThreeM[, 2:7], 
                      MARGIN=1, FUN=max,
                      na.rm =TRUE)
ThreeM$peak <- as.numeric(ThreeM$max > 200)



SixM$max <- apply(SixM[, 2:7], 
                    MARGIN=1, FUN=max,
                    na.rm =TRUE)


SixM$peak <- as.numeric(SixM$max > 200)


OneYear$max <- apply(OneYear[, 2:7], 
                    MARGIN=1, FUN=max,
                    na.rm =TRUE)


OneYear$peak <- as.numeric(OneYear$max > 200)


TwoY$max <- apply(TwoY[, 2:7], 
                     MARGIN=1, FUN=max,
                     na.rm =TRUE)


TwoY$peak <- as.numeric(TwoY$max > 200)

ThreeY$max <- apply(ThreeY[, 2:7], 
                     MARGIN=1, FUN=max,
                     na.rm =TRUE)

ThreeY$peak <- as.numeric(ThreeY$max > 200)


BP <- Baseline %>%  dplyr::select(study_id, max, peak) 
YP <- OneYear %>%  dplyr::select(study_id, max, peak) 
M3P <- ThreeM %>%  dplyr::select(study_id, max, peak)
M6P <- SixM %>%  dplyr::select(study_id, max, peak)
Y2P <- TwoY %>%  dplyr::select(study_id, max, peak)
Y3P <- ThreeY %>%  dplyr::select(study_id, max, peak)

PeakAUC <- rbind(BP, YP, M3P, M6P, Y2P, Y3P)

# check if we use time like this (in months) or more precise 
# (eg days exactly from visit db)
PeakAUC$time <- c(rep(0, nrow(B)), 
                    rep(12, nrow(Y)), 
                    rep(3, nrow(M3)), 
                    rep(6, nrow(M6)), 
                    rep(24, nrow(Y2)), 
                    rep(36, nrow(Y3)))

PeakAUC$timefactor <- c(rep('0', nrow(B)), 
                          rep('3', nrow(Y)), 
                          rep('1', nrow(M3)), 
                          rep('2', nrow(M6)), 
                          rep('4', nrow(Y2)), 
                          rep('5', nrow(Y3)))

ADSL <- read.table('td/ADSL.csv')

# Merge all into efficacy analysis set
EAB_peak <- left_join(ADSL, PeakAUC)
write.table(EAB_peak, 'ad/EAB_peak.csv')



###	Mean insulin dosage per kilo bodyweight at each visit
###	Number of severe and less severe hypoglycaemic events at each visit
Ins <- data$DiViDint_DIABETES_CARE_2022.08.31.csv
Ins <- Ins %>% filter(study_id %in% id_list)

#head(Ins)
# we need average daily dose and we have to divide it by weight
# we have the counts of events
# again, we can not add this to EAB because different timings

InsOutcomes <- Ins %>%  dplyr::select(study_id, visit, average_daily_dose, hypoglycemia_events_since_last_visit, less_severe_hypoglycemia_events_since_last_visit)
INS <- inner_join(InsOutcomes, ADSL)

INS <- INS %>% mutate(average_daily_dose_per_weight = average_daily_dose/weight)

write.table(INS, 'ad/INS_secondary.csv')




###	Absolute HbA1c at each visit

HbA <- data$DiViDint_SAMPLING_2022.08.31.csv
HbA <- HbA %>% filter(study_id %in% id_list)
#head(HbA)

HbAOutcomes <- HbA %>%  dplyr::select(study_id, visit, HbA1c_result, albumin_result) 
HBA <- inner_join(HbAOutcomes, ADSL)



###	Change from baseline in HbA1c at each visit 
HbABaseline <- HbAOutcomes %>% filter(visit == 'DIVIDInt Visit 2: Baseline') %>%  dplyr::select(study_id, HbA1c_result, albumin_result)
names(HbABaseline) <- c('study_id', 'HbA1c_baseline', 'Alb_baseline')

HBAFull <- inner_join(HBA, HbABaseline)
HBAFull <- HBAFull %>% mutate(HbA1c_Change = HbA1c_result - HbA1c_baseline) %>% mutate(Alb_Change = albumin_result - Alb_baseline)

write.table(HBAFull, 'ad/Hb_secondary.csv')










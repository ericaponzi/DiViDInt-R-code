##############################
# this script contains the code to analyse the primary outcome
# the primary outcome is the AUC derived from the MMTT test
# calculations of the AUC from the raw data are reported in the 
# file "raw.to.ad.R"
# the data used here is the "EAB" datafile
# (already processed)
##############################



##############################
# load libraries
##############################
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(emmeans)
library(margins)
library(gridExtra)
library(haven)
library(tidyverse)
library(broom)
library(knitr)
library(GGally)
library(purrr)


# load the data
EAB <- read.csv("../data/ad/EAB.csv", sep = '')
str(EAB)

# the main outcome to be used in the analysis is
# AUCimp (imputed in the script "raw.to.ad")
# information of compliance is in the variables 
# Compliance_Plec and Compliance_Rib
# definitions of per protocol
# PP80, PP70, PP60, depending on the compliance threshold

##############################
# some summary statistics and
# data visualization
##############################

# this filters out the two patients who were included but 
# dropped out immediately
# (no outcome)
EAB <- EAB %>% filter(!is.na(timefactor))

# this produces a table with the number of records for each visit
# I report both the complete records and the ones I have imputed
No = function(x) {sum(!is.na(x))}
N_visits <- EAB %>% group_by(timefactor) %>% summarise_at('AUC', No)
N_visits_imputed <- EAB %>% group_by(timefactor) %>% summarise_at('AUCimp', No)
N <- EAB %>% group_by(timefactor) %>% summarise_at('study_id', No)

temp <- merge(N_visits, N_visits_imputed, by = 'timefactor', all =TRUE)
temp <- merge(N, temp, by = 'timefactor', all =TRUE)

rownames(temp) <- c('Baseline', 'Three Months', 'Six Months', 'One Year', 'Two Years')
temp


# This produces longitudinal graphs of 
# the imputed log-transformed AUC C-peptide plotted by visit (up to 12 months) 
# I split the plots into 7 for better visualization

id_list <- unique(EAB$study_id)
EABfull <- EAB %>% filter(!is.na(timefactor))

EAB1 <- EABfull[EABfull$study_id %in% id_list[1:14], ]
EAB2 <- EABfull[EABfull$study_id %in% id_list[15:28], ]
EAB3 <- EABfull[EABfull$study_id %in% id_list[29:42], ]
EAB4 <- EABfull[EABfull$study_id %in% id_list[43:56], ]
EAB5 <- EABfull[EABfull$study_id %in% id_list[57:70], ]
EAB6 <- EABfull[EABfull$study_id %in% id_list[71:84], ]
EAB7 <- EABfull[EABfull$study_id %in% id_list[85:98], ]

p1 <-  ggplot(EAB1, aes(x = timefactor, y = AUCimp, col = study_id, group = study_id)) + geom_point() +
  geom_line()+ theme_bw() +
  #   theme(legend.position = "none") +
  xlab('Visit')+  ggtitle('Patient 1 to 14') +
  theme(plot.title = element_text(size = 10))

p2 <-  ggplot(EAB2, aes(x = timefactor, y = AUCimp, col = study_id, group = study_id)) + geom_point() +
  geom_line() + theme_bw()+
  #  theme(legend.position = "none")+
  xlab('Visit')+  ggtitle('Patient 15 to 28') +
  theme(plot.title = element_text(size = 10))


p3 <-  ggplot(EAB3, aes(x = timefactor, y = AUCimp, col = study_id, group = study_id)) + geom_point() +
  geom_line()+ theme_bw() +
  #  theme(legend.position = "none")+
  xlab('Visit')+  ggtitle('Patient 29 to 42') +
  theme(plot.title = element_text(size = 10))


p4 <-  ggplot(EAB4, aes(x = timefactor, y = AUCimp, col = study_id, group = study_id)) + geom_point() +
  geom_line()+ theme_bw() +
  #  theme(legend.position = "none")+
  xlab('Visit')+  ggtitle('Patient 43 to 56') +
  theme(plot.title = element_text(size = 10))

p5 <-  ggplot(EAB5, aes(x = timefactor, y = AUCimp, col = study_id, group = study_id)) + geom_point() +
  geom_line()+ theme_bw() +
  #  theme(legend.position = "none")+
  xlab('Visit')+  ggtitle('Patient 57 to 70') +
  theme(plot.title = element_text(size = 10))

p6 <-  ggplot(EAB6, aes(x = timefactor, y = AUCimp, col = study_id, group = study_id)) + geom_point() +
  geom_line()+ theme_bw() +
  #  theme(legend.position = "none")+
  xlab('Visit')+  ggtitle('Patient 71 to 84') +
  theme(plot.title = element_text(size = 10))

p7 <-  ggplot(EAB7, aes(x = timefactor, y = AUCimp, col = study_id, group = study_id)) + geom_point() +
  geom_line()+ theme_bw() +
  #  theme(legend.position = "none")+
  xlab('Visit')+  ggtitle('Patient 85 to 98') +
  theme(plot.title = element_text(size = 10))



# the main outcome is defined as the AUC at 12 months
# this produces summary statistics for the primary outcome
# N, mean and sd are reported overall and per treatment group

No = function(x) {sum(!is.na(x))}
Mean = function(x) {mean(x, na.rm = TRUE)}
Sd = function(x) {sd(x, na.rm = TRUE)}
stats = c('Mean', 'Sd','No')


summary_AUC_overall <- EAB %>% filter(timefactor == 3) %>% summarise_at('AUCimp', stats)
summary_AUC_bytreat <- EAB %>% filter(timefactor == 3) %>% group_by(treatment) %>% summarise_at('AUCimp', stats)

N <- c(summary_AUC_overall$No, summary_AUC_bytreat$No)
Summary <- c(paste(round(summary_AUC_overall$Mean,2), '(', 
                   round(summary_AUC_overall$Sd,2), ')', sep =''), 
             paste(round(summary_AUC_bytreat$Mean,2), '(', 
                   round(summary_AUC_bytreat$Sd,2), ')', sep =''))
tabAUC <- data.frame(N, Summary)
summarytableAUC <- cbind(tabAUC[1, ], tabAUC[3, ], tabAUC[2, ])
rownames(summarytableAUC) <- c('AUC at 12 months')
summarytableAUC



##############################
# primary analysis
# linear mixed model
# on FAS population
##############################


# we use a linear mixed model, including the AUC C-peptide values
# at baseline, at 3 months, 6 months and 12 months as longitudinal outcome. 
# The treatment effect will be 0 at baseline as this is prior to randomization. 
# Fixed effects: time, center, treatment and interaction time-treatment. 
# Time is categorical (levels = visits)
# Random effects: patient ID. 
# all patients are included (FAS)

# change the time into factor
# and use levels so that the reference value 
# is the baseline
EAB$timefactor <- as.factor(EAB$timefactor)
EAB$timefactornew <- paste('T', EAB$timefactor, sep ='')


# create a new treatment variable that is 0 at baseline
# and same as the treatment at the other visits
EAB.2treat <- EAB %>%
  mutate(trt2 = if_else(time == 0,0,as.numeric(trt))) %>%
  mutate(trt2f = factor(trt2))%>%
  mutate(center = factor(site))

# lmer
model.primary <- lmer(AUC ~  site + age + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                      data = EAB.2treat, REML=FALSE)

# the warning on dropping one column is expected (interaction term time and trt does not make sense at baseline!)
# look at model results
# coefficients and p-values

summary(model.primary)

coefs <- cbind(coef(summary(model.primary)),
               confint(model.primary)[-c(1,2), ])

coefs <- as.data.frame(coefs)
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$`t value`)))
coefs

# the main effect is only at 12 months!
# we need to extract the marginal effect
ME <- model.primary %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))

ME
# what we report from here is the third row (timefactornew = 3)
# the AME is the treatment effect (beta in treatment - beta ctrl).

# Plot marginal effect: 
  

p1 <- 
  EAB %>%
  mutate(trt2 = treatment) %>%
  mutate(trt2 = fct_expand(trt2,"0")) %>%
  mutate(trt2 = if_else(time != 0,trt2,factor(0) )) %>%
  lmer(AUC ~  site + timefactornew* trt2 +  (1|study_id),data=., REML=FALSE)  %>%
  ref_grid() %>%
  emmeans( ~ trt2 | timefactornew) %>%
  tidy(conf.int = TRUE) %>%
  filter((trt2 != 0 & timefactornew != 'T0')   | (trt2==0 & timefactornew == 'T0')) %>%
  bind_rows(.,mutate(filter(.,trt2==0 & timefactornew=='T0'),trt2 = fct_recode(trt2,"Placebo"="0"))) %>%
  mutate(trt2= fct_recode(trt2,"Active" = "0"))  


p1 %>%
  ggplot(aes(timefactornew, estimate, color=trt2, group=trt2)) +
  geom_point(position = position_dodge(0.04)) +
  geom_line() + 
  xlab("Time") + 
  ylab("Estimate") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width=0.1,
                position = position_dodge(0.04)) +
  theme_classic() + 
  theme(legend.position=c(0.9,0.9)) +
  scale_colour_brewer(palette = "Set1", name="Treatment")

# check that the two lines are exactly the same at baseline!
# check that the AME above is the difference between 
# the blue and red dot at T3
a <- p1 %>% filter(timefactornew == 'T3') %>% filter(trt2 == 'Active') %>% select(estimate)
b <- p1 %>% filter(timefactornew == 'T3') %>% filter(trt2 == 'Placebo') %>% select(estimate)
a-b
ME[3, ]$AME

# some checks on the main assumptions for the lmer
# normality of residuals

qqnorm(residuals(model.primary))


# presence of outliers
# visual check

ggplot(EAB %>% filter(timefactor == 3), aes(y = AUCimp, color = treatment))+ geom_boxplot()+ theme_bw() 

# there seems to be one outlier
# check with boxplot.stats

EAB1 <- EAB %>% filter(timefactor == 3)
value.out <- boxplot.stats(EAB1$AUCimp)$out
index.out <- which(EAB1$AUCimp %in% c(value.out))

##############################
# sensitivity analysis
# on PP population
##############################


# different definition of compliance to treatment
# first case: 70% of the treatment is taken
# we use average compliance 
# ribavirin + pleconaril
# and no withdrawal

EAB.2treat.comp <- EAB.2treat %>% filter(!is.na(Compliance_Rib)) %>% 
  mutate(Comp = (Compliance_Plec + Compliance_Rib)/2)

EAB.2treat.comp <- EAB.2treat.comp %>% mutate(PPP = ifelse(withdrawn_ever == 1, 0, 1))


EAB.2treat.70 <- EAB.2treat.comp %>% 
  filter(Comp > 0.7) %>% filter(PPP == 1)

ncompliant.70 <- sum(!duplicated(EAB.2treat.70$study_id))

# fit the same model as above
model.70 <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                 data = EAB.2treat.70, REML=FALSE)

# coefficients and p-values
coefs <- cbind(coef(summary(model.70)),
               confint(model.70)[-c(1,2), ])

coefs <- as.data.frame(coefs)

coefs$p.z <- 2 * (1 - pnorm(abs(coefs$`t value`)))
coefs

# marginal effects
ME.70 <- model.70 %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))
ME.70


# second case: 60% of the treatment is taken
# we use average compliance 
# ribavirin + pleconaril
# and no withdrawal


EAB.2treat.60 <- EAB.2treat.comp %>% 
  filter(Comp > 0.6) %>% filter(PPP == 1)

ncompliant.60 <- sum(!duplicated(EAB.2treat.60$study_id))

# fit the same model as above
model.60 <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                 data = EAB.2treat.60, REML=FALSE)

# coefficients and p-values
coefs <- cbind(coef(summary(model.60)),
               confint(model.60)[-c(1,2), ])

coefs <- as.data.frame(coefs)
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$`t value`)))

# marginal effects
ME.60 <- model.60 %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))

ME.60


# third case: 80% of the treatment is taken
# we use average compliance 
# ribavirin + pleconaril
# and no withdrawal

EAB.2treat.80 <- EAB.2treat.comp %>% 
  filter(Comp > 0.8) %>% filter(PPP == 1)

ncompliant.80 <- sum(!duplicated(EAB.2treat.80$study_id))

# run the same model as above
model.80 <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                 data = EAB.2treat.80, REML=FALSE)

# coefficients and pvalues
coefs <- cbind(coef(summary(model.80)),
               confint(model.80)[-c(1,2), ])

coefs <- as.data.frame(coefs)
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$`t value`)))

# marginal effects
ME.80 <- model.80 %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))

ME.80

##############################
# sensitivity analysis
# on outliers
##############################

# identify one outlier
# see above assumption checks
out <- EAB1[index.out, ]$study_id

# run the primary model on set without outlier
EAB.2treat.nooutliers <- EAB.2treat %>% filter(!is.element(study_id, out))
model.noout <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                    data = EAB.2treat.nooutliers, REML=FALSE)

# coefficients and p-values
coefs <- cbind(coef(summary(model.noout)),
               confint(model.noout)[-c(1,2), ])

coefs <- as.data.frame(coefs)
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$`t value`)))

# marginal effects
ME.noout <- model.noout %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))
ME.noout


##############################
# sensitivity analysis
# on missing data
##############################

# two missing data imputation methods are used
# 1. last observation carried forward (the last available AUC is used for the following visits)
# 2. worst observation (the lowest available AUC is used for the missing visits)


# LOCF

# create dataset for imputation
id <- c(rep(as.character(id_list),5))
visit <- rep(c(0, 3, 6, 12, 24), each = 96)
#timefactor <- paste('T', visit, sep ='')

dd_empty <- data.frame(study_id = id, time = visit)
ADSL = read.table(file="../data/td/ADSL.csv",sep="")
dd_empty <- full_join(dd_empty, ADSL)



dd0 <- full_join(dd_empty, EAB, by = c('study_id', 'time'))


dd0 <- dd0 %>%
  mutate(trt2 = if_else(timefactor == 0,0,as.numeric(trt.x))) %>%
  mutate(trt2f = factor(trt2))%>%
  mutate(center = factor(site.x))

dd0 = dd0 %>% filter(!is.na(time))
#dd0.1 = dd0 %>% filter(time !=  24)
dd = dd0  %>% group_by(study_id) %>% arrange(time) %>% 
  fill(AUCimp, .direction ='downup')



# some checks on imputation here: 
anymiss <- dd0 %>% filter(is.na(AUCimp)) %>% select(study_id)
anymiss.id <- unique(anymiss)

dd0 %>% filter(study_id == as.character(anymiss.id[18,1])) %>% select(AUCimp, time)
dd %>% filter(study_id == as.character(anymiss.id[18,1])) %>% select(AUCimp, time)
EAB%>% filter(study_id == as.character(anymiss.id[18,1]))
# this copied T1 into baseline (4 patient)


model.lf <- lmer(AUC ~  center + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                 data = dd, REML=FALSE)

# coefficients and pvalues
coefs <- cbind(coef(summary(model.lf)),
               confint(model.lf)[-c(1,2), ])

coefs <- as.data.frame(coefs)
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$`t value`)))

# marginal effects

ME.lf <- model.lf %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))
ME.lf


# WCI


dd2 = dd0  %>% group_by(study_id) %>% arrange(desc(AUCimp)) %>% 
  fill(AUCimp, .direction ='downup')


# some checks on imputation here: 
anymiss <- dd0 %>% filter(is.na(AUCimp)) %>% select(study_id)
anymiss.id <- unique(anymiss)

dd0 %>% filter(study_id == as.character(anymiss.id[34,1])) %>% select(AUCimp, time)
dd2 %>% filter(study_id == as.character(anymiss.id[34,1])) %>% select(AUCimp, time)




model.wc <- lmer(AUC ~  center + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                 data = dd2, REML=FALSE)

# coefficients and pvalues
coefs <- cbind(coef(summary(model.wc)),
               confint(model.wc)[-c(1,2), ])

coefs <- as.data.frame(coefs)
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$`t value`)))

# marginal effects
ME.wc <- model.wc %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))

ME.wc


##############################
# sensitivity analysis
# using Wilcoxon test
# on the change from baseline
##############################

# in this analysis I only use the 12 month measure 
# and compute the change from baseline
# then run a wilcoxon test on the change from baseline

baseline_AUC <- EAB %>% filter(timefactor == 0) %>% dplyr:::select('study_id', 'AUCimp')
names(baseline_AUC) <- c('study_id', 'AUC0')

EAB.change <- left_join(EAB, baseline_AUC)
EAB.change <- EAB.change %>% mutate(AUCchange = AUCimp-AUC0)


EAB.change.12 <- EAB.change %>% filter(timefactor == 3)

res <- wilcox.test(AUCchange ~ treatment, data = EAB.change.12,
                   exact = FALSE)

res


##############################
# subgroup analysis
##############################


# 1. compliance, above/below 70: 
# note that some do not have information on compliance and are here excluded from both models. 


df.above <- EAB.2treat.comp %>% filter(Comp > 0.7)
df.below <- EAB.2treat.comp %>% filter(Comp <= 0.7)

model.above <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                    data = df.above, REML=FALSE)

model.below <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                    data = df.below, REML=FALSE)

# Coefficients and p-values: 
  
coefs.a <- cbind(coef(summary(model.above)),
                 confint(model.above)[-c(1,2), ])

coefs.a <- as.data.frame(coefs.a)
coefs.a$p.z <- 2 * (1 - pnorm(abs(coefs.a$`t value`)))

coefs.b <- cbind(coef(summary(model.below)),
                 confint(model.below)[-c(1,2), ])

coefs.b <- as.data.frame(coefs.b)
coefs.b$p.z <- 2 * (1 - pnorm(abs(coefs.b$`t value`)))

coefs <- cbind(coefs.a, coefs.b)
coefs 

# marginal effects 

ME.a <- model.above %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))

ME.b <- model.below %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))



ME.a
ME.b


# 2. presence of enterovirus at baseline (present vs. not present)

# Not available yet.



# 3. included before/after the covid pandemic
# I defined the two groups as before/after march 2020

dates.db <- read.csv('../data/td/dates.csv', sep = ' ')

dates.db <- dates.db %>% mutate(year_inc = substr(as.character(date_inclusion), 7, 10)) %>% 
  mutate(month_inc = substr(as.character(date_inclusion), 4, 5))


dates.db <- dates.db %>% mutate(year_inc = as.numeric(year_inc)) %>%
  mutate(month_inc = as.numeric(month_inc)) %>% mutate(after = (month_inc > 2 & year_inc > 2019))

id.before <- dates.db %>% filter(after == FALSE) %>% select(study_id)
id.after <- dates.db %>% filter(after == TRUE) %>% select(study_id)

id.b <- as.character(id.before[,1])
id.a <- as.character(id.after[,1])

df.before <- EAB.2treat %>% filter(study_id %in% id.b)
df.after <- EAB.2treat %>% filter(study_id %in% id.a)

model.before <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                     data = df.before, REML=FALSE)

model.after <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                    data = df.after, REML=FALSE)

# Coefficients and p-values: 
  
coefs.a <- cbind(coef(summary(model.after)),
                 confint(model.after)[-c(1,2), ])

coefs.a <- as.data.frame(coefs.a)
coefs.a$p.z <- 2 * (1 - pnorm(abs(coefs.a$`t value`)))

coefs.b <- cbind(coef(summary(model.before)),
                 confint(model.before)[-c(1,2), ])

coefs.b <- as.data.frame(coefs.b)
coefs.b$p.z <- 2 * (1 - pnorm(abs(coefs.b$`t value`)))

coefs <- merge(coefs.b,coefs.a,by="row.names",all.x=TRUE)
coefs


# marginal effects
# note that there is no T4 for those included after the pandemic:
  

ME.a <- model.after %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))

ME.b <- model.before %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))



ME.a
ME.b


# 4. age at inclusion (below/above the age of 12)

df.above <- EAB.2treat %>% filter(age > 12)
df.below <- EAB.2treat %>% filter(age <= 12)

model.above <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                    data = df.above, REML=FALSE)

model.below <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                    data = df.below, REML=FALSE)

# Coefficients and p-values: 

coefs.a <- cbind(coef(summary(model.above)),
                 confint(model.above)[-c(1,2), ])

coefs.a <- as.data.frame(coefs.a)
coefs.a$p.z <- 2 * (1 - pnorm(abs(coefs.a$`t value`)))

coefs.b <- cbind(coef(summary(model.below)),
                 confint(model.below)[-c(1,2), ])

coefs.b <- as.data.frame(coefs.b)
coefs.b$p.z <- 2 * (1 - pnorm(abs(coefs.b$`t value`)))

coefs <- cbind(coefs.a, coefs.b)
coefs 

# marginal effects 

ME.a <- model.above %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))

ME.b <- model.below %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))



ME.a
ME.b



# 5. sex 

df.male <- EAB.2treat %>% filter(gender == 'MALE')
df.female <- EAB.2treat %>% filter(gender == 'FEMALE')

model.m <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                    data = df.male, REML=FALSE)

model.f <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                    data = df.female, REML=FALSE)

# Coefficients and p-values: 

coefs.a <- cbind(coef(summary(model.m)),
                 confint(model.m)[-c(1,2), ])

coefs.a <- as.data.frame(coefs.a)
coefs.a$p.z <- 2 * (1 - pnorm(abs(coefs.a$`t value`)))

coefs.b <- cbind(coef(summary(model.f)),
                 confint(model.f)[-c(1,2), ])

coefs.b <- as.data.frame(coefs.b)
coefs.b$p.z <- 2 * (1 - pnorm(abs(coefs.b$`t value`)))

coefs <- cbind(coefs.a, coefs.b)
coefs 

# marginal effects 

ME.a <- model.m %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))

ME.b <- model.f %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))



ME.a
ME.b

# 6. insulin dosage U/kg at baseline (below/above 0,5 U/kg)
INS <- read.csv('../data/ad/INS_secondary.csv', sep = '')
INS <- INS %>% filter(visit == 'DIVIDInt Visit 2: Baseline')

INS.above <- INS %>% filter(average_daily_dose_per_weight > 0.5)
INS.below <- INS %>% filter(average_daily_dose_per_weight <= 0.5)


id.above <- INS.above %>% select(study_id)
id.below <- INS.below %>% select(study_id)

id.a <- as.character(id.above[,1])
id.b <- as.character(id.below[,1])



df.below <- EAB.2treat %>% filter(study_id %in% id.b)
df.above <- EAB.2treat %>% filter(study_id %in% id.a)

model.below <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                     data = df.below, REML=FALSE)

model.above <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                    data = df.above, REML=FALSE)


# Coefficients and p-values: 

coefs.a <- cbind(coef(summary(model.above)),
                 confint(model.above)[-c(1,2), ])

coefs.a <- as.data.frame(coefs.a)
coefs.a$p.z <- 2 * (1 - pnorm(abs(coefs.a$`t value`)))

coefs.b <- cbind(coef(summary(model.below)),
                 confint(model.below)[-c(1,2), ])

coefs.b <- as.data.frame(coefs.b)
coefs.b$p.z <- 2 * (1 - pnorm(abs(coefs.b$`t value`)))

coefs <- cbind(coefs.a, coefs.b)
coefs 

# marginal effects 

ME.a <- model.above %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))

ME.b <- model.below %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))



ME.a
ME.b

# 7. HbA1c at baseline (below/above median) 


HBA <- read.csv('../data/ad/Hb_secondary.csv', sep = '')



HBA <- HBA %>% filter(visit == 'DIVIDInt Visit 2: Baseline')

m <- median(HBA$HbA1c_baseline, na.rm = TRUE)
HBA.above <- HBA %>% filter(HbA1c_baseline > m)
HBA.below <- HBA %>% filter(HbA1c_baseline <= m)


id.above <- HBA.above %>% select(study_id)
id.below <- HBA.below %>% select(study_id)

id.a <- as.character(id.above[,1])
id.b <- as.character(id.below[,1])



df.below <- EAB.2treat %>% filter(study_id %in% id.b)
df.above <- EAB.2treat %>% filter(study_id %in% id.a)

model.below <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                    data = df.below, REML=FALSE)

model.above <- lmer(AUC ~  site + trt2f + timefactornew + timefactornew* trt2f + (1|study_id), 
                    data = df.above, REML=FALSE)


# Coefficients and p-values: 

coefs.a <- cbind(coef(summary(model.above)),
                 confint(model.above)[-c(1,2), ])

coefs.a <- as.data.frame(coefs.a)
coefs.a$p.z <- 2 * (1 - pnorm(abs(coefs.a$`t value`)))

coefs.b <- cbind(coef(summary(model.below)),
                 confint(model.below)[-c(1,2), ])

coefs.b <- as.data.frame(coefs.b)
coefs.b$p.z <- 2 * (1 - pnorm(abs(coefs.b$`t value`)))

coefs <- cbind(coefs.a, coefs.b)
coefs 

# marginal effects 

ME.a <- model.above %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))

ME.b <- model.below %>%
  margins(model = .,variables = "trt2f", at= list(timefactornew=c("T1", "T2", "T3", "T4")), vce = "delta") %>%
  summary() %>%
  filter(!is.na(z)) %>% mutate_if(is.numeric, ~round(., 3))



ME.a
ME.b

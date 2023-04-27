## this script produces summary tables for adverse events
## in dividint 

#############################################
# load libraries 
library(dplyr)
library(tidyr)
library(readxl)

# load data
# open adverse event data
AE <- read.csv("../data/locked/DiViDint_Adverse_events_2022.08.22.csv", sep = ';')
#head(AE)

# merge with treatment allocation info 

# treatment_alloc <- read.table(file="../data/td/TreatmentAllocFake.csv",sep="")
AE <- inner_join(AE, treatment_alloc)

total_n <- 96

AE$se <- ifelse(AE$serious_adverse_event == 'YES', 1, 0)
AE$visit <- gsub("v", "V", AE$Reported.on.visit.)#reported or registered? 
#table(AE$visit)
#AE[which(AE$visit == 'Visit  4'), ]$visit <- 'Visit 4'
#AE[which(AE$visit == 'Visit  5'), ]$visit <- 'Visit 5'
AE$gr <- as.factor(AE$Grade..nr.)

# create AE tables

AE <- AE %>%
  mutate(anyae = if_else(is.na(ae_ctcae_v_4_0_term), 0, 1),
         sae = if_else(is.na(ae_ctcae_v_4_0_term), 0, se)
  ) %>%   group_by(visit, study_id) %>% 
  mutate(n_ae = sum(anyae),
         one_ae = n_ae == 1,
         two_ae = n_ae == 2,
         three_plus_ae = n_ae > 2,
         anysae = max(sae)) %>% 
  ungroup



tab0 <- AE %>% group_by(study_id, visit, treatment) %>% 
  summarise(n_aes = sum(anyae),
         one_ae = max(one_ae),
         two_ae = max(two_ae),
         three_plus_ae = max(three_plus_ae),
         anysae = sum(anysae),
        )


tab <- tab0 %>% group_by(visit, treatment) %>% 
  summarise(n_aes = sum(n_aes),
            one_ae = sum(one_ae),
            two_ae = sum(two_ae),
            three_plus_ae = sum(three_plus_ae),
            anysae = sum(anysae),
            n_pat = sum(!duplicated(study_id)), 
            pct = round(n_pat/total_n*100,digits = 2))
  


tab_SAE <- AE %>% filter(sae == 1) %>% group_by(visit, treatment) %>% 
  summarise(n_pat = sum(!duplicated(study_id)), 
            pct = round(n_pat/total_n*100,digits = 2))



s <- full_join(tab, tab_SAE, by = c('visit', 'treatment'))


#s <- t(s)
#rownames(s) = 



  
  

table <- list()

visits <- levels(as.factor(s$visit))
for (v in 1:length(visits)) {
  sv <- s %>% filter(visit == visits[v])
  table_t <- sv %>% filter(treatment == 'Pleconaril/Ribavirin')
  table_p <-  sv %>% filter(treatment == 'Placebo')
  n_aes <- c(table_t$n_aes, table_p$n_aes)
  n_pat <- c(table_t$n_pat.x, table_p$n_pat.x)
  perc_pat <- c(table_t$pct.x, table_p$pct.x)
  one_ae <- c(table_t$one_ae, table_p$one_ae)
  two_ae <- c(table_t$two_ae, table_p$two_ae)
  three_plus_ae <- c(table_t$three_plus_ae, table_p$three_plus_ae)
  any_sae <- c(table_t$anysae, table_p$anysae)
  n_pat_sae <- c(table_t$n_pat.y, table_p$n_pat.y)
  perc_pat_sae <- c(table_t$pct.y, table_p$pct.y)
  
  table[[v]] <- rbind(n_aes, n_pat, perc_pat,
                      one_ae, two_ae, three_plus_ae, 
                      any_sae, n_pat_sae, perc_pat_sae)
  rownames(table[[v]]) <- c('Number of AEs',  'Number of patients with any AEs', 'Percentage of patients with any AEs',
                            'Number of patients with one AE', 
                            'Number of patients with two AEs', 
                            'Number of patients with three or more AEs', 
                            'Number of SAEs',
                            'Number of patients with any SAEs',
                            'Percentage of patients with any SAEs')
  if (ncol(table[[v]]) == 2) colnames(table[[v]]) <- c('Pleconaril/Ribavirin', 'Placebo')
  else (colnames(table[[v]]) <- ifelse(nrow(table_t) == 1, 'Pleconaril/Ribavirin', 'Placebo'))
  
}



# cumulative up to 12 months

AE <- AE %>% filter(visit != 'Visit 8')
AE <- AE %>% filter(visit != 'Visit 9')

AE <- AE %>%
  mutate(anyae = if_else(is.na(ae_ctcae_v_4_0_term), 0, 1),
         sae = if_else(is.na(ae_ctcae_v_4_0_term), 0, se)
  ) %>%   group_by(visit, study_id) %>% 
  mutate(n_ae = sum(anyae),
         one_ae = n_ae == 1,
         two_ae = n_ae == 2,
         three_plus_ae = n_ae > 2,
         anysae = max(sae)) %>% 
  ungroup



tab0 <- AE %>% group_by(study_id, treatment) %>% 
  summarise(n_aes = sum(anyae),
            one_ae = max(one_ae),
            two_ae = max(two_ae),
            three_plus_ae = max(three_plus_ae),
            anysae = sum(anysae),
  )


tab <- tab0 %>% group_by(treatment) %>% 
  summarise(n_aes = sum(n_aes),
            one_ae = sum(one_ae),
            two_ae = sum(two_ae),
            three_plus_ae = sum(three_plus_ae),
            anysae = sum(anysae),
            n_pat = sum(!duplicated(study_id)), 
            pct = round(n_pat/total_n*100,digits = 2))



tab_SAE <- AE %>% filter(sae == 1) %>% group_by( treatment) %>% 
  summarise(n_pat = sum(!duplicated(study_id)), 
            pct = round(n_pat/total_n*100,digits = 2))



s <- full_join(tab, tab_SAE, by = c('treatment'))


# check hemolysis

hm <- AE %>% filter(ae_ctcae_v_4_0_term == 'Hemolysis')
hm2 <- AE %>% filter(ae_ctcae_v_4_0_term == 'hemolysis')

hmall <- rbind(hm, hm2)
hmall <- hmall[(which(!duplicated(hmall$study_id))), ]

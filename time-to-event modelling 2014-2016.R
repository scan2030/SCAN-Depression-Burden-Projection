library(dplyr)
library(lubridate)
library(tidyverse)
library(readxl)
library(fuzzyjoin)
library(survival)
library(survminer)
library(SurvRegCensCov)
library(flexsurv)

# 0.1 Steps ----
## 1. Identify all depression patients 2014-2016
## 2. Identify prevalent and incident patients in 2014-2016
## 3. Identify TRD development from incidence to December 2020
## 4. Identify age, sex and comorbidity status at incidence
## 5. Identify comorbidities acquisition after incident depression (for non-TRD) until Dec 2020
## 6. Identify post-TRD comorbidities acquisition after TRD (for TRD) until Dec 2020
## 7. Identify recorded deaths by cause by December 2020
## 8. Identify recovered patients who don't have any diagnosis after 5-year follow-up until Dec 2021

# 0.2 Read functions ----
calc_time <- function(refDate, endDate) {
    period <- as.period(interval(refDate, endDate), unit = "day")
    period$day + 1
}
calc_age<-function(birthDate, refDate) {
    period <- as.period(interval(birthDate, refDate), unit = "year")
    period$year
}
my_custom_name_repair <- function(x) tolower(gsub("\\.{1,}",'\\.',gsub("\n|  ", "", make.names(x))))
merge_files <- function(files){
    list_df <- lapply(files,readxl::read_xlsx,.name_repair=my_custom_name_repair) 
    df <- list_df %>% bind_rows(.)
    return(df)
}

## 1.1 Identify all depression patients from 2001 to 2019 ----
directory <- list.files("../../Depression/converted_data_dx/", full.names = T)
depression <- merge_files(directory)

# 1.2 Identify prevalent and incident patients in 2014-2016 ----
## Prevalent patients in 2014 = 15159
## Prevalent patients in 2015 = 15518
## Prevalent patients in 2016 = 15083
## Prevalent patients 2014-2016 = 41046
depression.14_16<-depression %>% 
    mutate(reference.date. = ymd(reference.date.),
           date.of.registered.death. = ymd(date.of.registered.death.),
           date.of.birth.yyyy.mm.dd. = as.Date(date.of.birth.yyyy.mm.dd.),
           date.of.birth.yyyy.mm.dd. = ymd(date.of.birth.yyyy.mm.dd.)) %>% 
    filter(reference.date.>= as.Date('2014-01-01') & reference.date.<= as.Date('2016-12-31')) %>% 
    arrange(reference.date.) %>% 
    distinct(reference.key., .keep_all=TRUE)

# Update death dates and causes of death of 2014-2016 patients to 2020 (current information is until 2019)
death.new.14.1<-read_excel('../../Depression/6531406_2014Depression_death-1_COM_conv_2020.xlsx')
death.new.14.2<-read_excel('../../Depression/6531407_2014Depression_death-2_COM_conv_2020.xlsx')
death.new.14<-rbind(death.new.14.1,death.new.14.2)

death.new.15.1<-read_excel('../../Depression/6531393_2015Depression_death-1_COM_conv_2020.xlsx')
death.new.15.2<-read_excel('../../Depression/6531395_2015Depression_death-2_COM_conv_2020.xlsx')
death.new.15<-rbind(death.new.15.1,death.new.15.2)

death.new.16.1<-read_excel('../../Depression/6531398_2016Depression_death-1_COM_conv_2020.xlsx')
death.new.16.2<-read_excel('../../Depression/6531396_2016Depression_death-2_COM_conv_2020.xlsx')
death.new.16<-rbind(death.new.16.1,death.new.16.2)

death.new.14_16<-rbind(death.new.14, death.new.15)
death.new.14_16<-rbind(death.new.14_16, death.new.16)
colnames(death.new.14_16)<-c('reference.key.', 'date.of.registered.death.', 'exact.date.of.death.','death.cause.main.cause.','death.cause.supplementary.cause.')

## Check inconsistent death dates between data files
temp<-death.new.14_16[duplicated(death.new.14_16$reference.key.), ]$reference.key.
temp2<-death.new.14_16 %>% filter(reference.key. %in% temp) %>% 
    distinct(reference.key., date.of.registered.death.)
temp<-temp2[duplicated(temp2$reference.key.), ]

death.new.14_16<-death.new.14_16 %>% distinct()

depression.14_16<-depression.14_16 %>% 
    inner_join(death.new.14_16, by = 'reference.key.') %>% 
    select(-date.of.registered.death..x, -exact.date.of.death..x, -death.cause.main.cause..x, -death.cause.supplementary.cause..x) %>% 
    rename(date.of.registered.death. = date.of.registered.death..y,
           exact.date.of.death. = exact.date.of.death..y,
           death.cause.main.cause. = death.cause.main.cause..y,
           death.cause.supplementary.cause. = death.cause.supplementary.cause..y) %>% 
    mutate(date.of.registered.death. = ymd(date.of.registered.death.))

## Prevalent patients who are aged 10 or above and valid 
### 2014 cohort = 15016
### 2015 cohort = 13466
### 2016 cohort = 12216
### 2014-2016 cohort = 40698
depression.14_16<-depression.14_16 %>% 
    mutate(follow.up.end.date = ymd(date.of.registered.death.)) %>% 
    replace_na(list(follow.up.end.date = as.Date('2020-12-31'))) %>% 
    mutate(follow.up.period = as.numeric(interval(reference.date., follow.up.end.date), unit = 'day')+1) %>% 
    filter(follow.up.period > 1) %>% 
    mutate(below.10 = if_else(reference.date.<= as.Date('2014-12-31') & date.of.birth.yyyy.mm.dd. > as.Date('2003-12-31') |
                                  reference.date.>= as.Date('2015-01-01') & reference.date.<= as.Date('2015-12-31') & date.of.birth.yyyy.mm.dd. > as.Date('2004-12-31') |
                                  reference.date.>= as.Date('2016-01-01') & reference.date.<= as.Date('2016-12-31') & date.of.birth.yyyy.mm.dd. > as.Date('2005-12-31'), 1, 0)) %>%
    filter(below.10 == 0) %>% 
    rename(first.dep.date = reference.date.) %>%  
    distinct(reference.key., .keep_all = T) %>% 
    select(-below.10)

## Incident patients ----
## 2014 = 8223
## 2015 = 8685 (valid incident pts = 8653)
## 2016 = 8354 (valid incident pts = 8314)
## 2014-2016 = 25190
inc.ref<-readRDS("inc.ref.rds")
inc.ref.2015<-readRDS("inc.ref.2015.rds")
inc.ref.2016<-readRDS("inc.ref.2016.rds")
depression.14_16.inc<-depression.14_16 %>% 
    filter(reference.key. %in% inc.ref | reference.key. %in% inc.ref.2015 | reference.key. %in% inc.ref.2016)

# 1.3 Identify TRD development from incidence to December 2020 ----
inc.third.date<-readRDS("inc.third.date.rds")
inc.third.date.2015<-readRDS("inc.third.date.2015.rds")
inc.third.date.2016<-readRDS("inc.third.date.2016.rds")
inc.third.date.14_16<-rbind(inc.third.date, inc.third.date.2015)
inc.third.date.14_16<-rbind(inc.third.date.14_16, inc.third.date.2016)

depression.14_16.inc<-depression.14_16.inc %>% 
    left_join(inc.third.date.14_16) %>% 
    rename(trd.date = start) %>% 
    replace_na(list(trd.date = as.Date('2022-12-31'))) %>% 
    mutate(trd.status = if_else(trd.date <= as.Date('2020-12-31'), 1, 0)) %>% 
    mutate(trd.date = na_if(trd.date, as.Date('2022-12-31')))

# 1.4 Identify age, sex, baseline comorbidities at incidence ----
# Age group
depression.14_16.inc <- depression.14_16.inc %>% 
    mutate(age = calc_age(date.of.birth.yyyy.mm.dd., first.dep.date),
           age.group = if_else(age <= 24, '10-24',
                               if_else(age>=25 & age<=40, '25-40',
                                       if_else(age>=41 & age<=65, '41-65',
                                               if_else(age>=66, '65+', '')))))

# Baseline comorbidities at incidence
comorbidities.code.table<-read_excel("comorbidities_icd9codes.xlsx")
comorbidities.code<-comorbidities.code.table %>% pull(icd.9.code)

# Combine all dx data from 1993 to 2020
directory <- list.files("../../Depression/converted_data_all_dx/", full.names = T)
directory <- directory[1:29]
dx.14 <- merge_files(directory)
dx.14<-dx.14 %>% mutate(date.of.birth.yyyy.mm.dd. = as.Date(date.of.birth.yyyy.mm.dd.),
                 reference.date. = as.Date(reference.date.))

directory <- list.files("../../Depression/2015_cohort/converted_data_all_dx/", full.names = T)
directory <- directory[1:29]
dx.15 <- merge_files(directory)
dx.15<-dx.15 %>% mutate(date.of.birth.yyyy.mm.dd. = as.Date(date.of.birth.yyyy.mm.dd.),
                 reference.date. = as.Date(reference.date.))

directory <- list.files("../../Depression/2016_cohort/converted_data_all_dx/", full.names = T)
directory <- directory[1:29]
dx.16 <- merge_files(directory)
dx.16<-dx.16 %>% mutate(date.of.birth.yyyy.mm.dd. = as.Date(date.of.birth.yyyy.mm.dd.),
                        reference.date. = as.Date(reference.date.))

dx.14_16<-rbind(dx.14, dx.15)
dx.14_16<-rbind(dx.14_16, dx.16)

dx.14_16.inc<-dx.14_16 %>% 
    mutate(reference.date. = ymd(reference.date.)) %>% 
    inner_join(depression.14_16.inc[c('reference.key.','first.dep.date','trd.date')]) %>%  # restrict data to incident patients
    distinct()

#saveRDS(dx.14_16.inc, 'dx.14_16.inc.rds')
    
dx.14_16.inc.com<-dx.14_16.inc %>% 
    filter(grepl(paste(comorbidities.code, collapse = "|"), all.diagnosis.code.icd9.)) %>% # restrict data to comorbidities of interest
    arrange(reference.key., reference.date.)

dx.14_16.inc.com<-fuzzy_join(dx.14_16.inc.com, comorbidities.code.table,
                       by = c('all.diagnosis.code.icd9.' = 'icd.9.code'),
                       match_fun = stringr::str_detect,
                       mode = "left") # merge dx data and code table by grepl

all.com.14_16<-dx.14_16.inc.com %>% 
    arrange(reference.key., reference.date.) %>% 
    distinct(reference.key., comorbidity.name, .keep_all = T) %>%  # keep only the first record of each comorbidity
    left_join(depression.14_16.inc[c('reference.key.','date.of.registered.death.','follow.up.end.date')], by = 'reference.key.')

baseline.somatic<-all.com.14_16 %>% 
    filter(reference.date. < first.dep.date) %>% # find baseline 
    filter(comorbidity.type == 'Somatic') %>% 
    distinct(reference.key., .keep_all = T) %>% # the first somatic disease per reference key
    select(reference.key., comorbidity.type)

baseline.mental<-all.com.14_16 %>% 
    filter(reference.date. < first.dep.date) %>% # find baseline 
    filter(comorbidity.type == 'Mental') %>% 
    distinct(reference.key., .keep_all = T) %>% # the first mental disease per reference key
    select(reference.key., comorbidity.type)

depression.14_16.inc<-depression.14_16.inc %>% 
    mutate(baseline.somatic = if_else(reference.key. %in% baseline.somatic$reference.key., 1, 0),
           baseline.mental = if_else(reference.key. %in% baseline.mental$reference.key., 1, 0),
           baseline.medical = if_else(baseline.somatic == 1 | baseline.mental == 1, 1, 0))

# 1.5 Identify death (external and natural) after incidence until Dec 2020 ----
## The classification of natural and external causes are TEMPORARILY at this stage!!!
depression.14_16.inc<-depression.14_16.inc %>% 
    replace_na(list(date.of.registered.death. = as.Date('2022-12-31'))) %>% 
    mutate(death = if_else(date.of.registered.death. <= as.Date('2020-12-31'), 1, 0),
           edeath = if_else(death == 1 & grepl('S|T|V|W|X|Y', death.cause.main.cause.) | 
                                death == 1 & grepl('S|T|V|W|X|Y', death.cause.supplementary.cause.), 1, 0),
           ndeath = if_else(death == 1 & edeath == 0, 1, 0)) %>% 
    mutate(date.of.registered.death. = na_if(date.of.registered.death., as.Date('2022-12-31')))

# 1.6 Identified recovery, i.e. at least 1-year absence of depression-related dx / AD prescription across settings until 2021 ----
## 1.6.1 Patients who still had depression-related OP dx (all ranks) in 2021 ----
## Combine OP dx data
### 2014 cohort
directory <- list.files("../../Depression/converted_data_op/", full.names = T)
directory<-directory[2:27]
op.data.14.1 <- merge_files(directory)
op.data.14.1$appointment.date.yyyy.mm.dd.<-ymd(op.data.14.1$appointment.date.yyyy.mm.dd.)

directory <- list.files("../../Depression/converted_data_op/", full.names = T)
directory<-directory[28:35]
op.data.14.2 <- merge_files(directory)
op.data.14.2$appointment.date.yyyy.mm.dd.<-ymd(op.data.14.2$appointment.date.yyyy.mm.dd.)
op.data.14<-rbind(op.data.14.1, op.data.14.2)

### 2015 cohort
directory <- list.files("../../Depression/2015_cohort/converted_data_op/", full.names = T)
op.data.15 <- merge_files(directory)

### 2016 cohort
directory <- list.files("../../Depression/2016_cohort/converted_data_op/", full.names = T)
op.data.16 <- merge_files(directory)

op.data.14_16<-rbind(op.data.14, op.data.15)
op.data.14_16<-rbind(op.data.14_16, op.data.16)

op.data.14_16<-op.data.14_16 %>% 
    filter(reference.key. %in% depression.14_16.inc$reference.key.) %>% 
    mutate(appointment.date.yyyy.mm.dd. = ymd(appointment.date.yyyy.mm.dd.)) %>% 
    mutate(depression.dx = if_else(grepl('296.2|296.3|300.4|311', 
                                         paste(most.recent.diagnosis.1.diagnosis.code.icd9.)), 1, 0)) %>% 
    distinct()

#saveRDS(op.data.14_16, 'op.data.14_16.rds')

unrecovered.op<-op.data.14_16 %>% 
    filter(appointment.date.yyyy.mm.dd. >= as.Date('2021-01-01')) %>% 
    filter(depression.dx == 1) %>% 
    distinct(reference.key.) %>% pull(reference.key.) ### 8282 still had depression OP dx in 2021

## Patients who no longer had depression-related OP dx (all ranks) in 2021
depression.14_16.inc<-depression.14_16.inc %>% 
    mutate(recovered.op = if_else(reference.key. %in% unrecovered.op, 0, 1)) # 16908 had no 2021 OP record

recovered.op<-depression.14_16.inc %>% 
    filter(recovered.op == 1) %>% pull(reference.key.) # count = 16908 who recovered based on only OP data

# Recovered patients' last depression-related OP dx date
recovery.op.date<-op.data.14_16 %>% 
    filter(reference.key. %in% recovered.op) %>% 
    filter(depression.dx == 1) %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    group_by(reference.key.) %>% slice(n()) %>% # among the recovered, find their last depression-related dx
    select(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    mutate(recovery.op.date = appointment.date.yyyy.mm.dd. + days(180)) %>% 
    rename(last.dep.op.date = appointment.date.yyyy.mm.dd.)
# assume the recovery date is 180 days after the last dx if no further dx

depression.14_16.inc<-depression.14_16.inc %>% left_join(recovery.op.date) 

## 1.6.2 Patients who still had depression-related AE dx (all ranks) in 2021 ----
## Load AE dx data
### 2014 cohort
directory <- list.files("../../Depression/converted_data_ae/", full.names = T)
directory <- directory[c(1:2)]
ae.data.14 <- merge_files(directory)
ae.data.14 <- ae.data.14[c(1:10, 13:41)]
ae.data.14<-ae.data.14 %>% 
    rename(district.of.residence.cluster.data.up.to.31.03.2015. = district.of.residence.cluster.,
           district.of.residence.district.data.up.to.31.03.2015. = district.of.residence.district.)

directory <- list.files("../../Depression/converted_data_ae/", full.names = T)
directory <- directory[c(3:9)]
ae.data.14.2 <- merge_files(directory)
ae.data.14.2 <- ae.data.14.2[c(1:10, 13:26, 43:57)]
ae.data.14.2<-ae.data.14.2[c(1,3,4,2,5:39)]
colnames(ae.data.14) == colnames(ae.data.14.2)
ae.data.14<-rbind(ae.data.14, ae.data.14.2)

ae.data.14<-ae.data.14[c(1:24)]

### 2015 cohort
directory <- list.files("../../Depression/2015_cohort/converted_data_ae/", full.names = T)
ae.data.15 <- merge_files(directory)
ae.data.15<-ae.data.15[c(1,3,4,2,5:9,10, 13:26)]

### 2016 cohort
directory <- list.files("../../Depression/2016_cohort/converted_data_ae/", full.names = T)
ae.data.16 <- merge_files(directory)
ae.data.16<-ae.data.16[c(1,3,4,2,5:9,10, 13:26)]

ae.data.14_16<-rbind(ae.data.14, ae.data.15)
ae.data.14_16<-rbind(ae.data.14_16, ae.data.16)

ae.data.14_16<-ae.data.14_16 %>% 
    filter(reference.key. %in% depression.14_16.inc$reference.key.) %>% 
    mutate(attendance.date.yyyy.mm.dd. = as.Date(attendance.date.yyyy.mm.dd.),
           depression.dx = if_else(grepl('296.2|296.3|300.4|311', 
                                         paste(a.e.to.ip.ward.principal.diagnosis., a.e.to.ip.ward.diagnosis.rank.2.,
                                               a.e.to.ip.ward.diagnosis.rank.3.,a.e.to.ip.ward.diagnosis.rank.4.,
                                               a.e.to.ip.ward.diagnosis.rank.5.,a.e.to.ip.ward.diagnosis.rank.6.,
                                               a.e.to.ip.ward.diagnosis.rank.7.,a.e.to.ip.ward.diagnosis.rank.8.,
                                               a.e.to.ip.ward.diagnosis.rank.9.,a.e.to.ip.ward.diagnosis.rank.10.,
                                               a.e.to.ip.ward.diagnosis.rank.11.,a.e.to.ip.ward.diagnosis.rank.12.,
                                               a.e.to.ip.ward.diagnosis.rank.13.,a.e.to.ip.ward.diagnosis.rank.14.,
                                               a.e.to.ip.ward.diagnosis.rank.15.)), 1, 0)) %>% 
    distinct()

#saveRDS(ae.data.14_16, 'ae.data.14_16.rds')

unrecovered.ae<-ae.data.14_16 %>% 
    filter(attendance.date.yyyy.mm.dd. >= as.Date('2021-01-01')) %>% 
    filter(depression.dx == 1) %>% 
    distinct(reference.key.) %>% pull(reference.key.) # count = 656 still had AE records

## Patients who no longer had depression-related AE dx (all ranks) in 2021
depression.14_16.inc<-depression.14_16.inc %>% 
    mutate(recovered.ae = if_else(reference.key. %in% unrecovered.ae, 0, 1))

recovered.ae<-depression.14_16.inc %>% 
    filter(recovered.ae == 1) %>% pull(reference.key.) # count = 24893 who recovered based on OP data

# Recovered patients' last depression-related AE dx date
recovery.ae.date<-ae.data.14_16 %>% 
    filter(reference.key. %in% recovered.ae) %>% 
    filter(depression.dx == 1) %>% 
    arrange(reference.key., attendance.date.yyyy.mm.dd.) %>% 
    group_by(reference.key.) %>% slice(n()) %>% # among the recovered, find their last depression-related dx
    select(reference.key., attendance.date.yyyy.mm.dd.) %>% 
    mutate(recovery.ae.date = attendance.date.yyyy.mm.dd. + days(180)) %>% 
    rename(last.dep.ae.date = attendance.date.yyyy.mm.dd.)
# assume the recovery date is 180 days after the last dx if no further dx

depression.14_16.inc<-depression.14_16.inc %>% left_join(recovery.ae.date) 

## 1.6.3 Patients who still had depression-related IP dx (all ranks) in 2021 ----
## Load IP dx data
### 2014 chort
directory <- list.files("../../Depression/converted_data_ip_ward/", full.names = T)
directory <- directory[c(1:12)]
ip.data.14 <- merge_files(directory)
directory <- list.files("../../Depression/converted_data_ip_ward/", full.names = T)
directory <- directory[c(13:14)]
ip.data.14.2 <- merge_files(directory)
ip.data.14<-rbind(ip.data.14, ip.data.14.2)

### 2015 cohort
directory <- list.files("../../Depression/2015_cohort/converted_data_ip_ward/", full.names = T)
ip.data.15 <- merge_files(directory)

### 2016 cohort
directory <- list.files("../../Depression/2016_cohort/converted_data_ip_ward/", full.names = T)
ip.data.16 <- merge_files(directory)

ip.data.14_16<-rbind(ip.data.14, ip.data.15)
ip.data.14_16<-rbind(ip.data.14_16, ip.data.16)

ip.data.14_16<-ip.data.14_16 %>% 
    filter(reference.key. %in% depression.14_16.inc$reference.key.) %>% 
    mutate(admission.date.yyyy.mm.dd. = as.Date(admission.date.yyyy.mm.dd.),
           discharge.date.yyyy.mm.dd. = as.Date(discharge.date.yyyy.mm.dd.)) %>% 
    mutate(depression.dx = if_else(grepl('296.2|296.3|300.4|311', 
                                         paste(principal.diagnosis.code., diagnosis.rank.2.,
                                               diagnosis.rank.3.,diagnosis.rank.4.,
                                               diagnosis.rank.5.,diagnosis.rank.6.,
                                               diagnosis.rank.7.,diagnosis.rank.8.,
                                               diagnosis.rank.9.,diagnosis.rank.10.,
                                               diagnosis.rank.11.,diagnosis.rank.12.,
                                               diagnosis.rank.13.,diagnosis.rank.14.,
                                               diagnosis.rank.15.)), 1, 0)) %>% 
    distinct()

#saveRDS(ip.data.14_16, 'ip.data.14_16.rds')

unrecovered.ip<-ip.data.14_16 %>% 
    filter(discharge.date.yyyy.mm.dd. >= as.Date('2021-01-01')) %>% 
    filter(depression.dx == 1) %>% 
    distinct(reference.key.) %>% pull(reference.key.) # count = 390 still had IP record in 2021

## Patients who no longer had depression-related IP dx (all ranks) in 2021
depression.14_16.inc<-depression.14_16.inc %>% 
    mutate(recovered.ip = if_else(reference.key. %in% unrecovered.ip, 0, 1))

recovered.ip<-depression.14_16.inc %>% 
    filter(recovered.ip == 1) %>% pull(reference.key.) # count = 24800 who recovered based on IP data

# Recovered patients' last depression-related IP dx date
recovery.ip.date<-ip.data.14_16 %>% 
    filter(reference.key. %in% recovered.ip) %>% 
    filter(depression.dx == 1) %>% 
    arrange(reference.key., discharge.date.yyyy.mm.dd.) %>% 
    group_by(reference.key.) %>% slice(n()) %>% # among the recovered, find their last depression-related dx
    select(reference.key., discharge.date.yyyy.mm.dd.) %>% 
    mutate(recovery.ip.date = discharge.date.yyyy.mm.dd. + days(180)) %>% 
    rename(last.dep.ip.date = discharge.date.yyyy.mm.dd.)
# assume the recovery date is 180 days after the last dx if no further dx

depression.14_16.inc<-depression.14_16.inc %>% left_join(recovery.ip.date) 

## 1.6.4 Patients who still had AD presciptions in 2021 ----
### 2014 cohort
#### Combine RX data between 2020 and 2021
rxdirectory <- list.files("../../Depression/converted_data_rx", full.names = T)
rxdirectory<-rxdirectory[119:141]
rx.14.2020.2021 <- merge_files(rxdirectory)
rx.14.2020.2021<-rx.14.2020.2021 %>% 
    mutate(dispensing.date.yyyy.mm.dd. = ymd(dispensing.date.yyyy.mm.dd.),
           prescription.start.date. = ymd(prescription.start.date.),
           prescription.end.date. = ymd(prescription.end.date.)) %>% 
    select(-drug.frequency....18, -drug.frequency....19, 
           -action.status....21, -therapeutic.classification.bnf.principal.,-dispensing.institution.) %>% 
    rename(drug.frequency. = drug.frequency....13, 
           action.status. = action.status....16)

rx.14.2020.2021<-rx.14.2020.2021[c(1:11,13,14,16,18,19,12,15,17,20:23)]

#### Combine RX data of 2014-2019 and 2020-2021
rxdirectory <- list.files("../../Depression/converted_data_rx", full.names = T)
rxdirectory<-rxdirectory[c(1:65, 97:109, 114)]
rx.14.2014.2019 <- merge_files(rxdirectory)
rx.14.2014.2019<-rx.14.2014.2019 %>% 
    mutate(dispensing.date.yyyy.mm.dd. = ymd(dispensing.date.yyyy.mm.dd.),
           prescription.start.date. = ymd(prescription.start.date.),
           prescription.end.date. = ymd(prescription.end.date.))

rx.14.ad<-rbind(rx.14.2014.2019, rx.14.2020.2021)

#### Confine RX data to Antidepressants
ad.list <- read_excel("../../Depression/ad_ap_list.xlsx", sheet = 1)
rx.14.ad <- rx.14.ad %>%
    filter(!(grepl("CETI01", drug.item.code.) & grepl("4.3.1", therapeutic.classification.bnf.))) %>%
    inner_join(ad.list, by = "drug.item.code.") 

### 2015 cohort
rxdirectory <- list.files("../../Depression/2015_cohort/converted_data_rx", full.names = T)
rx.15 <- merge_files(rxdirectory)
rx.15.ad <- rx.15 %>%
    filter(!(grepl("CETI01", drug.item.code.) & grepl("4.3.1", therapeutic.classification.bnf.))) %>%
    inner_join(ad.list, by = "drug.item.code.") 
rx.15.ad<-rx.15.ad[c(1:11,17,12,13,19,14,20,15,16,23,24,25,26,46)]
rx.15.ad<-rx.15.ad %>% 
    mutate(dispensing.date.yyyy.mm.dd. = ymd(dispensing.date.yyyy.mm.dd.),
           prescription.start.date. = ymd(prescription.start.date.),
           prescription.end.date. = ymd(prescription.end.date.)) %>% 
    rename(drug.frequency. = drug.frequency....13, 
           action.status. = action.status....16)

### 2016 cohort
rxdirectory <- list.files("../../Depression/2016_cohort/converted_data_rx", full.names = T)
rx.16 <- merge_files(rxdirectory)
rx.16.ad <- rx.16 %>%
    filter(!(grepl("CETI01", drug.item.code.) & grepl("4.3.1", therapeutic.classification.bnf.))) %>%
    inner_join(ad.list, by = "drug.item.code.") 
rx.16.ad<-rx.16.ad[c(1:11,17,12,13,18,14,20,15,16,23,24,25,26,40)]
rx.16.ad<-rx.16.ad %>% 
    mutate(dispensing.date.yyyy.mm.dd. = ymd(dispensing.date.yyyy.mm.dd.),
           prescription.start.date. = ymd(prescription.start.date.),
           prescription.end.date. = ymd(prescription.end.date.)) %>% 
    rename(drug.frequency. = drug.frequency....13, 
           action.status. = action.status....16)

rx.14_16.ad<-rbind(rx.14.ad, rx.15.ad)
rx.14_16.ad<-rbind(rx.14_16.ad, rx.16.ad)
rx.14_16.ad<-rx.14_16.ad %>% 
    filter(reference.key. %in% depression.14_16.inc$reference.key.) %>% 
    distinct()

#saveRDS(rx.14_16.ad, 'rx.14_16.ad.rds')

unrecovered.rx<-rx.14_16.ad %>% 
    filter(prescription.end.date. >= as.Date('2021-01-01')) %>% 
    distinct(reference.key.) %>% pull(reference.key.) # count = 12403 who still had AD rx in 2021

## Patients who no longer had AD rx in 2021
depression.14_16.inc<-depression.14_16.inc %>% 
    mutate(recovered.rx = if_else(reference.key. %in% unrecovered.rx, 0, 1))

recovered.rx<-depression.14_16.inc %>% 
    filter(recovered.rx == 1) %>% pull(reference.key.) # count = 12787 who recovered based on RX data

# Recovered patients' last prescription end date
recovery.rx.date<-rx.14_16.ad %>% 
    filter(reference.key. %in% recovered.rx) %>% 
    arrange(reference.key., prescription.end.date.) %>% 
    group_by(reference.key.) %>% slice(n()) %>% # among the recovered, find their last AD rx date
    select(reference.key., prescription.end.date.) %>% 
    rename(last.ad.rx.date = prescription.end.date.)

depression.14_16.inc<-depression.14_16.inc %>% left_join(recovery.rx.date) 

# Recovered patients using both RX and DX definitions
depression.14_16.inc<-depression.14_16.inc %>% 
    # Count recovery only if patients had no RX AD records and no DX records in 2021
    mutate(recovered = if_else(recovered.op == 0 | recovered.ip == 0 | recovered.ae == 0 | recovered.rx == 0, 0, 1)) %>% 
    # Change temporarily the recovery dates to unrealistic ones for NA
    replace_na(list(recovery.op.date = as.Date('2000-01-01'))) %>% 
    replace_na(list(recovery.ip.date = as.Date('2000-01-01'))) %>% 
    replace_na(list(recovery.ae.date = as.Date('2000-01-01'))) %>% 
    replace_na(list(last.ad.rx.date = as.Date('2000-01-01'))) %>% 
    # Choose recovery date: the later date of last AD RX data or last dx+180 date
    mutate(recovery.date = if_else(recovery.op.date > last.ad.rx.date, recovery.op.date, last.ad.rx.date),
           recovery.date = if_else(recovery.ip.date > recovery.date, recovery.ip.date, recovery.date),
           recovery.date = if_else(recovery.ae.date > recovery.date, recovery.ae.date, recovery.date)) %>% 
    mutate(first.dep.date180 = first.dep.date + days(180),
           recovery.date = if_else(last.ad.rx.date == as.Date('2000-01-01') & 
                                       recovery.op.date == as.Date('2000-01-01') & 
                                       recovery.ip.date == as.Date('2000-01-01') & 
                                       recovery.ae.date == as.Date('2000-01-01') & 
                                       recovered == 1, first.dep.date180, recovery.date)) %>% 
    # Change the unrealistic dates into 2022-12-31
    replace_na(list(date.of.registered.death. = as.Date('2022-12-31'))) %>% 
    mutate(recovery.date = na_if(recovery.date, as.Date('2000-01-01'))) %>% 
    replace_na(list(recovery.date = as.Date('2022-12-31'))) %>% 
    # For recovery dates on or later than 2021-01-01, treat them as unrecovered as of Dec 2020
    mutate(recovered = if_else(recovery.date >= as.Date('2021-01-01'), 0, recovered),
           recovery.date = if_else(recovery.date >= as.Date('2021-01-01'), as.Date('2022-12-31'), recovery.date)) %>% 
    # For recovery dates within 180 days after first depression dx, use 180 days after depression
    mutate(recovery.date = if_else(recovery.date < first.dep.date180,first.dep.date180, recovery.date)) %>% 
    # Death and recovery should be mutually exclusive, i.e. no death = 1 and recovery = 1
    # Two conditions of same recovery date and death date: 
    ## 1) death = 0 and recovery = 0
    ## 2) death = 1 and recovery = 1, which should be changed to death = 1 and recovery = 0
    mutate(death = if_else(recovery.date < date.of.registered.death., 0, death), # to make clear the death = 1 and recovered = 1
           recovered = if_else(date.of.registered.death. <= recovery.date, 0, recovered)) %>% # ditto, and correct those with same dates
    mutate(adj.death.date = if_else(death == 0, as.Date('2022-12-31'), date.of.registered.death.),
           adj.recovery.date = if_else(recovered == 0, as.Date('2022-12-31'), recovery.date)) %>% 
    # Change the unrealistic dates back to NAs
    mutate(date.of.registered.death. = na_if(date.of.registered.death., as.Date('2022-12-31')),
           recovery.date = na_if(recovery.date, as.Date('2022-12-31'))) %>% 
    # Adjust the info in causes of death
    mutate(ndeath = if_else(death == 0, 0, ndeath),
           edeath = if_else(death == 0, 0, edeath))

## Check
depression.14_16.inc %>% filter(trd.status == 1 & recovered == 1) %>% filter(trd.date >= adj.recovery.date) %>% count() # Count = 0
depression.14_16.inc %>% filter(trd.status == 1 & recovered == 1) %>% filter(trd.date >= adj.death.date) %>% count() # Count = 0
depression.14_16.inc %>% filter(trd.status == 1 & recovered == 1) %>% filter(first.dep.date >= adj.recovery.date) %>% count() # Count = 0
depression.14_16.inc %>% filter(trd.status == 1 & recovered == 1) %>% filter(first.dep.date >= adj.death.date) %>% count() # Count = 0
depression.14_16.inc %>% filter(trd.status == 1 & recovered == 1) %>% filter(first.dep.date180 > adj.recovery.date) %>% count() # Count = 0

# 1.7 Identify comorbidities acquisition after incident depression (for non-TRD) until Dec 2020 / TRD development / death / recovery ----
## Remember to revise if the diseases are after recovery!!
## Remember to check if anyone becomes TRD / acquire comorbidity on study end date!!
post.dep.medical<-all.com.14_16 %>% 
    replace_na(list(trd.date = as.Date('2020-12-31'))) %>% # No one turned TRD on 2020-12-31
    left_join(depression.14_16.inc[c('adj.death.date','adj.recovery.date','reference.key.')]) %>% 
    mutate(temp.censor.date = if_else(adj.death.date > adj.recovery.date, adj.recovery.date, adj.death.date),
           temp.censor.date = if_else(temp.censor.date > trd.date, trd.date, temp.censor.date)) %>% 
    filter(reference.date. >= first.dep.date & reference.date. < temp.censor.date) %>% # after depression but before TRD / study end / death / recovery
    distinct(reference.key., .keep_all = T) %>% # the first disease per reference key
    select(reference.key., comorbidity.type, reference.date.)

depression.14_16.inc<-depression.14_16.inc %>% 
    mutate(post.dep.medical = if_else(reference.key. %in% post.dep.medical$reference.key., 1, 0)) %>% 
    left_join(post.dep.medical, by ='reference.key.') %>% 
    replace_na(list(comorbidity.type = 'None')) %>% 
    rename(post.dep.medical.date = reference.date.,
           post.dep.medical.type = comorbidity.type)

# Check
depression.14_16.inc %>% filter(post.dep.medical.date < first.dep.date) %>% count() 
depression.14_16.inc %>% filter(post.dep.medical.date > trd.date) %>% count() 
depression.14_16.inc %>% filter(post.dep.medical.date > adj.death.date) %>% count() 
depression.14_16.inc %>% filter(post.dep.medical.date > adj.recovery.date) %>% count() 

#saveRDS(depression.14_16.inc, 'depression.14_16.inc.rds')

# 1.8 Identify comorbidities acquisition after TRD until Dec 2020  ----
## Remember to revise if the diseases are after recovery!!
## Remember to check if anyone becomes TRD / acquire comorbidity on study end date!!
post.trd.medical<-all.com.14_16 %>% 
    group_by(reference.key.) %>% 
    distinct(comorbidity.name, .keep_all = T) %>% 
    replace_na(list(trd.date = as.Date('2020-12-31'))) %>% 
    left_join(depression.14_16.inc[c('adj.death.date','adj.recovery.date','reference.key.')]) %>% 
    mutate(temp.censor.date = if_else(adj.death.date > adj.recovery.date, adj.recovery.date, adj.death.date)) %>% 
    filter(reference.date. >= trd.date & reference.date. <= as.Date('2020-12-31')) %>% # after TRD but before study end 
    filter(reference.date. >= trd.date & reference.date. <= temp.censor.date) %>% # after TRD but before recovery / death
    filter(trd.date != as.Date('2020-12-31')) %>% 
    distinct(reference.key., .keep_all = T) %>% # the first disease per reference key
    select(reference.key., comorbidity.type, reference.date.)

depression.14_16.inc<-depression.14_16.inc %>% 
    mutate(post.trd.medical = if_else(reference.key. %in% post.trd.medical$reference.key., 1, 0)) %>% 
    left_join(post.trd.medical, by ='reference.key.') %>% 
    replace_na(list(comorbidity.type = 'None')) %>% 
    rename(post.trd.medical.date = reference.date.,
           post.trd.medical.type = comorbidity.type)

# Check
depression.14_16.inc %>% filter(post.trd.medical.date < first.dep.date) %>% count() 
depression.14_16.inc %>% filter(post.trd.medical.date < trd.date) %>% count() 
depression.14_16.inc %>% filter(post.trd.medical.date > adj.death.date) %>% count() 
depression.14_16.inc %>% filter(post.trd.medical.date > adj.recovery.date) %>% count() 
depression.14_16.inc %>% filter(recovered == 1 & death == 1) %>% count() 

# 2.1 Regression for incidence to TRD ----
## Censoring = TRD (outcome) / study end / death / recovery / Pre-TRD comorbidity
temp<-depression.14_16.inc %>% 
    mutate(censor.date = trd.date) %>%  # time to TRD
    replace_na(list(censor.date = as.Date('2020-12-31'))) %>% # set to study end date for non-TRD
    replace_na(list(trd.date = as.Date('2022-12-31'))) %>% 
    mutate(temp.censor.date = if_else(adj.death.date > adj.recovery.date, adj.recovery.date, adj.death.date)) %>% 
    mutate(censor.date = if_else(censor.date > temp.censor.date, temp.censor.date, censor.date)) %>% # censor at death or recovery whichever earlier
    replace_na(list(post.dep.medical.date = as.Date('2022-12-31'))) %>% 
    mutate(censor.date = if_else(censor.date > post.dep.medical.date, post.dep.medical.date, censor.date)) %>% # censor at pre-TRD comorbidity
    mutate(time = calc_time(first.dep.date, censor.date),
           time = time/365.25,
           event = trd.status,
           event = if_else(trd.date > censor.date & trd.date != as.Date('2022-12-31'), 0, event))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose lognormal 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose lognormal 

fit<-survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')
summary(fit)

# Lognormal S(t) for all subgroups
# Age 10-24, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result<-x
    }else{
        result<-rbind(result,x)
    }
}
# Age 10-24, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 10-24, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 10-24, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)

#write.csv(result, file = 'st_ntrd2trd2.csv')
#write.csv(result, file = 'st_com2dead2.csv')
#write.csv(result, file = 'st_com2rec2.csv')
#write.csv(result, file = 'st_com2trd2.csv')
#write.csv(result, file = 'st_ntrd2rec2.csv')
#write.csv(result, file = 'st_com2trd2.csv')
#write.csv(result, file = 'st_tcom2rec2.csv')

# 2.2 Regression for NTRD to comorbidities ----
## Censoring =  pre-TRD comorbidity (outcome) / study end / death / recovery / TRD
temp<-depression.14_16.inc %>% 
    mutate(censor.date = post.dep.medical.date) %>%  # time to pre-TRD comorbidities
    replace_na(list(censor.date = as.Date('2020-12-31'))) %>% # set to study end date for no pre-TRD comorbidities
    replace_na(list(post.dep.medical.date = as.Date('2022-12-31'))) %>% 
    mutate(temp.censor.date = if_else(adj.death.date > adj.recovery.date, adj.recovery.date, adj.death.date)) %>% 
    mutate(censor.date = if_else(censor.date > temp.censor.date, temp.censor.date, censor.date)) %>% # censor at death or recovery whichever earlier
    replace_na(list(trd.date = as.Date('2022-12-31'))) %>% 
    mutate(censor.date = if_else(censor.date > trd.date, trd.date, censor.date)) %>% # censor at TRD
    mutate(time = calc_time(first.dep.date, censor.date),
           time = time/365.25,
           event = post.dep.medical,
           event = if_else(post.dep.medical.date > censor.date & post.dep.medical.date != as.Date('2022-12-31'), 0, event))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose weibull 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose weibull 

weibull<-WeibullReg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp)
weibull

fit<-(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) #27774.88
summary(fit)

# 2.3 Regression for NTRD to any deaths ----
## Censoring =  Death (outcome) / study end / recovery / TRD / Pre-TRD comorbidity
## One person had death date = study end date
temp<-depression.14_16.inc %>% 
    mutate(censor.date = adj.death.date) %>%  # time to death
    mutate(censor.date = if_else(censor.date == as.Date('2022-12-31'), as.Date('2020-12-31'), censor.date)) %>% # set to study end date for alive persons
    mutate(censor.date = if_else(censor.date > adj.recovery.date, adj.recovery.date, censor.date)) %>% # censor at recovery 
    replace_na(list(trd.date = as.Date('2022-12-31'))) %>% 
    mutate(censor.date = if_else(censor.date > trd.date, trd.date, censor.date)) %>% # censor at TRD
    replace_na(list(post.dep.medical.date = as.Date('2022-12-31'))) %>% 
    mutate(censor.date = if_else(censor.date > post.dep.medical.date, post.dep.medical.date, censor.date)) %>% # censor at pre-TRD comorbidity
    mutate(time = calc_time(first.dep.date, censor.date),
           time = time/365.25,
           event = death,
           event = if_else(adj.death.date > censor.date & adj.death.date != as.Date('2022-12-31'), 0, event))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
AIC(flexsurvreg(formula = Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = "gompertz")) 
## Choose loglogistic
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
BIC(flexsurvreg(formula = Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = "gompertz")) 
## Choose loglogistic

fit<-survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')
summary(fit)

# 2.4 Regression for NTRD to recovery ----
## Censoring =  Recovery (outcome) / study end / death / TRD / Pre-TRD comorbidity
## One person had recovery date = study end date
temp<-depression.14_16.inc %>% 
    mutate(censor.date = adj.recovery.date) %>%  # time to recovery
    mutate(censor.date = if_else(censor.date == as.Date('2022-12-31'), as.Date('2020-12-31'), censor.date)) %>% # set to study end date for unrecovered
    mutate(censor.date = if_else(censor.date > adj.death.date, adj.death.date, censor.date)) %>% # censor at death
    replace_na(list(trd.date = as.Date('2022-12-31'))) %>% 
    mutate(censor.date = if_else(censor.date > trd.date, trd.date, censor.date)) %>% # censor at TRD
    replace_na(list(post.dep.medical.date = as.Date('2022-12-31'))) %>% 
    mutate(censor.date = if_else(censor.date > post.dep.medical.date, post.dep.medical.date, censor.date)) %>% # censor at pre-TRD comorbidity
    mutate(time = calc_time(first.dep.date, censor.date),
           time = time/365.25,
           event = recovered,
           event = if_else(adj.recovery.date > censor.date & adj.recovery.date != as.Date('2022-12-31'), 0, event))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose lognormal
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose lognormal

fit<-survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')
summary(fit)

# Lognormal S(t) for all subgroups
# Age 10-24, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result<-x
    }else{
        result<-rbind(result,x)
    }
}
# Age 10-24, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 10-24, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 10-24, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)

#write.csv(result, file = 'st_ntrd2rec2.csv')

# 2.5 Regression for TRD to post-TRD comorbidities ----
## Censoring =  post-TRD comorbidity / study end / death / recovery
temp<-depression.14_16.inc %>% 
    filter(trd.status == 1) %>% 
    mutate(censor.date = post.trd.medical.date) %>%  # time to post-TRD comorbidities
    replace_na(list(censor.date = as.Date('2020-12-31'))) %>% # set to study end date 
    replace_na(list(post.trd.medical.date = as.Date('2022-12-31'))) %>% 
    mutate(temp.censor.date = if_else(adj.death.date > adj.recovery.date, adj.recovery.date, adj.death.date)) %>% 
    mutate(censor.date = if_else(censor.date > temp.censor.date, temp.censor.date, censor.date)) %>% # censor at death or recovery whichever earlier
    mutate(time = calc_time(trd.date, censor.date),
           time = time/365.25,
           event = post.trd.medical,
           event = if_else(post.trd.medical.date > censor.date & post.trd.medical.date != as.Date('2022-12-31'), 0, event))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose weibull
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose weibull

weibull<-WeibullReg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp)
weibull

# 2.6 Regression for TRD to any deaths ----
## Censoring =  Death (outcome) / study end / recovery / Post-TRD comorbidity
## One person had death date = study end date
temp<-depression.14_16.inc %>% 
    filter(trd.status == 1) %>% 
    mutate(censor.date = adj.death.date) %>%  # time to death
    mutate(censor.date = if_else(censor.date == as.Date('2022-12-31'), as.Date('2020-12-31'), censor.date)) %>% # set to study end date for alive persons
    mutate(censor.date = if_else(censor.date > adj.recovery.date, adj.recovery.date, censor.date)) %>% # censor at recovery 
    replace_na(list(post.trd.medical.date = as.Date('2022-12-31'))) %>% 
    mutate(censor.date = if_else(censor.date > post.trd.medical.date, post.trd.medical.date, censor.date)) %>% # censor at post-TRD comorbidity
    mutate(time = calc_time(trd.date, censor.date),
           time = time/365.25,
           event = death,
           event = if_else(adj.death.date > censor.date & adj.death.date != as.Date('2022-12-31'), 0, event)) %>% 
    mutate(above.65 = if_else(age >= 66, 1, 0),
           event = as.numeric(event),
           time = as.numeric(time))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
AIC(flexsurvreg(formula = Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = "gompertz")) 
## Choose loglogistic
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
BIC(flexsurvreg(formula = Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = "gompertz")) 
## Choose loglogistic

fit<-survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')
summary(fit)

# 2.7 Regression for TRD to recovery ----
## Censoring =  Recovery (outcome) / study end / death / Post-TRD comorbidity
## One person had recovery date = study end date
temp<-depression.14_16.inc %>% 
    filter(trd.status == 1) %>% 
    mutate(censor.date = adj.recovery.date) %>%  # time to recovery
    mutate(censor.date = if_else(censor.date == as.Date('2022-12-31'), as.Date('2020-12-31'), censor.date)) %>% # set to study end date for unrecovered
    mutate(censor.date = if_else(censor.date > adj.death.date, adj.death.date, censor.date)) %>% # censor at death
    replace_na(list(post.trd.medical.date = as.Date('2022-12-31'))) %>% 
    mutate(censor.date = if_else(censor.date > post.trd.medical.date, post.trd.medical.date, censor.date)) %>% # censor at post-TRD comorbidity
    mutate(time = calc_time(trd.date, censor.date),
           time = time/365.25,
           event = recovered,
           event = if_else(adj.recovery.date > censor.date & adj.recovery.date != as.Date('2022-12-31'), 0, event))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose weibull
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose weibull

weibull<-WeibullReg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp)
weibull

# 2.8 Regression for pre-TRD comorbidities to TRD ----
## Censoring =  TRD (outcome) / study end / recovery / death
temp<-depression.14_16.inc %>% 
    mutate(censor.date = trd.date) %>%  # time to TRD
    replace_na(list(censor.date = as.Date('2020-12-31'))) %>% # set to study end date for non-TRD
    replace_na(list(trd.date = as.Date('2022-12-31'))) %>% 
    mutate(temp.censor.date = if_else(adj.death.date > adj.recovery.date, adj.recovery.date, adj.death.date)) %>% 
    mutate(censor.date = if_else(censor.date > temp.censor.date, temp.censor.date, censor.date)) %>% # censor at death or recovery whichever earlier
    mutate(time = calc_time(post.dep.medical.date, censor.date),
           time = time/365.25,
           event = trd.status,
           event = if_else(trd.date > censor.date & trd.date != as.Date('2022-12-31'), 0, event))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose weibull
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose weibull

fit<-survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')
summary(fit)

# Lognormal S(t) for all subgroups
# Age 10-24, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result<-x
    }else{
        result<-rbind(result,x)
    }
}
# Age 10-24, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 10-24, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 10-24, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)

#write.csv(result, file = 'st_com2trd2.csv')

# 2.9 Regression for pre-TRD comorbidities to any deaths ----
## Censoring = death (outcome) / study end / recovery / TRD
## One person had death date = study end date
temp<-depression.14_16.inc %>% 
    filter(post.dep.medical == 1) %>% 
    mutate(censor.date = adj.death.date) %>%  # time to death
    mutate(censor.date = if_else(censor.date == as.Date('2022-12-31'), as.Date('2020-12-31'), censor.date)) %>% # set to study end date for alive persons
    mutate(temp.censor.date = adj.recovery.date) %>% 
    mutate(censor.date = if_else(censor.date > adj.recovery.date, adj.recovery.date, censor.date)) %>% # censor at recovery 
    replace_na(list(trd.date = as.Date('2022-12-31'))) %>% 
    mutate(censor.date = if_else(censor.date > trd.date, trd.date, censor.date)) %>% # censor at TRD
    mutate(time = calc_time(post.dep.medical.date, censor.date),
           time = time/365.25,
           event = death,
           event = if_else(adj.death.date > censor.date & adj.death.date != as.Date('2022-12-31'), 0, event))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
AIC(flexsurvreg(formula = Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = "gompertz")) 
## Choose lognormal
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
BIC(flexsurvreg(formula = Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = "gompertz")) 
## Choose lognormal

fit<-survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')
summary(fit)

# Lognormal S(t) for all subgroups
# Age 10-24, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result<-x
    }else{
        result<-rbind(result,x)
    }
}
# Age 10-24, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 10-24, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 10-24, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)

write.csv(result, file = 'st_com2dead3.csv')

# 2.10 Regression for pre-TRD comorbidity to recovery ----
## Censoring =  Recovery (outcome) / study end / death / TRD 
temp<-depression.14_16.inc %>% 
    filter(post.dep.medical == 1) %>% 
    mutate(censor.date = adj.recovery.date) %>%  # time to recovery
    mutate(censor.date = if_else(censor.date == as.Date('2022-12-31'), as.Date('2020-12-31'), censor.date)) %>% # set to study end date for unrecovered
    mutate(censor.date = if_else(censor.date > adj.death.date, adj.death.date, censor.date)) %>% # censor at death
    replace_na(list(trd.date = as.Date('2022-12-31'))) %>% 
    mutate(censor.date = if_else(censor.date > trd.date, trd.date, censor.date)) %>% # censor at TRD
    mutate(time = calc_time(first.dep.date, censor.date),
           time = time/365.25,
           event = recovered,
           event = if_else(adj.recovery.date > censor.date & adj.recovery.date != as.Date('2022-12-31'), 0, event))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) #11238.96
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) #11295.77
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) #11168.26
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) #11032.54
AIC(flexsurvreg(formula = Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = "gompertz")) #11144.96
## Choose loglogistic
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) #11285.52
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) #11335.68
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) #11214.82
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) #11079.1
BIC(flexsurvreg(formula = Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = "gompertz")) #11191.52
## Choose loglogistic

flexsurvreg(formula = Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = "llogis")
summary(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) #4584.708

# 2.11 Regression for post-TRD comorbidity to death ----
## Censoring = death (outcome) / study end / recovery
## One person had death date = study end date
temp<-depression.14_16.inc %>% 
    filter(post.trd.medical == 1) %>% 
    mutate(censor.date = adj.death.date) %>%  # time to death
    mutate(censor.date = if_else(censor.date == as.Date('2022-12-31'), as.Date('2020-12-31'), censor.date)) %>% # set to study end date for alive persons
    mutate(censor.date = if_else(censor.date > adj.recovery.date, adj.recovery.date, censor.date)) %>% # censor at recovery 
    mutate(time = calc_time(post.trd.medical.date, censor.date),
           time = time/365.25,
           event = death,
           event = if_else(adj.death.date > censor.date & adj.death.date != as.Date('2022-12-31'), 0, event))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) #11238.96
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) #11295.77
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) #11168.26
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) #11032.54
AIC(flexsurvreg(formula = Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = "gompertz")) #11144.96
## Choose weibull
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) #11285.52
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) #11335.68
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) #11214.82
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) #11079.1
BIC(flexsurvreg(formula = Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = "gompertz")) #11191.52
## Choose weibull

weibull<-WeibullReg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp)
weibull

# 2.12 Regression for post-TRD comorbidity to recovery ----
## Censoring =  Recovery (outcome) / study end / death
## One person had recovery date = study end date
temp<-depression.14_16.inc %>% 
    filter(post.trd.medical == 1) %>% 
    mutate(censor.date = adj.recovery.date) %>%  # time to recovery
    mutate(censor.date = if_else(censor.date == as.Date('2022-12-31'), as.Date('2020-12-31'), censor.date)) %>% # set to study end date for unrecovered
    mutate(censor.date = if_else(censor.date > adj.death.date, adj.death.date, censor.date)) %>% # censor at death
    mutate(time = calc_time(post.trd.medical.date, censor.date),
           time = time/365.25,
           event = recovered,
           event = if_else(adj.recovery.date > censor.date & adj.recovery.date != as.Date('2022-12-31'), 0, event))

AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
AIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose lognormal
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'w')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'exponential')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'loglogistic')) 
BIC(survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')) 
## Choose lognormal

fit<-survreg(Surv(time, event) ~ factor(age.group) + factor(sex.) + factor(baseline.medical), data = temp, dist = 'lognormal')
summary(fit)

# Lognormal S(t) for all subgroups
# Age 10-24, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result<-x
    }else{
        result<-rbind(result,x)
    }
}
# Age 10-24, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 10-24, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 10-24, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 25-40, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[2]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 41-65, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[3]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, F, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, F, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, M, history+ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[5]+fit$coefficients[6], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)
# Age 65+, M, history-ve
for (i in 1:10){
    x<-1-plnorm(i, meanlog = fit$coefficients[1]+fit$coefficients[4]+fit$coefficients[5], sdlog = fit$scale)
    if(i==1){
        result2<-x
    }else{
        result2<-rbind(result2,x)
    }
}
result<-cbind(result, result2)

#write.csv(result, file = 'st_tcom2rec2.csv')

# 3.1 Model validation (Numerical, at 4-year horizon) ----
## Event = death
## Censoring =  Death (outcome) / study end / recovery
## One person had death date = study end date
temp<-depression.14_16.inc %>% 
    mutate(censor.date = adj.death.date) %>%  # time to death
    mutate(censor.date = if_else(censor.date == as.Date('2022-12-31'), as.Date('2020-12-31'), censor.date)) %>% # set to study end date for alive persons
    mutate(temp.censor.date = if_else(adj.death.date > adj.recovery.date, adj.recovery.date, adj.death.date)) %>% 
    mutate(censor.date = if_else(censor.date > temp.censor.date, temp.censor.date, censor.date)) %>% # censor at death or recovery whichever earlier
    mutate(time = calc_time(first.dep.date, censor.date),
           #time = time/365.25,
           event = death,
           event = if_else(adj.death.date > censor.date & adj.death.date != as.Date('2022-12-31'), 0, event))

# Observed deaths after each cycle
for (i in 1:10){
    surv.group.kmfit<-survfit(Surv(time, event) ~ 1, data = temp)
    x<-sum(summary(surv.group.kmfit, times = c(1:365.25*i))$n.event) 
    if(i==1){
        result<-x
    }else{
        result<-rbind(result,x)
    }
}
result

# Observed deaths by age / sex / baseline co-morbidity after cycle 4
surv.group.kmfit<-survfit(Surv(time, event) ~ factor(age.group), data = temp)
x<-summary(surv.group.kmfit, times = c(1:365.25*4))$n.event
sum(x[1:(length(x)/2)])
sum(x[(length(x)/2+1):length(x)])

surv.group.kmfit<-survfit(Surv(time, event) ~ factor(sex.), data = temp)
sum(x[1:(length(x)/2)])
sum(x[(length(x)/2+1):length(x)])

surv.group.kmfit<-survfit(Surv(time, event) ~ factor(baseline.medical), data = temp)
sum(x[1:365])
sum(x[366:730])
sum(x[731:1095])
sum(x[1096:1460])

## Event = recovery
## Censoring =  recovery (outcome) / study end / death
## One person had death date = study end date
temp<-depression.14_16.inc %>% 
    mutate(censor.date = adj.recovery.date) %>%  # time to recovery
    mutate(censor.date = if_else(censor.date == as.Date('2022-12-31'), as.Date('2020-12-31'), censor.date)) %>% # set to study end date for unrecovered
    mutate(censor.date = if_else(censor.date > adj.death.date, adj.death.date, censor.date)) %>% # censor at death
    mutate(time = calc_time(first.dep.date, censor.date),
           #time = time/365.25,
           event = recovered,
           event = if_else(adj.recovery.date > censor.date & adj.recovery.date != as.Date('2022-12-31'), 0, event))

# Observed deaths after each cycle
for (i in 1:10){
    surv.group.kmfit<-survfit(Surv(time, event) ~ 1, data = temp)
    x<-sum(summary(surv.group.kmfit, times = c(1:365.25*i))$n.event) 
    if(i==1){
        result<-x
    }else{
        result<-rbind(result,x)
    }
}
result

# Observed deaths by age / sex / baseline co-morbidity after cycle 4
surv.group.kmfit<-survfit(Surv(time, event) ~ factor(age.group), data = temp)
x<-summary(surv.group.kmfit, times = c(1:365.25*4))$n.event
sum(x[1:(length(x)/2)])
sum(x[(length(x)/2+1):length(x)])

surv.group.kmfit<-survfit(Surv(time, event) ~ factor(sex.), data = temp)
sum(x[1:(length(x)/2)])
sum(x[(length(x)/2+1):length(x)])

surv.group.kmfit<-survfit(Surv(time, event) ~ factor(baseline.medical), data = temp)
sum(x[1:365])
sum(x[366:730])
sum(x[731:1095])
sum(x[1096:1460])

# 3.2 Model validation (Visual inspection) ----
## Plot by sex
s <- with(temp, Surv(time,event))
fKM <- survfit(s ~ sex.,data=temp)
sWei <- survreg(s ~ as.factor(sex.),dist='weibull',data=temp)

pred.sex1 = predict(sWei, newdata=list(sex.='F'),type="quantile",p=seq(.01,.99,by=.01))
pred.sex2 = predict(sWei, newdata=list(sex.='M'),type="quantile",p=seq(.01,.99,by=.01))

df = data.frame(y=seq(.99,.01,by=-.01), sex1=pred.sex1, sex2=pred.sex2)
df_long = gather(df, key= "sex", value="time", -y)

p = ggsurvplot(fKM, data = temp, risk.table = F)
p$plot = p$plot + geom_line(data=df_long, aes(x=time, y=y, group=sex))
p

## Plot by age
s <- with(temp, Surv(time,event))
fKM <- survfit(s ~ age.group,data=temp)
sWei <- survreg(s ~ as.factor(age.group),dist='lognormal',data=temp)

pred.sex1 = predict(sWei, newdata=list(age.group='10-24'),type="quantile",p=seq(.01,.99,by=.01))
pred.sex2 = predict(sWei, newdata=list(age.group='25-40'),type="quantile",p=seq(.01,.99,by=.01))
pred.sex3 = predict(sWei, newdata=list(age.group='41-65'),type="quantile",p=seq(.01,.99,by=.01))
pred.sex4 = predict(sWei, newdata=list(age.group='65+'),type="quantile",p=seq(.01,.99,by=.01))

df = data.frame(y=seq(.99,.01,by=-.01), sex1=pred.sex1, sex2=pred.sex2, sex3=pred.sex3, sex4=pred.sex4)
df_long = gather(df, key= "sex", value="time", -y)

p = ggsurvplot(fKM, data = temp, risk.table = F)
p$plot = p$plot + geom_line(data=df_long, aes(x=time, y=y, group=sex))
p

## Plot by baseline medical
s <- with(temp, Surv(time,event))
fKM <- survfit(s ~ baseline.medical,data=temp)
sWei <- survreg(s ~ as.factor(baseline.medical),dist='loglogistic',data=temp)

pred.sex1 = predict(sWei, newdata=list(baseline.medical='1'),type="quantile",p=seq(.01,.99,by=.01))
pred.sex2 = predict(sWei, newdata=list(baseline.medical='0'),type="quantile",p=seq(.01,.99,by=.01))

df = data.frame(y=seq(.99,.01,by=-.01), sex1=pred.sex1, sex2=pred.sex2)
df_long = gather(df, key= "sex", value="time", -y)
p = ggsurvplot(fKM, data = temp, risk.table = F)
p$plot = p$plot + geom_line(data=df_long, aes(x=time, y=y, group=sex))
p
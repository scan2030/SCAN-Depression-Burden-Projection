library(dplyr)
library(tidyverse)
library(lubridate)
library(readxl)
library(fuzzyjoin)
library(MatchIt)
library(MASS)
library(ggplot)
library(ggpubr)
library(fastDummies)
library(tableone)

# 0.1 Steps ----
#1 - Load files/data frames/merged datasets saved from R project "SCAN".
#2.1 - Perform PS matching 1:4 based on age, sex and baseline co-morbidity to get index date. Match 3 years separately.
#2.2 - Exclude patients who died before index date and check Table 1.
#3 - Get disease history after matching but before the index date: co-morbidity acquisition.
#4 - Get disease history after the index date: new-onset co-morbidity.
#5 - Identify inpatient, A&E, outpatient cost after index date.
#6 - Perform multiple negative binomial regression, adjust for: post-matching chx, TRD, new-onset co-morbidity.

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

# 1. Load files/data frames/merged datasets saved from R project "SCAN" ----
depression.14_16.inc<-readRDS('depression.14_16.inc.RDS') # Main dataframe containing all 2014-2016 incident patients
inc.ref<-readRDS('inc.ref.RDS') # 2014 incident patients (8223)
inc.ref.2015<-readRDS('inc.ref.2015.RDS') # 2015 incident patients 8685 (valid incident pts = 8653)
inc.ref.2016<-readRDS('inc.ref.2016.RDS') # 2016 incident patients 8354 (valid incident pts = 8314)
dx.14_16.inc<-readRDS('dx.14_16.inc.RDS') # All dx from 1993 to 2020 of all 2014-2016 incident patients
dx.14_16.inc.com<-readRDS('dx.14_16.inc.com.RDS') # Interested comorbidity dx from 1993 to 2020 of all 2014-2016 incident patients
all.com.14_16<-readRDS('all.com.14_16.RDS') # Interested comorbidity (1st record) dx from 1993 to 2020 of all 2014-2016 incident patients
ae.data.14_16<-readRDS('ae.data.14_16.RDS') # A&E dx from 1993 to 2021 of all 2014-2016 incident patients
ip.data.14_16<-readRDS('ip.data.14_16.RDS') # IP dx from 1993 to 2021 of all 2014-2016 incident patients

# 2. Perform PS matching 1:4 based on age, sex and baseline co-morbidity to get index date ----
## Match 3 years separately
### 2.1 PS matching for 2014 cohort ----
ca.14<-depression.14_16.inc %>% filter(reference.key. %in% inc.ref) # dataframe for cost analysis 2014
ps.model <- glm(trd.status ~ age + sex. + baseline.medical, data = ca.14, family=binomial)
ca.14$ps <- predict(ps.model, type='response')
head(ca.14)
m.out <- matchit(trd.status ~ ps, data = ca.14, method = "nearest", ratio = 4, caliper = 0.25)
summary(m.out)
ca.14<-match.data(m.out)

temp<-ca.14 %>% 
    filter(trd.status == 1) %>%
    dplyr::select(subclass, trd.date)
ca.14<-ca.14 %>% 
    left_join(temp, by = 'subclass') %>% 
  dplyr::select(-trd.date.x) %>% 
    rename(trd.date = trd.date.y) %>% 
    mutate(trd.date = if_else(trd.date < first.dep.date, first.dep.date, trd.date))

### 2.2 PS matching for 2015 cohort ----
ca.15<-depression.14_16.inc %>% filter(reference.key. %in% inc.ref.2015) # dataframe for cost analysis 2015
ps.model <- glm(trd.status ~ age + sex. + baseline.medical, data = ca.15, family=binomial)
ca.15$ps <- predict(ps.model, type='response')
head(ca.15)
m.out <- matchit(trd.status ~ ps, data = ca.15, method = "nearest", ratio = 4, caliper = 0.25)
summary(m.out)
ca.15<-match.data(m.out)

temp<-ca.15 %>% 
    filter(trd.status == 1) %>%
  dplyr::select(subclass, trd.date)
ca.15<-ca.15 %>% 
    left_join(temp, by = 'subclass') %>% 
  dplyr::select(-trd.date.x) %>% 
    rename(trd.date = trd.date.y) %>% 
    mutate(trd.date = if_else(trd.date < first.dep.date, first.dep.date, trd.date))

### 2.3 PS matching for 2016 cohort ----
ca.16<-depression.14_16.inc %>% filter(reference.key. %in% inc.ref.2016) # dataframe for cost analysis 2016
ps.model <- glm(trd.status ~ age + sex. + baseline.medical, data = ca.16, family=binomial)
ca.16$ps <- predict(ps.model, type='response')
head(ca.16)
m.out <- matchit(trd.status ~ ps, data = ca.16, method = "nearest", ratio = 4, caliper = 0.25)
summary(m.out)
ca.16<-match.data(m.out)

temp<-ca.16 %>% 
    filter(trd.status == 1) %>%
  dplyr::select(subclass, trd.date)
ca.16<-ca.16 %>% 
    left_join(temp, by = 'subclass') %>% 
  dplyr::select(-trd.date.x) %>% 
    rename(trd.date = trd.date.y) %>% 
    mutate(trd.date = if_else(trd.date < first.dep.date, first.dep.date, trd.date))

# 2.4 Combine matched data of 2014, 2015 and 2016 cohorts ----
ca.14_16 <- rbind(ca.14, ca.15)
ca.14_16 <- rbind(ca.14_16, ca.16)
ca.14_16 <- ca.14_16 %>% 
    rename(index.date = trd.date) %>%  # rename to avoid confusion
    mutate(follow.up.period = calc_time(index.date, follow.up.end.date)) # correct the FU period after matching
### NOTE: Some patients will have follow-up period <= 0. Normal. They will be excluded in main analysis.
#saveRDS(ca.14_16, 'ca.14_16.rds')
#ca.14_16<-readRDS('ca.14_16.rds')

# 3. Identify post-matching co-morbidity acquisition ----
## Definition = those who only became comorbid after depression but before index date (1/0)
all.com.14_16<-all.com.14_16 %>% 
    left_join(ca.14_16[c('reference.key.','index.date')], by ='reference.key.') %>% 
    replace_na(list(index.date = as.Date('2022-12-31'))) # hypothetical

post.match.somatic<-all.com.14_16 %>% 
    filter(reference.date. < index.date) %>% # find all wanted dx before index date  
    filter(comorbidity.type == 'Somatic') %>% 
    distinct(reference.key., .keep_all = T) %>% # the first somatic disease per reference key
  dplyr::select(reference.key., comorbidity.type)

post.match.mental<-all.com.14_16 %>% 
    filter(reference.date. < index.date) %>% # find all wanted dx before index date  
    filter(comorbidity.type == 'Mental') %>% 
    distinct(reference.key., .keep_all = T) %>% # the first mental disease per reference key
  dplyr::select(reference.key., comorbidity.type)

ca.14_16<-ca.14_16 %>% 
    mutate(post.match.somatic = if_else(reference.key. %in% post.match.somatic$reference.key., 1, 0),
           post.match.mental = if_else(reference.key. %in% post.match.mental$reference.key., 1, 0),
           post.match.medical = if_else(post.match.somatic == 1 | post.match.mental == 1, 1, 0))

# check
ca.14_16 %>% filter(baseline.medical == 1) %>% dplyr::select(baseline.medical, post.match.medical) %>% 
    summarise(b.sum = sum(baseline.medical), p.sum = sum(post.match.medical)) 
# all baseline comorbid pts should be found comorbid during post-match too

ca.14_16<-ca.14_16 %>% 
    mutate(post.match.comorbid = if_else(post.match.medical > baseline.medical, 1, 0)) # used for regression adjustment

# 4. Identify new-onset co-morbidity after index date ----
## Definition = acquire a new type of co-morbidity (increase in no. of comorbidity) compare with period before index date (1/0)
## Period = index date to study end / death
post.index.medical<-all.com.14_16 %>% 
    filter(reference.date. >= index.date & reference.date. <= follow.up.end.date) %>% # after depression but before follow-up end date
    distinct(reference.key., .keep_all = T) %>% # the first disease per reference key
  dplyr::select(reference.key., comorbidity.type, reference.date.)

ca.14_16<-ca.14_16 %>% 
    mutate(post.index.medical = if_else(reference.key. %in% post.index.medical$reference.key., 1, 0))

# 5. Identify inpatient, A&E, outpatient cost after index date ----
# 5.1 A&E utilization ----
ae.cost.14_16<-ae.data.14_16 %>% 
    inner_join(ca.14_16[c('reference.key.','index.date','follow.up.end.date','trd.status')], by = 'reference.key.') %>% 
    filter(attendance.date.yyyy.mm.dd. >= index.date & attendance.date.yyyy.mm.dd. <= follow.up.end.date) %>% 
    arrange(reference.key., attendance.date.yyyy.mm.dd., discharge.date.yyyy.mm.dd.) %>% 
    distinct() %>% 
    group_by(reference.key.) %>% 
    summarise(ae.episodes = n(), 
              ae.cost = ae.episodes*1230)

ca.14_16<-ca.14_16 %>% 
    left_join(ae.cost.14_16[c('reference.key.','ae.episodes','ae.cost')], by = 'reference.key.') %>% 
    replace_na(list(ae.episodes = 0)) %>% 
    replace_na(list(ae.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.ae.episodes.py = mean(ae.episodes/(follow.up.period/365.25)),
              s.ae.episodes.py = sd(ae.episodes/(follow.up.period/365.25)))

# 5.2 Inpatient utilization ----
# Sum up the LOS and cost by ward type
ip.cost.14_16<-ip.data.14_16 %>% 
    inner_join(ca.14_16[c('reference.key.','index.date','follow.up.end.date','trd.status')], by = 'reference.key.') %>% 
    filter(admission.date.yyyy.mm.dd. <= follow.up.end.date) %>% # restrict to data up to 2020
    mutate(los.of.ward.care.type.acute.general. = as.numeric(los.of.ward.care.type.acute.general.),
           los.of.ward.care.type.acute.general.acute. = as.numeric(los.of.ward.care.type.acute.general.acute.),
           los.of.ward.care.type.convalescent.rehabilitation.infirmary. = as.numeric(los.of.ward.care.type.convalescent.rehabilitation.infirmary.),
           los.of.ward.care.type.acute.general.high.dependency. = as.numeric(los.of.ward.care.type.acute.general.high.dependency.),
           los.of.ward.care.type.acute.general.intensive.care. = as.numeric(los.of.ward.care.type.acute.general.intensive.care.),
           los.of.ward.care.type.psychiatry.mentally.handicapped. = as.numeric(los.of.ward.care.type.psychiatry.mentally.handicapped.)) %>% 
    replace_na(list(los.of.ward.care.type.acute.general.acute. = 0)) %>% 
    replace_na(list(los.of.ward.care.type.acute.general. = 0)) %>% 
    replace_na(list(los.of.ward.care.type.acute.general.high.dependency. = 0)) %>% 
    replace_na(list(los.of.ward.care.type.acute.general.intensive.care. = 0)) %>% 
    replace_na(list(los.of.ward.care.type.convalescent.rehabilitation.infirmary. = 0)) %>% 
    replace_na(list(los.of.ward.care.type.psychiatry.mentally.handicapped. = 0)) %>% 
    filter(!discharge.date.yyyy.mm.dd. < index.date) %>% # don't count those who were discharged before index date
    arrange(reference.key., admission.date.yyyy.mm.dd., discharge.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., admission.date.yyyy.mm.dd., discharge.date.yyyy.mm.dd., .keep_all = T) %>%
    rename(los.acute = los.of.ward.care.type.acute.general.acute.,
           los.convalescent = los.of.ward.care.type.convalescent.rehabilitation.infirmary.,
           los.high = los.of.ward.care.type.acute.general.high.dependency.,
           los.icu = los.of.ward.care.type.acute.general.intensive.care.,
           los.psy = los.of.ward.care.type.psychiatry.mentally.handicapped.) %>% 
    mutate(los.general = los.acute + los.convalescent, 
           los.total = los.general + los.high + los.icu + los.psy) %>% 
    mutate(lead.days = calc_time(admission.date.yyyy.mm.dd., index.date), # hospital days before index date
           lead.days = if_else(lead.days <= 0, 0, lead.days),
           lag.days = calc_time(follow.up.end.date, discharge.date.yyyy.mm.dd.),
           lag.days = if_else(lag.days <= 0, 0, lag.days), # hospital days after follow-up end date
           adm.to.disc = calc_time(admission.date.yyyy.mm.dd., discharge.date.yyyy.mm.dd.), # calculated period from admission to discharge
           adj.los.general = los.general*(1 - (lead.days+lag.days)/adm.to.disc), # remove bed days out of follow-up period 
           adj.los.psy = los.psy*(1 - (lead.days+lag.days)/adm.to.disc),
           adj.los.high = los.high*(1 - (lead.days+lag.days)/adm.to.disc),
           adj.los.icu = los.icu*(1 - (lead.days+lag.days)/adm.to.disc),
           adj.los.total = adj.los.general+adj.los.psy+adj.los.high+adj.los.icu) %>% 
    group_by(reference.key.) %>% 
    summarise(los.general = sum(adj.los.general),
              los.high = sum(adj.los.high),
              los.icu = sum(adj.los.icu),
              los.psy = sum(adj.los.psy),
              los.total = sum(adj.los.total)) %>% 
    mutate(ip.general.cost = los.general*5100,
           ip.psy.cost = los.psy*2340,
           ip.high.cost = los.high*13650,
           ip.icu.cost = los.icu*24400,
           ip.cost = ip.general.cost + ip.psy.cost + ip.high.cost + ip.icu.cost)

ca.14_16<-ca.14_16 %>% 
    left_join(ip.cost.14_16[c('los.general','los.psy','los.high','los.icu','los.total',
                              'ip.general.cost','ip.psy.cost','ip.high.cost','ip.icu.cost','ip.cost','reference.key.')], by = 'reference.key.') %>% 
    replace_na(list(los.general = 0)) %>% 
    replace_na(list(los.psy = 0)) %>% 
    replace_na(list(los.high = 0)) %>% 
    replace_na(list(los.icu = 0)) %>% 
    replace_na(list(los.total = 0)) %>% 
    replace_na(list(ip.general.cost = 0)) %>% 
    replace_na(list(ip.psy.cost = 0)) %>% 
    replace_na(list(ip.high.cost = 0)) %>% 
    replace_na(list(ip.icu.cost = 0)) %>% 
    replace_na(list(ip.cost = 0))

ca.14_16 %>% 
    filter(follow.up.period > 0) %>% 
    mutate(follow.up.period = calc_time(index.date, follow.up.end.date)/365.25) %>% 
    group_by(trd.status) %>% 
    summarise(los.general.py = mean(los.general/follow.up.period),
              los.psy.py = mean(los.psy/follow.up.period),
              los.high.py = mean(los.high/follow.up.period),
              los.icu.py = mean(los.icu/follow.up.period),
              los.high.icu.py = mean((los.high+los.icu)/follow.up.period),
              los.total.py = mean(los.total/follow.up.period))

# 5.3 Outpatient utilization ----
## Combine OP dx data
### 2014 cohort
directory <- list.files("../../Depression/converted_data_op_serv/", full.names = T)
directory<-directory[1:72]
op.serv.data.14.1 <- merge_files(directory)
op.serv.data.14.1$appointment.date.yyyy.mm.dd.<-ymd(op.serv.data.14.1$appointment.date.yyyy.mm.dd.)

directory <- list.files("../../Depression/converted_data_op_serv/", full.names = T)
directory<-directory[73:84]
op.serv.data.14.2 <- merge_files(directory)
op.serv.data.14.2$appointment.date.yyyy.mm.dd.<-ymd(op.serv.data.14.2$appointment.date.yyyy.mm.dd.)
op.serv.data.14<-rbind(op.serv.data.14.1, op.serv.data.14.2)

### 2015 cohort
directory <- list.files("../../Depression/2015_cohort/converted_data_op_serv/", full.names = T)
op.serv.data.15 <- merge_files(directory)

### 2016 cohort
directory <- list.files("../../Depression/2016_cohort/converted_data_op_serv/", full.names = T)
op.serv.data.16 <- merge_files(directory)

op.serv.data.14_16<-rbind(op.serv.data.14, op.serv.data.15)
op.serv.data.14_16<-rbind(op.serv.data.14_16, op.serv.data.16)
op.serv.data.14_16<-op.serv.data.14_16 %>% 
    filter(reference.key. %in% depression.14_16.inc$reference.key.) %>% 
    mutate(appointment.date.yyyy.mm.dd. = ymd(appointment.date.yyyy.mm.dd.)) %>% 
    distinct()

#saveRDS(op.serv.data.14_16, 'op.serv.data.14_16.rds')
op.serv.data.14_16<-readRDS('op.serv.data.14_16.rds')

op.serv.data.14_16<-op.serv.data.14_16 %>% 
    inner_join(ca.14_16[c('reference.key.','index.date','follow.up.end.date')], by ='reference.key.') %>% 
    filter(appointment.date.yyyy.mm.dd. >= index.date) %>% 
    filter(appointment.date.yyyy.mm.dd. <= follow.up.end.date) 

# 5.3.1 OP Psychiatric day hospital (PDH) --------
pdh.cost.14_16<-op.serv.data.14_16 %>% 
    filter(grepl('^AHPD$|^PDH$', service.type.code.eis.)) %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., appointment.date.yyyy.mm.dd., .keep_all = TRUE) %>% 
    group_by(reference.key.) %>% 
    summarise(pdh.episodes = n(),
              pdh.cost = pdh.episodes*1260)

ca.14_16<-ca.14_16 %>% 
    left_join(pdh.cost.14_16[c('reference.key.','pdh.episodes','pdh.cost')], by = 'reference.key.') %>% 
    replace_na(list(pdh.episodes = 0)) %>% 
    replace_na(list(pdh.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(pdh.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(pdh.cost/(follow.up.period/365.25)))

# 5.3.2 OP Geriatric day hospital (GDH) --------
gdh.cost.14_16<-op.serv.data.14_16 %>% 
    filter(grepl('^AHGD$|^GDH', service.type.code.eis.)) %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., appointment.date.yyyy.mm.dd., .keep_all = TRUE) %>% 
    group_by(reference.key.) %>% 
    summarise(gdh.episodes = n(),
              gdh.cost = gdh.episodes*1960)

ca.14_16<-ca.14_16 %>% 
    left_join(gdh.cost.14_16[c('reference.key.','gdh.episodes','gdh.cost')], by = 'reference.key.') %>% 
    replace_na(list(gdh.episodes = 0)) %>% 
    replace_na(list(gdh.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(gdh.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(gdh.cost/(follow.up.period/365.25)))

# 5.3.3 OP Rehab day hospital (RDH) --------
rdh.cost.14_16<-op.serv.data.14_16 %>% 
    filter(grepl('^AHRD$|^RDP$|^RDPC$|^RDPN$|^RDPP$', service.type.code.eis.)) %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., appointment.date.yyyy.mm.dd., .keep_all = TRUE) %>% 
    group_by(reference.key.) %>% 
    summarise(rdh.episodes = n(),
              rdh.cost = rdh.episodes*1320)

ca.14_16<-ca.14_16 %>% 
    left_join(rdh.cost.14_16[c('reference.key.','rdh.episodes','rdh.cost')], by = 'reference.key.') %>% 
    replace_na(list(rdh.episodes = 0)) %>% 
    replace_na(list(rdh.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(rdh.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(rdh.cost/(follow.up.period/365.25)))

# 5.3.4 OP General Outpatient Clinic (GOPC) --------
gopc.cost.14_16<-op.serv.data.14_16 %>% 
    filter(grepl('^GOP$|^GOP IPC$|^NAHC|^RAMP|^Smoking', service.group.eis.) | grepl('^GOAD$', specialty.code.opas.)) %>% 
    filter(!grepl('Dressing|Injection', service.group.eis.)) %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., appointment.date.yyyy.mm.dd., .keep_all = TRUE) %>% 
    group_by(reference.key.) %>% 
    summarise(gopc.episodes = n(),
              gopc.cost = gopc.episodes*445)

ca.14_16<-ca.14_16 %>% 
    left_join(gopc.cost.14_16[c('reference.key.','gopc.episodes','gopc.cost')], by = 'reference.key.') %>% 
    replace_na(list(gopc.episodes = 0)) %>% 
    replace_na(list(gopc.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(gopc.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(gopc.cost/(follow.up.period/365.25)))

# 5.3.5 OP Specialist Outpatient Clinic (SOPC), non-Psychiatric --------
sopc.gen.cost.14_16<-op.serv.data.14_16 %>% 
    filter(grepl('^FM$|^NUR|^SOP', service.group.eis.)) %>% 
    filter(!op.eis.specialty. == "PSY") %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., appointment.date.yyyy.mm.dd., .keep_all = TRUE) %>% 
    group_by(reference.key.) %>% 
    summarise(sopc.gen.episodes = n(),
              sopc.gen.cost = sopc.gen.episodes*1190)

ca.14_16<-ca.14_16 %>% 
    left_join(sopc.gen.cost.14_16[c('reference.key.','sopc.gen.episodes','sopc.gen.cost')], by = 'reference.key.') %>% 
    replace_na(list(sopc.gen.episodes = 0)) %>% 
    replace_na(list(sopc.gen.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(sopc.gen.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(sopc.gen.cost/(follow.up.period/365.25)))

# 5.3.6 OP Specialist Outpatient Clinic (SOPC), Psychiatric --------
sopc.psy.cost.14_16<-op.serv.data.14_16 %>% 
    filter(grepl('^NUR|^SOP', service.group.eis.)) %>% 
    filter(op.eis.specialty. == "PSY") %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., appointment.date.yyyy.mm.dd., .keep_all = TRUE) %>% 
    group_by(reference.key.) %>% 
    summarise(sopc.psy.episodes = n(),
              sopc.psy.cost = sopc.psy.episodes*1190)

ca.14_16<-ca.14_16 %>% 
    left_join(sopc.psy.cost.14_16[c('reference.key.','sopc.psy.episodes','sopc.psy.cost')], by = 'reference.key.') %>% 
    replace_na(list(sopc.psy.episodes = 0)) %>% 
    replace_na(list(sopc.psy.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(sopc.psy.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(sopc.psy.cost/(follow.up.period/365.25)))

# 5.3.7 OP Specialist Outpatient Clinic (SOPC), Allied health --------
sopc.ah.cost.14_16<-op.serv.data.14_16 %>% 
    filter(grepl('^AHJ$|^AHOO$|^AHOP$|^AHTC$', service.type.code.eis.)) %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., appointment.date.yyyy.mm.dd., .keep_all = TRUE) %>% 
    group_by(reference.key.) %>% 
    summarise(sopc.ah.episodes = n(),
              sopc.ah.cost = sopc.ah.episodes*1190)

ca.14_16<-ca.14_16 %>% 
    left_join(sopc.ah.cost.14_16[c('reference.key.','sopc.ah.episodes','sopc.ah.cost')], by = 'reference.key.') %>% 
    replace_na(list(sopc.ah.episodes = 0)) %>% 
    replace_na(list(sopc.ah.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(sopc.ah.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(sopc.ah.cost/(follow.up.period/365.25)))

# 5.3.8 OP Community services, General ----
com.gen.cost.14_16<-op.serv.data.14_16 %>% 
    filter(grepl('^CEOL$|^CGAT$|^CVMO$|^ICMD$|^ICMN$', service.type.code.eis.) |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "CGAS" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "DH" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "MPMG" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "MPMR" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "PC" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "PCR" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "PCNC" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "CPLD" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "DPGC" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "DPGN" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "ICM" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "POAH") %>% 
    filter(!grepl('^VITO$|^VITP$|^VOS$|^DHOT$|^DHPT$', sub.specialty.opas.)) %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., appointment.date.yyyy.mm.dd., .keep_all = TRUE) %>% 
    group_by(reference.key.) %>% 
    summarise(com.gen.episodes = n(),
              com.gen.cost = com.gen.episodes*535)

ca.14_16<-ca.14_16 %>% 
    left_join(com.gen.cost.14_16[c('reference.key.','com.gen.episodes','com.gen.cost')], by = 'reference.key.') %>% 
    replace_na(list(com.gen.episodes = 0)) %>% 
    replace_na(list(com.gen.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(com.gen.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(com.gen.cost/(follow.up.period/365.25)))

# 5.3.8 OP Community services, Psychiatric ----
com.psy.cost.14_16<-op.serv.data.14_16 %>% 
    filter(grepl('^Community Psychiatric Services$|^PG Outreach$', service.group.eis.) |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "CSAC" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "DPPC" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "DPPN" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "PGER" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "PPSY" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "PSY" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "SPSY" |
               service.group.eis. == 'C - Community Service' & specialty.code.opas. == "STUC") %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., appointment.date.yyyy.mm.dd., .keep_all = TRUE) %>% 
    group_by(reference.key.) %>% 
    summarise(com.psy.episodes = n(),
              com.psy.cost = com.psy.episodes*1550)

ca.14_16<-ca.14_16 %>% 
    left_join(com.psy.cost.14_16[c('reference.key.','com.psy.episodes','com.psy.cost')], by = 'reference.key.') %>% 
    replace_na(list(com.psy.episodes = 0)) %>% 
    replace_na(list(com.psy.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(com.psy.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(com.psy.cost/(follow.up.period/365.25)))

# 5.3.9 OP Community services, Allied health ----
com.ah.cost.14_16<-op.serv.data.14_16 %>% 
    filter(grepl('^AHC$|^ICMA$', service.type.code.eis.) |
               service.group.eis. == 'C - Community Service' & sub.specialty.opas. == "OT" |
               service.group.eis. == 'C - Community Service' & sub.specialty.opas. == "PT" |
               service.group.eis. == 'C - Community Service' & sub.specialty.opas. == "DHOT" |
               service.group.eis. == 'C - Community Service' & sub.specialty.opas. == "DHPT" |
               service.group.eis. == 'C - Community Service' & sub.specialty.opas. == "VITO" |
               service.group.eis. == 'C - Community Service' & sub.specialty.opas. == "VITP" |
               service.group.eis. == 'C - Community Service' & sub.specialty.opas. == "VOS") %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., appointment.date.yyyy.mm.dd., .keep_all = TRUE) %>% 
    group_by(reference.key.) %>% 
    summarise(com.ah.episodes = n(),
              com.ah.cost = com.ah.episodes*1730)

ca.14_16<-ca.14_16 %>% 
    left_join(com.ah.cost.14_16[c('reference.key.','com.ah.episodes','com.ah.cost')], by = 'reference.key.') %>% 
    replace_na(list(com.ah.episodes = 0)) %>% 
    replace_na(list(com.ah.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(com.ah.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(com.ah.cost/(follow.up.period/365.25)))

# 5.3.10 OP Dressing/Injection ----
inj.cost.14_16<-op.serv.data.14_16 %>% 
    filter(grepl('Dressing|Injection', service.group.eis.))  %>% 
    arrange(reference.key., appointment.date.yyyy.mm.dd.) %>% 
    distinct(reference.key., appointment.date.yyyy.mm.dd., .keep_all = TRUE) %>% 
    group_by(reference.key.) %>% 
    summarise(inj.episodes = n(),
              inj.cost = inj.episodes*100)

ca.14_16<-ca.14_16 %>% 
    left_join(inj.cost.14_16[c('reference.key.','inj.episodes','inj.cost')], by = 'reference.key.') %>% 
    replace_na(list(inj.episodes = 0)) %>% 
    replace_na(list(inj.cost = 0))

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(inj.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(inj.cost/(follow.up.period/365.25)))

# 5.3.11 Combine all OP costs ----
ca.14_16<-ca.14_16 %>% 
    mutate(op.episodes = sopc.psy.episodes + sopc.gen.episodes + sopc.ah.episodes + gopc.episodes + pdh.episodes + gdh.episodes + rdh.episodes + 
               com.gen.episodes + com.psy.episodes + com.ah.episodes + inj.episodes,
           op.cost = sopc.psy.cost + sopc.gen.cost + sopc.ah.cost + gopc.cost + pdh.cost + gdh.cost + rdh.cost + 
               com.gen.cost + com.psy.cost + com.ah.cost + inj.cost) 

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    summarise(m.episodes.py = mean(op.episodes/(follow.up.period/365.25)),
              m.cost.py = mean(op.cost/(follow.up.period/365.25)))

# 5.4 Combine AE + IP + OP costs ----
ca.14_16<-ca.14_16 %>% 
    mutate(overall.cost = ae.cost + ip.cost + op.cost)

# 5.5 Calculate psychiatric-only cost ----
ca.14_16<-ca.14_16 %>% 
  mutate(op.psy.cost = sopc.psy.cost + pdh.cost + com.psy.cost,
         overall.psy.cost = ip.psy.cost + op.psy.cost)

ca.14_16 %>% 
    group_by(trd.status) %>% 
    filter(follow.up.period > 0) %>% 
    mutate(follow.up.period = follow.up.period/365.25) %>% 
    summarise(m.cost.py = mean(overall.cost/follow.up.period),
              m.cost = mean(overall.cost))
# 2014 cohort: 54558 vs 29737

# 6. GLM with negative binomial ----
## REMEMBER to exclude FU period <= 0 !!!
temp<-ca.14_16 %>% 
    filter(follow.up.period >= 0) %>% 
    mutate(follow.up.period = follow.up.period/365.25) %>% 
    mutate(overall.cost.py = overall.cost/follow.up.period) %>% 
  mutate(ip.cost.py = ip.cost/follow.up.period) %>% 
  mutate(op.cost.py = op.cost/follow.up.period) %>% 
  mutate(ae.cost.py = ae.cost/follow.up.period) %>% 
  mutate(ip.psy.cost.py = ip.psy.cost/follow.up.period) %>% 
  mutate(op.psy.cost.py = op.psy.cost/follow.up.period) %>% 
  mutate(overall.psy.cost.py = overall.psy.cost/follow.up.period)
  
# Unadjusted
glm1<-glm.nb(formula = overall.cost~offset(log(follow.up.period)), data = temp)
glm1<-glm.nb(formula = overall.cost.py~1, data = temp)

# Adjusted
glm1<-glm.nb(formula = ip.cost~offset(log(follow.up.period))
             + trd.status
             #+ post.match.comorbid
             + age.group
             + sex.
             + baseline.medical
             + post.index.medical
             , maxit=100
             , data = temp)

glm1<-glm.nb(formula = overall.psy.cost.py~trd.status
             #+ post.match.comorbid
             + age.group
             + sex.
             + baseline.medical
             + post.index.medical
             , maxit=100
             , data = temp)
summary(glm1)
result<-round(exp(cbind(coef(glm1), confint(glm1))), 3)
write.csv(result, 'result.csv')
deviance(glm1)/df.residual(glm1)

newdata<-read_excel('../Cost/Summary of cost.xlsx', sheet = 'new_data')
newdata<-data.frame(newdata)

tmp <- predict(glm1, newdata=newdata, se.fit=TRUE, type = "link")
temp<-t(rbind(rbind(exp(tmp$fit), exp(tmp$fit - 2*tmp$se.fit)), exp(tmp$fit + 2*tmp$se.fit)))
temp<-round(temp)
write.csv(temp, file = 'subgroup_cost.csv')

newdata<-read_excel('incidence_results.xlsx', sheet = 9)
newdata<-data.frame(newdata)




exp(coef(glm1)[1]) 
exp(confint(glm1)[1]) 
exp(confint(glm1)[2]) # Unadjusted

exp(coef(glm1)[1]) 
exp(confint(glm1)[1,1]) 
exp(confint(glm1)[1,2]) # NTRD cost ppy = $22284 [21655, 22936]

exp(coef(glm1)[1]+coef(glm1)[2]) 
exp(confint(glm1)[1,1]+confint(glm1)[2,1]) 
exp(confint(glm1)[1,2]+confint(glm1)[2,2]) # TRD cost ppy = $43355 [39736, 47347]

exp(coef(glm1)[1]+coef(glm1)[3]) 
exp(confint(glm1)[1,1]+confint(glm1)[3,1]) 
exp(confint(glm1)[1,2]+confint(glm1)[3,2]) # NTRD comorbid cost ppy = $70482 [64052, 77660]

exp(coef(glm1)[1]+coef(glm1)[2]+coef(glm1)[3]) 
exp(confint(glm1)[1,1]+confint(glm1)[2,1]+confint(glm1)[3,1]) 
exp(confint(glm1)[1,2]+confint(glm1)[2,2]+confint(glm1)[3,2]) # TRD comorbid cost ppy = $137129 [117534, 160313]

# 7. Plot graphs ----
cost.summary<-read_excel('../Cost/summary of cost.xlsx', sheet = 3)
cost.summary<-data.frame(cost.summary)
cost.summary<-cost.summary %>% 
  arrange(Type, State, Setting) %>% 
  mutate(Setting = factor(Setting, levels = c('Accident & Emergency', 'Outpatient', 'Inpatient','Overall')),
         State = factor(State, levels = c('NTRD','TRD','NTRD-comorbid','TRD-comorbid')))
cost.summary.overall<-cost.summary %>% filter(Setting == 'Overall')
cost.summary.each<-cost.summary %>% filter(!Setting == 'Overall')

p<-ggplot() +
  geom_bar(data=cost.summary.each, aes(x = State, y = Cost, fill=Setting), stat ="identity", position="stack", width =0.8)+
  geom_errorbar(data=cost.summary.overall, aes(x = State, y = Cost, ymin=Cost.lower, ymax=Cost.upper), width=.3, size = 0.3) +
  geom_text(data=cost.summary.overall, aes(x = State, y = Cost.upper), 
            label = format(cost.summary.overall$Cost, big.mark = ",", scientific = FALSE), 
            size = 3, position = position_dodge(width = 0.9), vjust = -1) +
  scale_y_continuous(name = "Cost per patient year (HKD)", label = scales::comma, limits = c(0, 175000)) +
  scale_x_discrete(name = 'Markov health state') +
  facet_grid(~Type)+
  scale_fill_brewer(palette = "PuBuGn") +
  theme_pubclean()
p

# 8.1 Table one (After matching) ----
temp<-dx.14_16.inc.com %>% 
  filter(reference.date. < first.dep.date) %>% 
  select(reference.key., comorbidity.name) %>% 
  distinct() %>% 
  dummy_cols(select_columns = 'comorbidity.name')
colnames(temp)<-c('reference.key.','comorbidity.name',
                  'ADHD','AIDS','Anxiety','Tumors','Autism','Bipolar','CVD','CKD','COPD','CHF',
                  'CTD','Dementia','DM','Eating','Epilepsy','Hemiplegia','Leukemia','Liver','Lymphoma',
                  'MI','OCD','PUD','PVD','Personality','Psychosis','Substance','Suicidal')

temp<-temp %>% select(-comorbidity.name) %>% group_by(reference.key.) %>% summarise_all(funs(max(as.character(.)))) 

temp<-ca.14_16 %>% 
  filter(follow.up.period > 0) %>% 
  left_join(temp) %>% 
  mutate_at(95:121, ~replace_na(.,0))

V<-c('age','age.group','sex.','baseline.somatic',
     'AIDS','CVD','CKD','COPD','CHF','CTD','Dementia','DM',
     'Hemiplegia','Leukemia','Liver','Lymphoma','MI','PUD','PVD','Tumors',
     'baseline.mental',
     'ADHD','Anxiety','Eating','Epilepsy','Autism','Bipolar','OCD','Personality','Psychosis','Substance','Suicidal',
     'follow.up.period')
cV<-c('age.group','sex.','baseline.somatic','baseline.mental',
      'ADHD','AIDS','Anxiety','Tumors','Autism','Bipolar','CVD','CKD','COPD','CHF',
      'CTD','Dementia','DM','Eating','Epilepsy','Hemiplegia','Leukemia','Liver','Lymphoma',
      'MI','OCD','PUD','PVD','Personality','Psychosis','Substance','Suicidal')

Table.1<-CreateTableOne(vars=V, 
                        strata="trd.status",
                        data=temp, factorVars= cV)
Table.1<-print(Table.1,quote=FALSE, noSpaces=TRUE, smd = T)
write.csv(Table.1, 'Table.1.csv')

# 8.2 Table one (Before matching) ----
temp<-dx.14_16.inc.com %>% 
  filter(reference.date. < first.dep.date) %>% 
  dplyr::select(reference.key., comorbidity.name) %>% 
  distinct() %>% 
  dummy_cols(select_columns = 'comorbidity.name')
colnames(temp)<-c('reference.key.','comorbidity.name',
                  'ADHD','AIDS','Anxiety','Tumors','Autism','Bipolar','CVD','CKD','COPD','CHF',
                  'CTD','Dementia','DM','Eating','Epilepsy','Hemiplegia','Leukemia','Liver','Lymphoma',
                  'MI','OCD','PUD','PVD','Personality','Psychosis','Substance','Suicidal')

temp<-temp %>% 
  dplyr::select(-comorbidity.name) %>%
  group_by(reference.key.) %>% 
  summarise_all(funs(max(as.character(.)))) 

temp<-depression.14_16.inc %>% 
  filter(follow.up.period > 0) %>% 
  left_join(temp) %>% 
  mutate_at(49:75, ~replace_na(.,0))

V<-c('age','age.group','sex.','baseline.somatic',
     'AIDS','CVD','CKD','COPD','CHF','CTD','Dementia','DM',
     'Hemiplegia','Leukemia','Liver','Lymphoma','MI','PUD','PVD','Tumors',
     'baseline.mental',
     'ADHD','Anxiety','Eating','Epilepsy','Autism','Bipolar','OCD','Personality','Psychosis','Substance','Suicidal',
     'follow.up.period')
cV<-c('age.group','sex.','baseline.somatic','baseline.mental',
      'ADHD','AIDS','Anxiety','Tumors','Autism','Bipolar','CVD','CKD','COPD','CHF',
      'CTD','Dementia','DM','Eating','Epilepsy','Hemiplegia','Leukemia','Liver','Lymphoma',
      'MI','OCD','PUD','PVD','Personality','Psychosis','Substance','Suicidal')

temp2<-temp %>% filter(age <= 25)

Table.1<-CreateTableOne(vars=V, 
                        strata="trd.status",
                        data=temp2, factorVars= cV)
Table.1<-print(Table.1,quote=FALSE, noSpaces=TRUE, smd = T)
write.csv(Table.1, 'Table.1.csv')



library(survival)
library(survminer)
library(tidyverse)
library(dplyr)
library(patchwork)
library(rstanarm)
library(brms)
library(rstan)
library(lme4)

setwd("~/STA723/Case Study 5/Case_Study_5")
data <- read.csv(file = 'data/case_01_main.csv')
data <- data[,-1]
colnames(data) <- c("PATID","treatment","time","event")


# EDA 

test <- data %>% 
  group_by(treatment, event) %>% 
  summarise(mean(time))

tmp1 <- data %>% 
  filter(event == 1)

tmp2 <- data %>% 
  filter(event == 0)


p1 <- tmp1 %>%
  ggplot(mapping = aes(x = treatment, y = time)) +
  geom_boxplot(color = "black", fill = "red", alpha = .8) +
  labs(x = "Treatments", y = "Relapse Time (Days)") +
  ggtitle("Relapse Time by Treatments") +
  theme_bw(base_size = 16)

p2 <- tmp2 %>%
  ggplot(mapping = aes(x = treatment, y = time)) +
  geom_boxplot(color = "black", fill = "yellow", alpha = .8) +
  labs(x = "Treatments", y = "Censored Time (Days)") +
  ggtitle("Censored Time by Treatments") +
  theme_bw(base_size = 16)

p1 + p2

# the relapse time of XR-NTX is a little strange

time1 <- data %>% 
  filter(treatment == "XR-NTX") %>% 
  filter(event == 1)


ggplot(time, mapping = aes(x = time)) +
  geom_histogram(binwidth = 5, fill = "brown") +
  labs(x = "Frequencies", y = "Relapse Time (Days)") +
  ggtitle("Relapse Time of XR-NTX (Exclude Censored)") +
  theme_bw(base_size = 16)


# fit survival data using the Kaplan-Meier method"BUP-NX"

surv_object <- Surv(time = data$time, event = data$event)
surv_object

fit1 <- survfit(surv_object ~ treatment, data = data)
summary(fit1)

ggsurvplot(fit1, data = data, pval = TRUE)


# import other data into a single one

setwd("~/STA723/Case Study 5/Case_Study_5/data")
vis <- read.csv(file = 'case_01_vis.csv')
qol <- read.csv(file = 'case_01_qol.csv')
qlp <- read.csv(file = 'case_01_qlp.csv')
mhx <- read.csv(file = 'case_01_mhx.csv')
inf <- read.csv(file = 'case_01_inf.csv')
fnd <- read.csv(file = 'case_01_fnd.csv')
dsm <- read.csv(file = 'case_01_dsm.csv')
dem <- read.csv(file = 'case_01_dem.csv')
asp <- read.csv(file = 'case_01_asp.csv')
asi <- read.csv(file = 'case_01_asi.csv')
asf <- read.csv(file = 'case_01_asf.csv')
ase <- read.csv(file = 'case_01_ase.csv')
asd <- read.csv(file = 'case_01_asd.csv')

test <- left_join(data, 
                  vis, by = "PATID")
test <- left_join(test, 
                  qol, by = "PATID")
test <- left_join(test, 
                  qlp, by = "PATID")
test <- left_join(test, 
                  mhx, by = "PATID")
test <- left_join(test, 
                  inf, by = "PATID")
test <- left_join(test, 
                  fnd, by = "PATID")
test <- left_join(test, 
                  dsm, by = "PATID")
test <- left_join(test, 
                  dem, by = "PATID")
test <- left_join(test, 
                  asp, by = "PATID")
test <- left_join(test, 
                  asi, by = "PATID")
test <- left_join(test, 
                  asf, by = "PATID")
test <- left_join(test, 
                  ase, by = "PATID")
test <- left_join(test, 
                  asd, by = "PATID")
test2 <- test

test <- test2[-c(49,64,89,93,97,207,209,224,235,304,419,494,496,502,516),]
v <- colnames(test)

# Fit a Cox proportional hazards model

data1 <- test %>% 
  filter(treatment == "BUP-NX")

surv_object1 <- Surv(time = data1$time, event = data1$event)

# QLHMLESS + ADOPI30D + VIBMI +  QLANXDEP

fit.coxph <- coxph(surv_object1 ~  VIBMI +  QLANXDEP + QLHMLESS + ADOPI30D + MHOPIWDL + ADMDPLFT + ADTHC30D + ADTHCLFT, 
                   data = data1)

p3 <- ggforest(fit.coxph, data = data1, main = "BUP-NX Hazard ratio",
               fontsize = 1) 




# Fit a Cox proportional hazards model

data2 <- test %>% 
  filter(treatment == "XR-NTX")

surv_object2 <- Surv(time = data2$time, event = data2$event)

# QLHMLESS + ADOPI30D + MHOPIWDL + ADMDPLFT + ADTHC30D + ADTHCLFT

fit.coxph2 <- coxph(surv_object2 ~ VIBMI +  QLANXDEP + QLHMLESS + ADOPI30D +  MHOPIWDL + ADMDPLFT + ADTHC30D + ADTHCLFT, 
                   data = data2)

p4 <- ggforest(fit.coxph2, data = data2, main = "XR-NTX Hazard ratio",
               fontsize = 1)

p3 + p4


################ per protocol experiment


# Fit a Cox proportional hazards model

dat1 <- test %>% 
  mutate(induction_failure = if_else(is.na(INFAILDT), 0, 1))%>% 
  filter(induction_failure == 0) %>% 
  filter(treatment == "BUP-NX")

surv_object1 <- Surv(time = dat1$time, event = dat1$event)

# QLHMLESS + ADOPI30D + VIBMI +  QLANXDEP

fit.coxph <- coxph(surv_object1 ~  VIBMI +  QLANXDEP + QLHMLESS + ADOPI30D + MHOPIWDL + ADMDPLFT + ADTHC30D + ADTHCLFT, 
                   data = dat1)

p3 <- ggforest(fit.coxph, data = dat1, main = "Per Protocol BUP-NX Hazard ratio",
               fontsize = 1) 




# Fit a Cox proportional hazards model

dat2 <- test %>% 
  mutate(induction_failure = if_else(is.na(INFAILDT), 0, 1))%>% 
  filter(induction_failure == 0) %>% 
  filter(treatment == "XR-NTX")

surv_object2 <- Surv(time = dat2$time, event = dat2$event)

# QLHMLESS + ADOPI30D + MHOPIWDL + ADMDPLFT + ADTHC30D + ADTHCLFT

fit.coxph2 <- coxph(surv_object2 ~ VIBMI +  QLANXDEP + QLHMLESS + ADOPI30D +  MHOPIWDL + ADMDPLFT + ADTHC30D + ADTHCLFT, 
                    data = dat2)

p4 <- ggforest(fit.coxph2, data = dat2, main = "Per Protocol XR-NTX Hazard ratio",
               fontsize = 1)

p4







# Bayesian logistic model

# EDA
tmp <- test %>% 
  mutate(induction_failure = if_else(is.na(INFAILDT), 0, 1))

tmp$induction_failure <- as.factor(tmp$induction_failure)


ttt <- tmp %>% 
  group_by(induction_failure, treatment) %>% 
  dplyr::summarise(n=n())


len <- c(268,203,12,72)

ggplot(data=ttt, aes(x=treatment, y=n, fill=induction_failure)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=len), vjust=-0.7, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  labs(x = "Treatments", y = "Total Number") +
  ggtitle("Total Number of Induction Failure by Treatments") +
  theme_minimal()


# proportion of successful induction 

dat1 <- test %>% 
  mutate(induction_failure = if_else(is.na(INFAILDT), 0, 1))%>% 
  filter(treatment == "BUP-NX")

# ADOPI30D

fit_logit <- brm(formula = induction_failure ~ ADOPI30D + MHOPIWDL + VIBPSYS1+ QLHMLESS,  
                 data = dat1, 
                  family = bernoulli(link = "logit"),
                  warmup = 500, 
                  iter = 4000, 
                  chains = 2, 
                  inits= "0", 
                  cores= 2,
                  seed = 123)


summary(fit_logit)

stanplot(fit_logit, 
         type = "trace")

pp_check(fit_logit, nsamples = 500)

predict(fit_logit)




dat2 <- test %>% 
  mutate(induction_failure = if_else(is.na(INFAILDT), 0, 1))%>% 
  filter(treatment == "XR-NTX")
 
# ADOPI30D + MHOPIWDL + VIBPSYS1+ QLHMLESS

fit_logit2 <- brm(formula = induction_failure ~ ADOPI30D + MHOPIWDL + VIBPSYS1+ QLHMLESS,  
                 data = dat2, 
                 family = bernoulli(link = "logit"),
                 warmup = 500, 
                 iter = 4000, 
                 chains = 2, 
                 inits= "0", 
                 cores= 2,
                 seed = 123)


summary(fit_logit2)

stanplot(fit_logit2, 
         type = "trace")

pp_check(fit_logit2, nsamples = 500)



# Bayesian logistic model

tmp <- test %>% 
  mutate(induction_failure = if_else(is.na(INFAILDT), 0, 1))

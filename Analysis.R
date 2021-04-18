library(survival)
library(survminer)
library(tidyverse)
library(dplyr)
library(ggplot2)

setwd("~/STA723/Case Study 5/Case_Study_5")
data <- read.csv(file = 'data/case_01_main.csv')


# EDA 

data_BUP <- data %>% 
  group_by(TRTSHOWN, event) %>% 
  summarise(mean(time))


# Fit survival data using the Kaplan-Meier method

surv_object <- Surv(time = data$time, event = data$event)
surv_object

fit1 <- survfit(surv_object ~ TRTSHOWN, data = data)
summary(fit1)

ggsurvplot(fit1, data = data, pval = TRUE)


# import other data into a single one

setwd("~/STA723/Case Study 5/Case_Study_5/data")
vis <- read.csv(file = 'case_01_vis.csv')
qol <- read.csv(file = 'case_01_qol.csv')
qlp <- read.csv(file = 'case_01_qlp.csv')
fnd <- read.csv(file = 'case_01_fnd.csv')
asp <- read.csv(file = 'case_01_asp.csv')
asf <- read.csv(file = 'case_01_asf.csv')
ase <- read.csv(file = 'case_01_ase.csv')


test <- left_join(data, 
                  vis, by = "PATID")
test <- left_join(test, 
                  qol, by = "PATID")
test <- left_join(test, 
                  qlp, by = "PATID")

# Fit a Cox proportional hazards model
fit.coxph <- coxph(surv_object ~ VIBMI + VIHGTIN + QLHMLESS + QLMOBIL + QLSLFCAR + QLACTIVE + QLPAIN + QLANXDEP, 
                   data = test)

ggforest(fit.coxph, data = ovarian)




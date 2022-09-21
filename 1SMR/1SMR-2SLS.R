# Overview
# This file contains the R script of the main analysis of one-sample Mendelian randomization (1SMR) with two-stage least squares (2SLS) instrumental variable analyses
# Here we take accelerometer-derived (AcD) sleep duration vs HbA1c as an example (R codes of other sleep traits are similar)
# For non-fasting glucose, fasting time and dilution factor are additionally adjusted for 
# Sensitivity analyses excluding participants with diabetes (HbA1c >= 48 mmol/mol or diabetic status defined by the Eastwood algorithm )
# For HbA1c, it is needed to divide the estimates by 0.14 log mmol/mol to obtain the SD estimates; for non-fasting glucose, it is needed to divide the estimates by 0.16 log mmol/l to obtain the SD unit

# Install related R packages
library(AER)
library(ff)

# Read the sub-sample of UK Biobank data with AcD sleep measures
data<- fread(file="data.csv",
             fill=T, header = T)

#continuous covariates (age at recruitment + principle components (PC) 1 - 40) 
rc_numeric<-cbind(data$age_recruitment,
                  data$PC01,data$PC02,data$PC03,data$PC04,data$PC05,data$PC06,data$PC07,data$PC08,data$PC09,data$PC010,
                  data$PC011,data$PC012,data$PC013,data$PC014,data$PC015,data$PC016,data$PC017,data$PC018,data$PC019,data$PC020,
                  data$PC021,data$PC022,data$PC023,data$PC024,data$PC025,data$PC026,data$PC027,data$PC028,data$PC029,data$PC030,
                  data$PC031,data$PC032,data$PC033,data$PC034,data$PC035,data$PC036,data$PC037,data$PC038,data$PC039,data$PC040)

# AcD sleep duration vs HbA1c (log unit)
# Adjusted for sex, age, assessment centre, chip, PC1 - PC40
mr<-ivreg(log_hba1c ~ AcD_sleepduration  + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric 
          |AcD_sleepduration_allelescore + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric,  
          data = data)
                  
AcD_sleepduration_hba1c<-coef(summary(mr))

AcD_sleepduration_hba1c_2sls_b<-AcD_sleepduration_hba1c[2,1]
AcD_sleepduration_hba1c_2sls_lci<-AcD_sleepduration_hba1c[2,1] - 1.96 * AcD_sleepduration_hba1c[2,2]
AcD_sleepduration_hba1c_2sls_uci<-AcD_sleepduration_hba1c[2,1] + 1.96 * AcD_sleepduration_hba1c[2,2]
AcD_sleepduration_hba1c_2sls_p<-AcD_sleepduration_hba1c[2,4]

n_AcD_sleepduration_hba1c<-summary(mr)$'df'[1] + summary(mr)$'df'[2]







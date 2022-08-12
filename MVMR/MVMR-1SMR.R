# Overview
# This file contains the R script of multivariable Mendelian randomization in one-sample setting (MVMR-1SMR)
# Here we take accelerometer-derived (AD) sleep duration + sleep efficiency vs HbA1c as an example (R codes of other sleep traits are similar)
# For non-fasting glucose, fasting time and dilution factor are additionally adjusted for 
# For HbA1c, it is needed to divide the estimates by 0.14 log mmol/mol to obtain the SD estimates; for non-fasting glucose, it is needed to divide the estimates by 0.16 log mmol/l to obtain the SD unit

# Install related R packages
library(ff)
library(AER)
library(ivpack)

rm(list = ls())
data<- fread(file="data.csv",
             fill=T, header = T)
# efficiency is a % unit, it shows 0-1 in the dataset, it should be converted into % 0-100 for better interpretation
data$efficiency <- data$acc_sleep_eff_AD_mn * 100

# continuous covariates (baseline age + principle compoments (PC) 1 - 40) 
rc_numeric<-cbind(data$age_recruitment,
                  data$PC1,data$PC2,data$PC3,data$PC4,data$PC5,data$PC6,data$PC7,data$PC8,data$PC9,data$PC10,
                  data$PC11,data$PC12,data$PC13,data$PC14,data$PC15,data$PC16,data$PC17,data$PC18,data$PC19,data$PC20,
                  data$PC21,data$PC22,data$PC23,data$PC24,data$PC25,data$PC26,data$PC27,data$PC28,data$PC29,data$PC30,
                  data$PC31,data$PC32,data$PC33,data$PC34,data$PC35,data$PC36,data$PC37,data$PC38,data$PC39,data$PC40)


# 1st stage to predict  AD sleepduration 
# (linear regression: adjusted for sex, age, assessment centre, chip, 40 genetic principal components)
data$pre_ad_sleepduration <- predict(glm(ad_sleepduration ~ ad_sleepduration_allelescore + efficiency_allelescore + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric,
                                         na.action = na.exclude, data=data))

# 1st stage to predict efficiency  
# (linear regression: adjusted for sex, age, assessment centre, chip, 40 genetic principal components)
data$pre_efficiency <- predict(glm(efficiency ~ ad_sleepduration_allelescore + efficiency_allelescore + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric,
                                   na.action = na.exclude, data=data))


# 2nd stage
# (linear regression: adjusted for sex, age, assessment centre, chip, 40 genetic principal components)
twoStage <- glm(log_hba1c ~ pre_ad_sleepduration + pre_efficiency  + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric, 
                data=data)

ad_sleepduration_efficiency_hba1c<-coef(summary(twoStage))
n_ad_sleepduration_efficiency_hba1c_mvmr<-summary(twoStage)$'df'[1] + summary(twoStage)$'df'[2]

# AD sleep duration
ad_sleepduration_hba1c_mvmr_b<-ad_sleepduration_efficiency_hba1c[2,1]
ad_sleepduration_hba1c_mvmr_lci<-(ad_sleepduration_efficiency_hba1c[2,1] - 1.96 * ad_sleepduration_efficiency_hba1c[2,2])
ad_sleepduration_hba1c_mvmr_uci<-(ad_sleepduration_efficiency_hba1c[2,1] + 1.96 * ad_sleepduration_efficiency_hba1c[2,2])
ad_sleepduration_hba1c_mvmr_p<-ad_sleepduration_efficiency_hba1c[2,4]

# Efficiency
efficiency_hba1c_mvmr_b<-ad_sleepduration_efficiency_hba1c[3,1]
efficiency_hba1c_mvmr_lci<-(ad_sleepduration_efficiency_hba1c[3,1] - 1.96 * ad_sleepduration_efficiency_hba1c[3,2])
efficiency_hba1c_mvmr_uci<-(ad_sleepduration_efficiency_hba1c[3,1] + 1.96 * ad_sleepduration_efficiency_hba1c[3,2])
efficiency_hba1c_mvmr_p<-ad_sleepduration_efficiency_hba1c[3,4]


# Sanderson–Windmeijer conditional F statistics of the allelescore allele score https://doi.org/10.1093/ije/dyy262
# select data to calcuate the conditional F statistics 
data_fstatistic<-select(data, ad_sleepduration, efficiency, ad_sleepduration_allelescore,efficiency_allelescore,
                        sex, assessment_centre, chip, age_recruitment,
                        PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                        PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,
                        PC21,PC22,PC23,PC24,PC25,PC26,PC27,PC28,PC29,PC30,
                        PC31,PC32,PC33,PC34,PC35,PC36,PC37,PC38,PC39,PC40)

data_fstatistic<-data_fstatistic[complete.cases(data_fstatistic), ]

# continuous covariates (baseline age + principle compoments (PC) 1 - 40) 
rc_numeric<-cbind(data_fstatistic$age_recruitment,
                  data_fstatistic$PC1,data_fstatistic$PC2,data_fstatistic$PC3,data_fstatistic$PC4,data_fstatistic$PC5,data_fstatistic$PC6,data_fstatistic$PC7,data_fstatistic$PC8,data_fstatistic$PC9,data_fstatistic$PC10,
                  data_fstatistic$PC11,data_fstatistic$PC12,data_fstatistic$PC13,data_fstatistic$PC14,data_fstatistic$PC15,data_fstatistic$PC16,data_fstatistic$PC17,data_fstatistic$PC18,data_fstatistic$PC19,data_fstatistic$PC20,
                  data_fstatistic$PC21,data_fstatistic$PC22,data_fstatistic$PC23,data_fstatistic$PC24,data_fstatistic$PC25,data_fstatistic$PC26,data_fstatistic$PC27,data_fstatistic$PC28,data_fstatistic$PC29,data_fstatistic$PC30,
                  data_fstatistic$PC31,data_fstatistic$PC32,data_fstatistic$PC33,data_fstatistic$PC34,data_fstatistic$PC35,data_fstatistic$PC36,data_fstatistic$PC37,data_fstatistic$PC38,data_fstatistic$PC39,data_fstatistic$PC40)


# conditional F statistics for AD sleep duration
reg1 <- ivreg(ad_sleepduration ~ efficiency + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric|ad_sleepduration_allelescore + efficiency_allelescore + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric, data = data_fstatistic)
data_fstatistic$residuals <- residuals(reg1)

res1 <-(lm(residuals ~ sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric, data = data_fstatistic))
res2 <- (lm(residuals ~ ad_sleepduration_allelescore + efficiency_allelescore + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric, data = data_fstatistic))
F_ad_sleepduration <- anova(res1, res2)$F[2]*(anova(res1,res2)$Df[2]/(anova(res1,res2)$Df[2]-1))



# Conditional F statistics for efficiency
reg1 <- ivreg(efficiency ~ ad_sleepduration + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric|ad_sleepduration_allelescore + efficiency_allelescore + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric, data = data_fstatistic)
data_fstatistic$residuals <- residuals(reg1)

res1 <-(lm(residuals ~ sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric, data = data_fstatistic))
res2 <- (lm(residuals ~ ad_sleepduration_allelescore + efficiency_allelescore + sex + as.factor(assessment_centre) + as.factor(chip) + rc_numeric, data = data_fstatistic))

F_efficiency <- anova(res1, res2)$F[2]*(anova(res1,res2)$Df[2]/(anova(res1,res2)$Df[2]-1))




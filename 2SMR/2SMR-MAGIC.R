# Overview
# This file contains the R script of the two-sample Mendelian randomization analyses (2SMR) applying to the MAGIC (2SMR-MAGIC)
# (inverse-variance weighted (IVW),  weighted median (WM),  MR-Egger, MR-Egger_SiMEX)
# Here we take accelerometer-derived (AcD) sleep duration vs HbA1c as an example  (R codes of other sleep traits and with glucose are similar)
# For HbA1c, it is needed to divide the estimates by 0.41% to obtain the SD estimates; for non-fasting glucose, it is needed to divide the estimates by 0.84 mmol/l to obtain the SD unit



# Install related R packages
library(TwoSampleMR)
library(dplyr)
library(simex)


# Exposure:  The summary statistics of AcD sleep duration and other AcD sleep traits were from the GWAS (Jones S et al, PMID:0952852)
## The discovery GWAS presented the SNP-AcD sleep duration association in SD unit, where 1SD is equal to 0.86 hours. Therefore, We transformed it into per hour change unit.
AcD_sleepduration <- read.csv("AcD_sleepduration.csv")


# Outcomes: SNPs - HbA1c associations are from the GWAS (Chen J et al 2021, PMID: 34059833), with the "ebi-a-GCST90002244" id in the MRbase
# SNPs - fasting-glucose associations are from the same GWAS in above, with the "ebi-a-GCST90002232" id in the MRbase
hba1c <- extract_outcome_data(AcD_sleepduration$SNP, c("ebi-a-GCST90002244"))  


# To harmonise 
harmonised_AcD_sleepduration_hba1c <- harmonise_data(AcD_sleepduration, hba1c, action=2)

# To run MR analyses
AcD_sleepduration_hba1c_dat_mr <- mr(harmonised_AcD_sleepduration_hba1c) 
nsnp<-AcD_sleepduration_hba1c_dat_mr[1,6]

# IVW
ivw_b<-AcD_sleepduration_hba1c_dat_mr[3,7]
ivw_lci<-ivw_b - 1.96*AcD_sleepduration_hba1c_dat_mr[3,8]
ivw_uci<-ivw_b + 1.96*AcD_sleepduration_hba1c_dat_mr[3,8]
ivw_p<-AcD_sleepduration_hba1c_dat_mr[3,9]

# WM
wm_b<-AcD_sleepduration_hba1c_dat_mr[2,7]
wm_lci<-wm_b - 1.96*AcD_sleepduration_hba1c_dat_mr[2,8]
wm_uci<-wm_b + 1.96*AcD_sleepduration_hba1c_dat_mr[2,8]
wm_p<-AcD_sleepduration_hba1c_dat_mr[2,9]

# MR-Egger
egger_b<-AcD_sleepduration_hba1c_dat_mr[1,7]
egger_lci<-egger_b - 1.96*AcD_sleepduration_hba1c_dat_mr[1,8]
egger_uci<-egger_b + 1.96*AcD_sleepduration_hba1c_dat_mr[1,8]
egger_p<-AcD_sleepduration_hba1c_dat_mr[1,9]
egger_intercept_p<-mr_pleiotropy_test(harmonised_AcD_sleepduration_hba1c)[1,7]
mr_heterogeneity_ivw<-mr_heterogeneity(harmonised_AcD_sleepduration_hba1c)[2,8]

# MR-Egger with Simex
BetaYG <- harmonised_AcD_sleepduration_hba1c$beta.outcome
BetaXG <- harmonised_AcD_sleepduration_hba1c$beta.exposure
seBetaYG <- harmonised_AcD_sleepduration_hba1c$se.outcome
seBetaXG <- harmonised_AcD_sleepduration_hba1c$se.exposure
Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}
set.seed(12345)
Fit2 = lm(BetaYG~BetaXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)
IsqGX = Isq(BetaXG,seBetaXG)  
mod.sim <- simex(Fit2,B=1000,
                 measurement.error = seBetaXG,
                 SIMEXvariable="BetaXG",fitting.method ="quad",asymptotic="FALSE")
coef<-summary(mod.sim)$coefficients$jackknife
egger_simex_b <- coef[2,1] 
egger_simex_lci <- coef[2,1] - 1.96*coef[2,2]
egger_simex_uci <- coef[2,1] + 1.96*coef[2,2]
egger_simex_p <- coef[2,4]
egger_simex_intercept_p <- coef[1,4]



# Overview
# This file contains the R script for the two-sample Mendelian randomization analyses (2SMR) applying to the UK Biobank (2SMR-UKB)
# Here we take accelerometer-derived (AD) sleep duration vs HbA1c as an example  (R codes of other sleep traits and with glucose are similar)
# Sensitivity analyses excluding participants with diabetes (HbA1c >= 48 mmol/mol or diabetic status defined by the Eastwood algorithm )
## For HbA1c, it is needed to divide the estimates by 0.15 log mmol/mol to obtain the SD estimates; for non-fasting glucose, it is needed to divide the estimates by 0.18 mmol/l to obtain the SD unit

# Install related R packages
library(stringr)
library(ff)
library(dplyr)
library(simex)
library(TwoSampleMR)

# Step 1: prepare for the individual data to generate the summary statistics of SNP-Outcome association in the sub-sample of UKB without AcD sleep data
## Import the individual genotype data of AcD sleep duration from the UK Biobank (We extracted the individual genotype data of AcD sleep duration previously)
## AcD sleep duration genotype data
AcD_sleepduration_genotype<-read.csv("AcD_sleepduration_genotype.csv")

## Import the UKB data of HbA1c/glucose and other related covariates (age, sex, assessment centre, genotyping chip, PC1 - PC40) (would additionally include fasting time and dilution factor for non-fasting glucose)
data_ukb<- read.csv("ukb.csv")
data_ukb_hba1c<-select(data_ukb, projectid, log_hba1c, 
                      sex, age_recruitment, assessment_centre, chip,
                      PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                      PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,
                      PC21,PC22,PC23,PC24,PC25,PC26,PC27,PC28,PC29,PC30,
                      PC31,PC32,PC33,PC34,PC35,PC36,PC37,PC38,PC39,PC40)

## Merge the two datasets  
data2_ukb_hba1c<-merge(data_ukb_hba1c, AcD_sleepduration_genotype, by = "projectid")
data2_ukb_hba1c<-na.omit(data2_ukb_hba1c)


## Exclude participants with AcD sleep data 
data_AcD_gwas<-read.csv("data_AcD_gwas.csv") # dataset with AcD sleep traits
sample2_AcD_sleepduration_hba1c <- data2_ukb_hba1c[ !(data2_ukb_hba1c$projectid %in% c(data_AcD_gwas$projectid)), ]




# Step 2: generate the summary statistics for SNP-Outcome association in the sub-sample of UKB without AcD sleep data
n <- nrow(sample2_AcD_sleepduration_hba1c)
G1.2 <- as.matrix(sample2_AcD_sleepduration_hba1c[,-(1:46)])
Y <- as.matrix(sample2_AcD_sleepduration_hba1c[,2])

## covariates
sex<-as.matrix(sample2_AcD_sleepduration_hba1c[,3])
age<-as.matrix(sample2_AcD_sleepduration_hba1c[,4])
centre<-as.matrix(sample2_AcD_sleepduration_hba1c[,5])
chip<-as.matrix(sample2_AcD_sleepduration_hba1c[,6])

pc1 <- as.matrix(sample2_AcD_sleepduration_hba1c[,7])     
pc2 <- as.matrix(sample2_AcD_sleepduration_hba1c[,8]) 
pc3 <- as.matrix(sample2_AcD_sleepduration_hba1c[,9]) 
pc4 <- as.matrix(sample2_AcD_sleepduration_hba1c[,10]) 
pc5 <- as.matrix(sample2_AcD_sleepduration_hba1c[,11]) 
pc6 <- as.matrix(sample2_AcD_sleepduration_hba1c[,12]) 
pc7 <- as.matrix(sample2_AcD_sleepduration_hba1c[,13]) 
pc8 <- as.matrix(sample2_AcD_sleepduration_hba1c[,14]) 
pc9 <- as.matrix(sample2_AcD_sleepduration_hba1c[,15]) 
pc10 <- as.matrix(sample2_AcD_sleepduration_hba1c[,16]) 

pc11 <- as.matrix(sample2_AcD_sleepduration_hba1c[,17])     
pc12 <- as.matrix(sample2_AcD_sleepduration_hba1c[,18]) 
pc13 <- as.matrix(sample2_AcD_sleepduration_hba1c[,19]) 
pc14 <- as.matrix(sample2_AcD_sleepduration_hba1c[,20]) 
pc15 <- as.matrix(sample2_AcD_sleepduration_hba1c[,21]) 
pc16 <- as.matrix(sample2_AcD_sleepduration_hba1c[,22]) 
pc17 <- as.matrix(sample2_AcD_sleepduration_hba1c[,23]) 
pc18 <- as.matrix(sample2_AcD_sleepduration_hba1c[,24]) 
pc19 <- as.matrix(sample2_AcD_sleepduration_hba1c[,25]) 
pc20 <- as.matrix(sample2_AcD_sleepduration_hba1c[,26]) 

pc21 <- as.matrix(sample2_AcD_sleepduration_hba1c[,27])     
pc22 <- as.matrix(sample2_AcD_sleepduration_hba1c[,28]) 
pc23 <- as.matrix(sample2_AcD_sleepduration_hba1c[,29]) 
pc24 <- as.matrix(sample2_AcD_sleepduration_hba1c[,30]) 
pc25 <- as.matrix(sample2_AcD_sleepduration_hba1c[,31]) 
pc26 <- as.matrix(sample2_AcD_sleepduration_hba1c[,32]) 
pc27 <- as.matrix(sample2_AcD_sleepduration_hba1c[,33]) 
pc28 <- as.matrix(sample2_AcD_sleepduration_hba1c[,34]) 
pc29 <- as.matrix(sample2_AcD_sleepduration_hba1c[,35]) 
pc30 <- as.matrix(sample2_AcD_sleepduration_hba1c[,36]) 

pc31 <- as.matrix(sample2_AcD_sleepduration_hba1c[,37])     
pc32 <- as.matrix(sample2_AcD_sleepduration_hba1c[,38]) 
pc33 <- as.matrix(sample2_AcD_sleepduration_hba1c[,39]) 
pc34 <- as.matrix(sample2_AcD_sleepduration_hba1c[,40]) 
pc35 <- as.matrix(sample2_AcD_sleepduration_hba1c[,41]) 
pc36 <- as.matrix(sample2_AcD_sleepduration_hba1c[,42]) 
pc37 <- as.matrix(sample2_AcD_sleepduration_hba1c[,43]) 
pc38 <- as.matrix(sample2_AcD_sleepduration_hba1c[,44]) 
pc39 <- as.matrix(sample2_AcD_sleepduration_hba1c[,45]) 
pc40 <- as.matrix(sample2_AcD_sleepduration_hba1c[,46]) 


## SNP-Outcome association adjusting for assessment centre and 40 genetic principal components, baseline age, sex, and genotyping chip (fasting time and dilution factor were additionally adjusted for in the case of non-fasting glucose)
YGdata        = data.frame(Y, G1.2,
                           sex, age, as.factor(centre), chip, 
                           pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                           pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                           pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                           pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40)   

FIT3           = summary(lm(YGdata))                          
BetaYG         = FIT3$coef[-1,1]                             
BetaYG         = head(BetaYG, n= -63)                           
seBetaYG       = FIT3$coef[-1,2]                             
seBetaYG         = head(seBetaYG, n= -63)

BetaYG_data<-as.matrix(BetaYG)
colnames(BetaYG_data) <- c("BetaYG")
seBetaYG_data<-as.matrix(seBetaYG)
colnames(seBetaYG_data) <- c("seBetaYG")

data_sample2<-cbind(BetaYG, seBetaYG)
SNP<-as.matrix(row.names(data_sample2))
colnames(SNP) <- c("SNP_ukb")
data_sample2<-cbind(SNP, BetaYG, seBetaYG)
rownames(data_sample2) <- c()
data_sample2<-as.data_sample2.frame(data_sample2)

data_sample2$SNP<-as.factor(str_sub(data_sample2$SNP_ukb, 0, -3))
data_sample2$effect_allele<-as.factor(str_sub(data_sample2$SNP_ukb, -1))
data_sample2$outcome<-"log_hba1c"
data_sample2$samplesize<-n



# Step 3: harmonise the data to run the 2SMR analyses (inverse-variance weighted (IVW),  weighted median (WM),  MR-Egger, MR-Egger_SiMEX)
## Exposure:  The summary statistics of AcD sleep duration and other AcD sleep traits were from the GWAS (Jones S et al, PMID:0952852)
## The discovery GWAS presented the SNP-AcD sleep duration association in SD unit, where 1SD is equal to 0.86 hours. Therefore, We transformed it into per hour change unit.
AcD_sleepduration <- read.csv("AcD_sleepduration.csv")



## Outcome: SNP - HbA1c association, from data_sample2 as obtained in Step 2
### Assign the other allele to the data_sample2 because the individual genotype data has no information of other allele, where the information is from the discovery GWAS summary data (Jones S et al, PMID:0952852)
oa<-select(AcD_sleepduration, SNP, other_allele.exposure)
data_sample2_full<-merge(data_sample2, oa, by = "SNP")

### format outcome data
data_sample2_format<-format_data(data_sample2_full, type = "outcome",
                         snp_col="SNP", beta_col="BetaYG", phenotype_col="phenotype",
                         se_col="seBetaYG",
                         effect_allele_col="effect_allele", other_allele_col="other_allele.exposure")


### Harmonise data, action = 1. 
### Assuming all alleles are presented on the same strand because the SNP-exposure and SNP-outcome associations are both from the UKB). 
### As such, no palindromic SNP is removed
harmonised_AcD_sleepduration_hba1c<-harmonise_data(exposure_dat = AcD_sleepduration, outcome_dat = data_sample2_format, action = 1) 



# Step 4: Run the MR analyses (inverse-variance weighted (IVW),  weighted median (WM),  MR-Egger, MR-Egger_SiMEX)
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


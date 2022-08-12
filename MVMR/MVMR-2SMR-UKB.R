# Overview
# This file contains the R script of the multivariable Mendelian randomization in a two-sample design applying to the UK Biobank (MVMR-2SMR-UKB)
# (inverse-variance weighted (IVW),  weighted median (WM),  MR-Egger)
# Here we take accelerometer-derived (AD) sleep duration + sleep efficiency vs HbA1c as an example (R codes of other sleep traits are similar)
# For HbA1c, it is needed to divide the estimates by  0.15 log mmol/mol to obtain the SD estimates; for non-fasting glucose, it is needed to divide the estimates by 0.18 log mmol/l to obtain the SD unit



# Install related R packages
library(TwoSampleMR)
library(dplyr)
library(stringr)
library(MVMR)

# sample 1 - SNP-Exposures
## Before running the MVMR analyses, the author has extracted the summary statistics (e.g., BetaXG, seBetaXG, effect allele, other allele et al) of all the related sleep traits (e.g., AD sleep druation, efficiency et al) from the Sleep Disorder Knowledge Portal https://sleep.hugeamp.org/.

sample1_all_sleeptraits<-read.csv("sample1_all_sleeptraits.csv")

## Extract the related summary statistics
sample1_ad_sleepduration_efficiency<-select(sample1_all_sleeptraits, 
                                            SNP, effect_allele, other_allele, 
                                            beta_ad_sleepduration, se_ad_sleepduration, eaf_ad_sleepduration,
                                            beta_efficiency, se_efficiency, eaf_efficiency)

# The summary data of the AD sleep traits GWAS provides effect allele frequency (EAF) of each AD sleep trait, which are close to each other. 
# We took EAF of the ad_sleepduration as the reference for EAF of the two sleep traits 
sample1_ad_sleepduration_efficiency <- format_data(sample1_ad_sleepduration_efficiency,
                                                   snp_col = "SNP", 
                                                   beta_col = "beta_ad_sleepduration", se_col = "se_ad_sleepduration", eaf = "eaf_ad_sleepduration",
                                                   effect_allele = "effect_allele", other_allele = "other_allele")





# Sample 2 - SNP-Outcome - obtained from a sub-sample of UKB who did not participate into the AD sleep GWAS
## The summary statistics of SNP-Outcome associations are from the multivariable adjusted linear model accounting for assessment centre and 40 genetic principal components, baseline age, sex, and genotyping chip, fasting time and dilution factor (for glucose only)
sample2_hba1c<-read.csv("C:/Users/hx19384/OneDrive - University of Bristol/Desktop/Sleep non-linear project/data/For MVMR/outcome specific summary statistics ukb/ad_sleepduration_efficiency_hba1c_subukb.csv")



# Combine Sample 1 and Sample 2 for the MVMR analyses 
## Assume all alleles are presented on the forward strand because the SNP-exposure and SNP-outcome associations are both from the UKB
sampleCB<-merge(sample1_ad_sleepduration_efficiency, sample2_hba1c, by = "SNP")

# Format data
F.data <- format_mvmr(BXGs = sampleCB[,c(4,7)],
                      BYG = sampleCB[,11],
                      seBXGs = sampleCB[,c(5,8)],
                      seBYG = sampleCB[,12],
                      RSID = sampleCB[,1])

beta_accel_sleepduration<-F.data$betaX1
beta_efficiency<-F.data$betaX2
se_accel_sleepduration<-F.data$sebetaX1
se_efficiency<-F.data$sebetaX2
BetaYG<-F.data$betaYG
seBetaYG<-F.data$sebetaYG

# IVW
ivw_estimates<-mr_mvivw(mr_mvinput(bx = cbind(beta_accel_sleepduration, beta_efficiency), bxse = cbind(se_accel_sleepduration, se_efficiency),
                                   by = BetaYG, byse = seBetaYG))

# WM
wm_estimates<-mr_mvmedian(mr_mvinput(bx = cbind(beta_accel_sleepduration, beta_efficiency), bxse = cbind(se_accel_sleepduration, se_efficiency),
                                     by = BetaYG, byse = seBetaYG))

# MR-Egger 
egger_estimates<-mr_mvegger(mr_mvinput(bx = cbind(beta_accel_sleepduration, beta_efficiency), bxse = cbind(se_accel_sleepduration, se_efficiency),
                                       by = BetaYG, byse = seBetaYG))





# Calculate the conditional F statistics and Q statistics 
## Before so, it is needed to generate a dataset which contains the individual genotype and phenotype data. 
## The genetic data works for the conditional F statistics when covariance is based on the covariance across the SNPs
individual_geno_pheno_data_ad_sleepduration_efficiency<-fread(file="individual_geno_pheno_data_ad_sleepduration_efficiency.csv",
                                                              fill=T, header = T)

individual_geno_pheno_data_ad_sleepduration_efficiency<-select(individual_geno_pheno_data_ad_sleepduration_efficiency,c(,1:2), contains(sampleCB$SNP)) # The available SNPs as applied in above

## Format data 
F.data <- format_mvmr(BXGs = sampleCB[,c(4,7)],
                      BYG = sampleCB[,11],
                      seBXGs = sampleCB[,c(5,8)],
                      seBYG = sampleCB[,12],
                      RSID = sampleCB[,1])


# conditional F statistics and Q statistics are obtained based on the covariance across the phenotype 
phenocov_mvmr_ad_sleepduration_efficiency<-phenocov_mvmr(cor(individual_geno_pheno_data_ad_sleepduration_efficiency[,1:2]), sampleCB[,c(5,8)])
sres_phenocov <- strength_mvmr(r_input = F.data, gencov = phenocov_mvmr_ad_sleepduration_efficiency)
pres_phenocov<- pleiotropy_mvmr(r_input = F.data, gencov = phenocov_mvmr_ad_sleepduration_efficiency)


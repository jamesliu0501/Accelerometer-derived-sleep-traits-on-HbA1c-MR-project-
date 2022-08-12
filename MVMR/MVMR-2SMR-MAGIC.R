# Overview
# This file contains the R script of the multivariable Mendelian randomization in a two-sample design applying to the MAGIC (MVMR-2SMR-MAGIC)
# (inverse-variance weighted (IVW),  weighted median (WM),  MR-Egger)
# Here we take accelerometer-derived (AD) sleep duration + sleep efficiency vs HbA1c as an example (R codes of other sleep traits are similar)
# For HbA1c, it is needed to divide the estimates by 0.41 % to obtain the SD estimates; for non-fasting glucose, it is needed to divide the estimates by 0.84 mmol/l to obtain the SD unit

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

# Sample 2 - SNP-Outcome  (HbA1c / fasting-glucose (Chen J et al 2021, PMID: 34059833))
# Fasting-glucose was obtained from the same GWAS (Chen J et al 2021, PMID: 34059833)) with "ebi-a-GCST90002232" id in the MRbase
sample2_hba1c <- extract_outcome_data(sample1_ad_sleepduration_efficiency$SNP, c("ebi-a-GCST90002244")) 


# To harmonise 
harmonised <- harmonise_data(sample1_ad_sleepduration_efficiency, sample2_hba1c, action=2)


# Sort two samples for MVMR -- exclude the ambiguous SNPs
snp_sorted<-subset(harmonised, ambiguous == FALSE)
sample1_ad_sleepduration_efficiency_sorted<- subset(sample1_ad_sleepduration_efficiency,  sample1_ad_sleepduration_efficiency$SNP %in% snp_sorted$SNP)
sample2_hba1c_sorted<-select(harmonised, SNP, effect_allele.outcome, other_allele.outcome, eaf.outcome,
                       beta.outcome, se.outcome)
sample_analysis<-merge(sample1_ad_sleepduration_efficiency_sorted, sample2_hba1c_sorted)



# To run the MR analyses
beta_ad_sleepduration<-sample_analysis$beta_ad_sleepduration
beta_efficiency<-sample_analysis$beta_efficiency
se_ad_sleepduration<-sample_analysis$se_ad_sleepduration
se_efficiency<-sample_analysis$se_efficiency
BetaYG<-sample_analysis$beta.outcome
seBetaYG<-sample_analysis$se.outcome

## IVW
ivw_estimates<-mr_mvivw(mr_mvinput(bx = cbind(beta_ad_sleepduration, beta_efficiency), bxse = cbind(se_ad_sleepduration, se_efficiency),
                                   by = BetaYG, byse = seBetaYG))

## WM 
wm_estimates<-mr_mvmedian(mr_mvinput(bx = cbind(beta_ad_sleepduration, beta_efficiency), bxse = cbind(se_ad_sleepduration, se_efficiency),
                                     by = BetaYG, byse = seBetaYG))


## MR-Egger 
egger_estimates<-mr_mvegger(mr_mvinput(bx = cbind(beta_ad_sleepduration, beta_efficiency), bxse = cbind(se_ad_sleepduration, se_efficiency),
                                       by = BetaYG, byse = seBetaYG))





# Calculate the conditional F statistics and Q statistics 
## Before so, it is needed to generate a dataset which contains the individual genotype and phenotype data. 
## The genetic data works for the conditional F statistics when covariance is based on the covariance across the SNPs
individual_geno_pheno_data_ad_sleepduration_efficiency<-fread(file="individual_geno_pheno_data_ad_sleepduration_efficiency.csv",
                                              fill=T, header = T)

individual_geno_pheno_data_ad_sleepduration_efficiency<-select(individual_geno_pheno_data_ad_sleepduration_efficiency,c(,1:2), contains(sample_analysis$SNP)) # The available SNPs as applied in above

## Format data 
F.data <- format_mvmr(BXGs = sample_analysis[,c(4,7)],
                      BYG = sample_analysis[,13],
                      seBXGs = sample_analysis[,c(5,8)],
                      seBYG = sample_analysis[,14],
                      RSID = sample_analysis[,1])

# conditional F statistics and Q statistics are obtained based on the covariance across the phenotype 
phenocov_mvmr_ad_sleepduration_efficiency<-phenocov_mvmr(cor(individual_geno_pheno_data_ad_sleepduration_efficiency[,1:2]), sample_analysis[,c(5,8)])
sres_phenocov <- strength_mvmr(r_input = F.data, gencov = phenocov_mvmr_ad_sleepduration_efficiency)
pres_phenocov<- pleiotropy_mvmr(r_input = F.data, gencov = phenocov_mvmr_ad_sleepduration_efficiency)


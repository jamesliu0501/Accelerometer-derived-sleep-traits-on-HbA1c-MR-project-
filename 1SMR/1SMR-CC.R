# Overview
# This file contains the R script of the collider-correction in a one-sample design applying to the UK Biobank (1SMR-CC)
# (inverse-variance weighted (IVW),  and least absolute deviation regression (LADreg) ,  MR-Egger)
# Here we take accelerometer-derived (AcD) sleep duration vs HbA1c as an example (R codes of other sleep traits are similar)
# Sensitivity analyses excluding participants with diabetes (HbA1c >= 48 mmol/mol or diabetic status defined by the Eastwood algorithm )
# For HbA1c, it is needed to divide the estimates by 0.14 log mmol/mol to obtain the SD estimates; for non-fasting glucose, it is needed to divide the estimates by 0.16 log mmol/l to obtain the SD unit


# Install related R packages
library(foreign)
library(simex)
library(L1pack)
library(mr.raps)
library(coda)
library(xtable)
library(stringr)
library(data.table)
library(dplyr)


# Step 1: Generate the summary statistics for collider-correction estimates
## The following dataset is generated from the UKB with genotype data, exposure (AcD sleep duration), outcome (HbA1c), age, sex, chip, assessment centre, and 40 genetic principal components, (dilution factor and fasting time were also included in the case of non-fasting glucose)
data1 <- read.csv("data1.csv") 

n <- nrow(data1)
size <- ncol(data1[,-(1:47)])
G1.2 <- as.matrix(data1[,-(1:47)])
X.2 <- as.matrix(data1[,2])
Y <- as.matrix(data1[,3])

#covariates
sex<-as.matrix(data1[,5])
age<-as.matrix(data1[,4])
centre<-as.matrix(data1[,6])
chip<-as.matrix(data1[,7])

pc1 <- as.matrix(data1[,8])     
pc2 <- as.matrix(data1[,9]) 
pc3 <- as.matrix(data1[,10]) 
pc4 <- as.matrix(data1[,11]) 
pc5 <- as.matrix(data1[,12]) 
pc6 <- as.matrix(data1[,13]) 
pc7 <- as.matrix(data1[,14]) 
pc8 <- as.matrix(data1[,15]) 
pc9 <- as.matrix(data1[,16]) 
pc10 <- as.matrix(data1[,17]) 

pc11 <- as.matrix(data1[,18])     
pc12 <- as.matrix(data1[,19]) 
pc13 <- as.matrix(data1[,20]) 
pc14 <- as.matrix(data1[,21]) 
pc15 <- as.matrix(data1[,22]) 
pc16 <- as.matrix(data1[,23]) 
pc17 <- as.matrix(data1[,24]) 
pc18 <- as.matrix(data1[,25]) 
pc19 <- as.matrix(data1[,26]) 
pc20 <- as.matrix(data1[,27]) 

pc21 <- as.matrix(data1[,28])     
pc22 <- as.matrix(data1[,29]) 
pc23 <- as.matrix(data1[,30]) 
pc24 <- as.matrix(data1[,31]) 
pc25 <- as.matrix(data1[,32]) 
pc26 <- as.matrix(data1[,33]) 
pc27 <- as.matrix(data1[,34]) 
pc28 <- as.matrix(data1[,35]) 
pc29 <- as.matrix(data1[,36]) 
pc30 <- as.matrix(data1[,37]) 

pc31 <- as.matrix(data1[,38])     
pc32 <- as.matrix(data1[,39]) 
pc33 <- as.matrix(data1[,40]) 
pc34 <- as.matrix(data1[,41]) 
pc35 <- as.matrix(data1[,42]) 
pc36 <- as.matrix(data1[,43]) 
pc37 <- as.matrix(data1[,44]) 
pc38 <- as.matrix(data1[,45]) 
pc39 <- as.matrix(data1[,46]) 
pc40 <- as.matrix(data1[,47]) 


## The summary statistics were obtained from the linear regression adjusted for age, sex, chip, assessment centre, and 40 genetic principal components. Dilution factor and fasting time were additionally adjusted for in the case of non-fasting glucose.
## X~G
XGdata2        = data.frame(X.2, G1.2,
                            sex, age, as.factor(centre), chip, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                            pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40)   

FIT2           = summary(lm(XGdata2))                    
BetaXG         = FIT2$coef[-1,1]                               
BetaXG         = head(BetaXG, n= -63)                          
seBetaXG       = FIT2$coef[-1,2]                                
seBetaXG         = head(seBetaXG, n= -63)                           
Fbar           = mean((BetaXG^2)/(seBetaXG^2))                 

## Y~G 
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

## Y~G+X 
YXGdata        = data.frame(Y, X.2, G1.2,
                            sex, age, as.factor(centre), chip, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25, pc26, pc27, pc28, pc29, pc30,
                            pc31, pc32, pc33, pc34, pc35, pc36, pc37, pc38, pc39, pc40) 


FIT4             = summary(lm(YXGdata))                       
alphahatstar     = FIT4$coef[-c(1,2),1]                      
alphahatstar         = head(alphahatstar, n= -63)                           
se.alphahatstar  = FIT4$coef[-c(1,2),2]
se.alphahatstar         = head(se.alphahatstar, n= -63)                           
betastar         = FIT4$coef[2,1] 
se.betastar  =  FIT4$coef[2,2] 


## summary statistics
## X- G
BetaXG_data<-as.matrix(BetaXG)
colnames(BetaXG_data) <- c("BetaXG")
seBetaXG_data<-as.matrix(seBetaXG)
colnames(seBetaXG_data) <- c("seBetaXG")

## Y - G
BetaYG_data<-as.matrix(BetaYG)
colnames(BetaYG_data) <- c("BetaYG")
seBetaYG_data<-as.matrix(seBetaYG)
colnames(seBetaYG_data) <- c("seBetaYG")

## Y - X - G
alphahatstar_data<-as.matrix(alphahatstar)
colnames(alphahatstar_data) <- c("alphahatstar")
se.alphahatstar_data<-as.matrix(se.alphahatstar)
colnames(se.alphahatstar_data) <- c("se.alphahatstar")

betastar_data<-as.matrix(betastar)
colnames(betastar_data) <- c("betastar")
se.betastar_data<-as.matrix(se.betastar)
colnames(se.betastar_data) <- c("se.betastar")

data_summary_cc<-cbind(BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)

SNP<-as.matrix(row.names(data_summary_cc))
colnames(SNP) <- c("SNP_ukb")

data_summary_cc<-cbind(SNP, BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)
rownames(data_summary_cc) <- c()


# Step 2 
set.seed(888)
#import "ExactQ.R" function file (see the uploaded file)
source("ExactQ.R") 
## For sleep efficiency, predicted by 5 SNPs, the "ExactQ_parametric_bootstrap. R" file should be imported.  
## For the sleep efficiency analysis only, standard errors for the collider-corrected IVW approach we obtained using a parametric bootstrap. In this case the non-parametric bootstrap (used without issue in other analyses) provided unstable estimates due to the small number of SNPs used.

# function of generating the IGx2 
Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

##  The SNP applied to collider-correction is required a positive BetaXG value obtained in the UKB (i.e., align the effect allele corresponding to the effet allele with a positive BetaXG value in the AcD sleep trait GWAS)
data_summary_cc <- data_summary_cc[which(data_summary_cc$BetaXG>0),] 

BetaXG<-data_summary_cc$BetaXG 
seBetaXG<-data_summary_cc$seBetaXG
BetaYG<-data_summary_cc$BetaYG 
seBetaYG<-data_summary_cc$seBetaYG
alphahatstar<-data_summary_cc$alphahatstar 
se.alphahatstar<-data_summary_cc$se.alphahatstar  
betastar<-data_summary_cc$betastar 
se.betastar<-data_summary_cc$se.betastar

## F-statistic
F = BetaXG^2/seBetaXG^2
plot(BetaXG,F)

## Implement IVW using Modified weights 
set.seed(888)
Results = weightedIVW(BetaXG,alphahatstar,seBetaXG,se.alphahatstar,tol=0.00001)
IVW_mw  = Results$RESULTS[5,1]
Qexact  = Results$QStats[4,1]
Qp      = Results$QStats[4,2]

betaIVW = betastar[1] + IVW_mw
names(Results)


## MR-Egger
set.seed(888)
IsqGX         = Isq(BetaXG,seBetaXG)                                                   
betahat       = summary(lm(alphahatstar~ BetaXG,weights=1/se.alphahatstar^2))$coef    
MREgger       = betastar[1] + betahat[2,1]  
## Collider correction + SiMEX 
Fit           = lm(alphahatstar~BetaXG,x=TRUE,y=TRUE,weights=1/se.alphahatstar^2)        
mod.sim2      = simex(Fit,B=500,measurement.error=seBetaXG,                              
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE")
bSIMEX        = summary(mod.sim2)$coef$jackknife   
MREggersimex  = betastar[1] + bSIMEX[2,1]
EggerSE       = bSIMEX[2,2]    
Egger_intercept = bSIMEX[1,4]


## LAD regression
set.seed(888)
betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]   
LAD           = betastar[1] + betahat  

## Collider correction + SiMEX 
set.seed(888)
Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)         
mod.sim3      = simex(Fit,B=500,measurement.error=seBetaXG,       
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE",jackknife.estimation = FALSE)
bSIMEX        = mod.sim3$coef               
LADsimex      = betastar[1] + bSIMEX

## obtain bootstrap SE for LAD regression 
Ests = NULL
for(i in 1:100){
  L     = length(BetaXG)  
  d     = sample(L,L,replace=TRUE)  
  data3 = data_summary_cc[d,]   
  
  BetaXG          = data3$BetaXG
  seBetaXG        = data3$seBetaXG
  alphahatstar    = data3$alphahatstar
  se.alphahatstar = data3$se.alphahatstar
  
  # LAD regression
  betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]      
  LAD           = betastar[1] + betahat                              
  
  Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)        
  mod.sim2      = simex(Fit,B=200,measurement.error=seBetaXG,
                        SIMEXvariable="BetaXG",fitting.method="quad",
                        asymptotic="FALSE",jackknife.estimation = FALSE)
  Ests[i]       = mod.sim2$coef
  print(i)               
}
seLAD = sd(Ests) 

## sort the estimates
BetaXG          = data_summary_cc$BetaXG 
seBetaXG        = data_summary_cc$seBetaXG
alphahatstar    = data_summary_cc$alphahatstar
se.alphahatstar = data_summary_cc$se.alphahatstar

Fbar=mean(F)
Stats               = data.frame(Fbar,IsqGX,Qexact,Qp)
Estimates           = c(betaIVW,MREggersimex,LADsimex)

ColliderCorrections = Estimates-betastar[1]

SEs                 = sqrt(c(Results$RESULTS[5,2],EggerSE,seLAD)^2 + (se.betastar[1])^2)

LCIs                = Estimates - 1.96 * SEs

UCIs                = Estimates + 1.96 * SEs

pval                = 2*(1-pnorm(abs(Estimates/SEs)))

Final = data.frame(Estimates,SEs,pval,
                   row.names=c("IVW","MR-Egger","LADreg"))

fsNigx<-c(Fbar,IsqGX,NA)
intercept<-c(NA, Egger_intercept, NA)

Sorted_estimate = data.frame(Estimates, LCIs, UCIs, pval,fsNigx, intercept,
                                           row.names=c("IVW","MR-Egger","LADreg"))







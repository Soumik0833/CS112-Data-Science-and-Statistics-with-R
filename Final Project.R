##Code adopted and modifed from https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/25680/A9MKGP&version=2.0 


## THIS CODE IMPLEMENTS GENETIC MATCHING AND CREATES A MATCHED DATA SET
## CALLED 'data.mahal' (3664 OBSERVATIONS) WHICH IS USED FOR MY ANALYSIS.
## THE CODE ALSO CREATES A BALANCE TABLE SHOWING MEAN VALUES FOR OBESE AND NON-OBESE GROUPS
## BEFORE AND AFTER MATCHING.

set.seed(123)
############################################
### INSTALL AND CALL NECESSARY LIBRARIES ###
############################################
#install.packages("sandwich")
#install.packages("foreign")
#install.packages("stargazer")
#install.packages("stats")
#install.packages("lmtest")
#install.packages("weights")
#install.packages("cem")
#install.packages("MatchIt")
#install.packages("optmatch")
#install.packages("plyr")


library(sandwich)
library(foreign)
library(stargazer)
library(stats)
library(lmtest)
library(weights)
library(cem)
library(MatchIt)
library(optmatch)
library(plyr)
library(Matching)
source("AutoCluster4.R")


###################
### IMPORT DATA ###
###################

nlsy_norandom<-read.dta("maintmp.dta") #ORIGINAL DATA

#data with random samples taken to improve the speed of genetic matching
nlsy <- nlsy_norandom[sample(1:nrow(nlsy_norandom), 10000,replace=FALSE),] #


########################
### DATA PREPARATION ###
########################
#create list of control variables
controls<-c("d_sexf","childany","childf","d_race_b","d_race_o","d_marrnever","d_marroth","age","age2","d_educ9_12",
            "d_educ13_up","d_AFQT_0_25","d_AFQT_25_50","d_AFQT_50_75","d_AFQT_75_100",
            "d_urban_res","d_tenure0_1","d_tenure1_3","d_tenure3_6","d_tenure6_up",
            "d_emp0_9","d_emp10_24","d_emp25_49","d_emp50_999","d_emp1000_up",
            "d_year1989","d_year1990","d_year1992","d_year1993","d_year1994","d_year1996","d_year1998",
            "d_year2000","d_year2002","d_ind_ag","d_ind_for","d_ind_mining","d_ind_const","d_ind_mfrg",
            "d_ind_transp","d_ind_wtrade","d_ind_rtrade","d_ind_finance","d_ind_bus_svc","d_ind_pers_svc",
            "d_ind_entert","d_ind_prof_svc","d_ind_pub_ad","d_occ_mgmt","d_occ_tech","d_occ_admin","d_occ_svc",
            "d_occ_farming","d_occ_prodxn","d_occ_operators","d_occ_military")

#subset data to include only variables used in analyses
varlist<-c("person_id","ln_smwage","CPS_hourly_rec","d_hinsEMP","d_hinsEMPOTH","d_hinsINDCOV","d_hinsPUBLIC","d_hinsNONE",
           "d_hinsUNKNOWN","overwt","d_obese","d_obese1","d_obese2","d_overwt","overwtinsEMP",
           "d_obese1insEMP","d_obese2insEMP","d_obesinsEMP", controls)
nlsy<-nlsy[,c(varlist,"sample_wt","srvy_yr","birth_year","AFQTrevised",
              "educ","tenure","BMI","NumberEmp")]

#######################################
### IMPLEMENT PROPRIETARY FUNCTIONS ###
#######################################

## Calculate cluster-robust standard errors (CRSEs)
# Implement a linear model CRSE function - via Mahmood Arai @ Univ. Stockholm
clx <- function(fm, cluster){
  M <- length(unique(cluster))
  N <- length(cluster)
  dfc <- (M/(M-1))*((N-1)/(N-fm$rank))
  u <- apply(estfun(fm),2,function(x) tapply(x, cluster, sum))
  vcovCL <- dfc*sandwich(fm, meat=crossprod(u)/N)
  coeftest(fm, vcovCL) }


## Implement a weighted Welch 2 sample t-test function
welch <- function(nx,ny,mx,my,vx,vy){
  # mx = weighted mean of x, nx = number of observations in x vector, vx = weighted variance of x
  # http://en.wikipedia.org/wiki/Student's_t-test#Equal_or_Unequal_sample_sizes.2C_unequal_variances
  
  t.stat <- (mx - my)/(sqrt(vx/nx+vy/ny))
  df <- (vx/nx + vy/ny)^2/
    ((vx/nx)^2/(nx-1))+((vy/ny)^2/(ny-1)) ##NB may need to round this for it to work
  dt(t.stat, df=df)
}

# Implement a function to calculate weighted SD 
weighted.var <- function(x, w) {
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2)
}

################################################################
## DATA PRE-PROCESSING: GENETIC MATCHING ###
################################################################

#to initiate parallel computation for genetic matching to run faster.
# if your computer is quad-cored, please change the number to 2 instead of 4
cl <- NCPUS(4) 

mahal <- matchit(formula = d_obese~d_sexf+childany+childf+age+d_urban_res+
                   srvy_yr+AFQTrevised+educ+tenure+NumberEmp+
                   d_race_b+d_race_o+d_marrnever+d_marroth+
                   d_ind_ag+d_ind_for+d_ind_mining+d_ind_const+d_ind_mfrg+
                   d_ind_transp+d_ind_wtrade+d_ind_rtrade+d_ind_finance+
                   d_ind_bus_svc+d_ind_pers_svc+d_ind_entert+d_ind_prof_svc+
                   d_occ_mgmt+d_occ_tech+d_occ_admin+d_occ_svc+
                   d_occ_farming+d_occ_prodxn+d_occ_operators,
                 data = nlsy,
                 method = "genetic",
                 distance = "logit", max.generations=25,
                 wait.generations=3, hard.generation.limit=TRUE, cluster = cl, nboots = 500, 
                 discard="control")

data.mahal<-match.data(mahal)

mahal.imb<-imbalance(data.mahal$d_obese,data=data.mahal,drop=c("d_obese","CPS_hourly_rec","distance","weights",
                                                               "sample_wt","person_id","ln_smwage","CPS_hourly_rec",
                                                               "d_educ9_12","d_educ13_up",
                                                               "d_AFQT_0_25","d_AFQT_25_50","d_AFQT_50_75","d_AFQT_75_100",
                                                               "d_tenure0_1","d_tenure1_3","d_tenure3_6","d_tenure6_up",
                                                               "d_emp0_9","d_emp10_24","d_emp24_49","d_emp50_999","d_emp1000_up",
                                                               "d_hinsEMP","d_hinsEMPOTH","d_hinsINDCOV",
                                                               "d_hinsPUBLIC","d_hinsNONE","d_hinsUNKNOWN","overwt",
                                                               "d_obese","d_obese1","d_obese2","d_overwt",
                                                               "overwtinsEMP","d_obese1insEMP","d_obese2insEMP",
                                                               "d_obesinsEMP","BMI","age2"))

mahal.imb

################################################
## DATA PRE-PROCESSING: CREATE BALANCE TABLE ###
################################################
# ABBREVIATIONS KEY:
# oo - original obese; on - original non-obese; mo - matched obese; mn - matched non-obese


### Create balance tables pre- and post-matching
## Create vector of weighted mean value of variables for obese and non-obese subsets in pre- & post-match data
oo <- apply(nlsy[which (nlsy$d_obese==1),],2,weighted.mean,w=nlsy$sample_wt[nlsy$d_obese==1])
on <- apply(nlsy[which (nlsy$d_obese==0),],2,weighted.mean,w=nlsy$sample_wt[nlsy$d_obese==0])
mo <- apply(data.mahal[which (data.mahal$d_obese==1),],2,weighted.mean,w=data.mahal$sample_wt[data.mahal$d_obese==1])
mn <- apply(data.mahal[which (data.mahal$d_obese==0),],2,weighted.mean,w=data.mahal$sample_wt[data.mahal$d_obese==0])

## Calculate number of observations in each grouop
oo.n <- length(nlsy$person_id[nlsy$d_obese==1])
on.n <- length(nlsy$person_id[nlsy$d_obese==0])
mo.n <- length(data.mahal$person_id[data.mahal$d_obese==1])
mn.n <- length(data.mahal$person_id[data.mahal$d_obese==0])

## Calculate weighted variances for vectors of means
oo.var <- apply(nlsy[which (nlsy$d_obese==1),],2,weighted.var,w=nlsy$sample_wt[nlsy$d_obese==1])
on.var <- apply(nlsy[which (nlsy$d_obese==0),],2,weighted.var,w=nlsy$sample_wt[nlsy$d_obese==0])
mo.var <- apply(data.mahal[which (data.mahal$d_obese==1),],2,weighted.var,w=data.mahal$sample_wt[data.mahal$d_obese==1])
mn.var <- apply(data.mahal[which (data.mahal$d_obese==0),],2,weighted.var,w=data.mahal$sample_wt[data.mahal$d_obese==0])

## Trim last 2 columns of post-match means and variances (weights and distance)
mo <- mo[1:82]
mn <- mn[1:82]
mo.var <- mo.var[1:82]
mn.var <- mn.var[1:82]

## Calculate p-values for weighted Welch t-test of difference in means between obese and non-obese
o.pvalues <- welch(on.n,oo.n,on,oo,on.var,oo.var)
m.pvalues <- welch(mn.n,mo.n,mn,mo,mn.var,mo.var)

## Compile into a table
match.table <- matrix(nrow=length(colnames(nlsy)), ncol=7)
match.table[,1] <- round(t(on),2)
match.table[,2] <- round(t(oo),2)
match.table[,3] <- round(t(o.pvalues),3)
match.table[,4] <- round(t(mn),2)
match.table[,5] <- round(t(mo),2)
match.table[,6] <- round(t(m.pvalues),3)
match.table[,7] <- ifelse(match.table[,3]<0.05 & match.table[,6]>0.05, "Yes", "No")
# Drop variables not of interest
match.table <- match.table[c(3:4, 8, 11:14, 18, 19:26, 28:74,81),]
# Label table
var.mt <- c("Hourly wage*", "Employer coverage in own name", "Uninsured", "Obese (BMI>30)", 
            "Mildly obese (BMI>30 and BMI<35)", "Morbidly obese (BMI>35)", "Overweight",
            "Obese * employer coverage (own)", "Female", "Any children in household", 
            "Female with children in household","Race - Black", "Race - Other", "Never married", 
            "Formerly married","Age", "Education: 9-12", "Education: 13 and over", "AFQT: 0-25", 
            "AFQT: 25-50", "AFQT: 50-75", "AFQT: 75-100", "Urban residence", "Job tenure: 0-1 years", 
            "Job tenure: 1-3 years", "Job tenure: 3-6 years", "Job tenure: 6+ years","Employer size: 0-9",
            "Employer size: 10-24","Employer size: 25-49", "Employer size: 50-999","Employer size: 1000+",
            "Survey year: 1989", "Survey year: 1990", "Survey year: 1992", "Survey year: 1993",
            "Survey year: 1994", "Survey year: 1996", "Survey year: 1998", "Survey year: 2000", 
            "Survey year: 2002", "Industry: Agriculture", "Industry: Forestry", "Industry: Mining",
            "Industry: Construction", "Industry: Manufacturing", "Industry: Transport", 
            "Industry: Wholesale trade", "Industry: Retail trade", "Industry: Finance", 
            "Industry: Business services", "Industry: Personal services", "Industry: Entertainment",
            "Industry: Professional services", "Industry: Public administration", 
            "Occupation: Management", "Occupation: Technical", "Occupation: Administrative", 
            "Occupation: Service", "Occupation: Farming", "Occupation: Production", 
            "Occupation: Operators", "Occupation: Military", "BMI")
rownames(match.table) <- var.mt
colnames(match.table) <- c("Pre-match non-obese", "Pre-match obese", "Pre-match p-value","Post-match non-obese", 
                           "Post-match obese", "Matched p-value", "Improved balance with matching")
View(match.table)      # Inspect table
count(match.table[,7]) # Count number of variables where balance improved (NB - needs 'plyr' package)
## Matching significantly improved balance on 26 out of 64 variables


### CREATE A TABLE THAT SUMMARIZES EFFECT OF BALANCING
improve.table <- matrix(nrow=nrow(match.table), ncol=3)
rownames(improve.table) <- rownames(match.table)
colnames(improve.table) <- c("Balanced before matching", "Balanced after matching", 
                             "Balance improved by matching")
improve.table[,1] <- ifelse(match.table[,3]>0.05, "Yes", "No")
improve.table[,2] <- ifelse(match.table[,6]>0.05, "Yes", "No")
improve.table[,3] <- match.table[,7]
improve.table <- improve.table[c(1:3,9:63),] #Drop weight-related variables
View(improve.table)
count(improve.table)
#26 out of 32 variables that vere not balanced were balanced and their balance improved
## 7 variables remained unbalanced after matching
## 26 variables became balanced after matching
## 1 variable became unbalanced after matching
## 25 variables remained balanced after matching


## Clean up
rm(on.n)
rm(oo.n)
rm(on)
rm(oo)
rm(on.var)
rm(oo.var)
rm(mn.n)
rm(mo.n)
rm(mn)
rm(mo)
rm(mn.var)
rm(mo.var)
rm(var.mt)

#### Original results for comparison - Table 3 Model 1 ####
## Create model 1: Main study sample
t3.m1.orig<-lm(CPS_hourly_rec~d_obese+d_hinsEMP+d_obesinsEMP+d_sexf+childany+childf+d_race_b+d_race_o+
            d_marrnever+d_marroth+age+age2+d_educ9_12+d_educ13_up+d_AFQT_0_25+d_AFQT_25_50+d_AFQT_50_75+
            d_AFQT_75_100+d_urban_res+d_tenure0_1+d_tenure1_3+d_tenure3_6+d_tenure6_up+d_emp0_9+d_emp10_24+
            d_emp25_49+d_emp50_999+d_emp1000_up+d_year1989+d_year1990+d_year1992+d_year1993+d_year1994+
            d_year1996+d_year1998+d_year2000+d_year2002+d_ind_ag+d_ind_for+d_ind_mining+d_ind_const+
            d_ind_mfrg+d_ind_transp+d_ind_wtrade+d_ind_rtrade+d_ind_finance+d_ind_bus_svc+d_ind_pers_svc+
            d_ind_entert+d_ind_prof_svc+d_ind_pub_ad+d_occ_mgmt+d_occ_tech+d_occ_admin+d_occ_svc+d_occ_farming+
            d_occ_prodxn+d_occ_operators+d_occ_military, weight=sample_wt, data=nlsy_norandom)

# Calculate cluster-robust SEs for model 1
m1.orig<-clx(t3.m1.orig,nlsy_norandom$person_id)

#### Recreate main results with matched sample ####

## Create model 1: Main study sample
t3.m1<-lm(CPS_hourly_rec~d_obese+d_hinsEMP+d_obesinsEMP+d_sexf+childany+childf+d_race_b+d_race_o+
            d_marrnever+d_marroth+age+age2+d_educ9_12+d_educ13_up+d_AFQT_0_25+d_AFQT_25_50+d_AFQT_50_75+
            d_AFQT_75_100+d_urban_res+d_tenure0_1+d_tenure1_3+d_tenure3_6+d_tenure6_up+d_emp0_9+d_emp10_24+
            d_emp25_49+d_emp50_999+d_emp1000_up+d_year1989+d_year1990+d_year1992+d_year1993+d_year1994+
            d_year1996+d_year1998+d_year2000+d_year2002+d_ind_ag+d_ind_for+d_ind_mining+d_ind_const+
            d_ind_mfrg+d_ind_transp+d_ind_wtrade+d_ind_rtrade+d_ind_finance+d_ind_bus_svc+d_ind_pers_svc+
            d_ind_entert+d_ind_prof_svc+d_ind_pub_ad+d_occ_mgmt+d_occ_tech+d_occ_admin+d_occ_svc+d_occ_farming+
            d_occ_prodxn+d_occ_operators+d_occ_military, weight=sample_wt, data=data.mahal)

# Calculate cluster-robust SEs for model 1
m1<-clx(t3.m1,data.mahal$person_id)

## Compile into table 3 and output aS LaTeX code 
## which can then be converted on overleaf.com to produce well structured PDF.

stargazer(m1.orig, m1, title="Estimates of the obesity wage offset for health insurance, Original and Matched Samples", 
          align=TRUE,  no.space=TRUE,column.labels=c("Original", "Matched"),
          covariate.labels=c("Obese","Insured","Obese x Insured"),
          keep=c("d_obese", "d_hinsEMP","d_obesinsEMP","Constant"))

_______________________________________________________________________________________________________________________________
#genmatch output

NOTE: HARD MAXIMUM GENERATION LIMIT HIT

Solution Lexical Fitness Value:
0.000000e+00  1.268498e-06  2.095827e-06  2.095827e-06  3.608674e-06  3.608674e-06  2.570082e-05  2.570082e-05  2.912801e-05  8.991836e-05  6.018550e-04  2.240167e-03  2.240167e-03  1.424662e-02  1.424662e-02  1.424662e-02  1.424662e-02  2.400000e-02  3.091702e-02  3.136872e-02  3.136872e-02  3.729570e-02  3.729570e-02  6.800000e-02  7.800000e-02  9.459862e-02  9.459862e-02  1.796701e-01  1.796701e-01  2.513309e-01  2.513309e-01  2.580000e-01  2.640000e-01  3.000000e-01  3.173106e-01  3.173106e-01  3.173106e-01  3.173106e-01  3.173106e-01  3.173106e-01  3.173106e-01  3.173106e-01  3.457852e-01  3.457852e-01  3.650408e-01  3.650408e-01  5.127202e-01  5.127202e-01  5.398699e-01  5.840000e-01  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  

Parameters at the Solution:

 X[ 1] :	9.088441e+02
 X[ 2] :	4.141232e+02
 X[ 3] :	5.307889e+02
 X[ 4] :	5.055697e+01
 X[ 5] :	6.022365e+02
 X[ 6] :	6.672560e+02
 X[ 7] :	5.354206e+02
 X[ 8] :	7.866402e+02
 X[ 9] :	2.821236e+01
 X[10] :	6.919897e+01
 X[11] :	4.891194e+01
 X[12] :	1.326077e+02
 X[13] :	8.213863e+02
 X[14] :	5.633170e+02
 X[15] :	3.987195e+02
 X[16] :	5.313649e+02
 X[17] :	8.846210e+02
 X[18] :	5.301348e+02
 X[19] :	5.703647e+02
 X[20] :	9.850266e+02
 X[21] :	6.499526e+02
 X[22] :	6.822988e+02
 X[23] :	4.744775e+01
 X[24] :	9.162633e+00
 X[25] :	4.980530e+02
 X[26] :	1.211878e+02
 X[27] :	5.272579e+02
 X[28] :	7.274917e+02
 X[29] :	1.804035e+02
 X[30] :	3.909940e+02
 X[31] :	7.472372e+02
 X[32] :	1.923706e+02
 X[33] :	6.702722e+01
 X[34] :	9.138133e+02
 X[35] :	3.199009e+02

Solution Found Generation 25
Number of Generations Run 25

Thu Dec 20 19:37:35 2018
Total run time : 1 hours 2 minutes and 47 seconds

___________________________________________________________________________________________________________________
 #AutoCluster4.R (Helps with parallel computing for GenMatch so that results can be obtained faster).
 
 # Jasjeet S. Sekhon
# UC Berkeley
# nov 13, 2013
#
# version 4.1, updated for parallel (instead of the older snow package)
# Version 3, fixes a warning message on OS X
# Version 2, runs also on Windows but requires special version of snow, but no cygwin
#

library(rgenoud)
library(parallel)

NCPUS <- function(nchips=FALSE, ...)
{
  nchips <- as.integer(nchips)
  if(!is.integer(nchips))
    stop("'nchips' must be an integer.\n")
  
  if(nchips < 0)
    stop("'nchips' must be a positive integer if it is specified.\n")
  
  
  #should we figure out how many cpus to use?
  if(nchips==FALSE)
  {
    if(.Platform$OS.type=="unix")
    {
      sys <- as.list(Sys.info())$sysname
      
      if(is.null(sys))
      {
        stop("as.list(Sys.info())$sysname returned a NULL")            
      }
      
      if(sys=="Linux" | sys=="linux")
      {
        foo <- paste("cat /proc/cpuinfo | grep processor | wc | awk '{print $1}'")
        nchips <- as.integer(system(foo, intern=TRUE))
        
      } else if (sys=="Darwin" | sys=="darwin") #this should also work for any BSD system?
      {
        foo <- paste("sysctl hw.ncpu | head -1 | awk '{print $2}'")
        nchips <- as.integer(system(foo, intern=TRUE))            
      } else {
        stop("unix system not found. Please explicitly specify the number of chips you want to use.")
      }
      
    } else {
      stop("I think you are using Windows. For this OS please explicitly specify the number of chips you want to use.") 
    }
  }
  
  cat("...using",nchips,"CPUs on this computer.\n")
  
  if(.Platform$OS.type != "windows")
  {
    cl <- makeCluster(rep("localhost", nchips), type="PSOCK")
  } else {
    warning("it may not work on Windows")
    cl <- makeCluster(rep("localhost", nchips), type="SOCK")
  }
  return(cl)
}
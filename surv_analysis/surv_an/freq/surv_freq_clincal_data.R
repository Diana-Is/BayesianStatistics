setwd("C:/Users/Diana/Desktop/Bayes_project/")
library(survival)
library(survminer)


data <- read.csv("data_cleaned_add.csv", header=T)
data$status<-data$death_days
data$status[which(!is.na(data$status))]<-1
data$status[which(is.na(data$status))]<-0

data$pkyr<-(data$pkyr-mean(data$pkyr))/sd(data$pkyr)
data$age<-(data$age-mean(data$age))/sd(data$age)
data$progression_num<-(data$progression_num-mean(data$progression_num))/sd(data$progression_num)
data$lesionsize<-(data$lesionsize-mean(data$lesionsize))/sd(data$lesionsize)
data$BMI<-(data$BMI-mean(data$BMI))/sd(data$BMI)
data$lung_dang_work<-(data$lung_dang_work-mean(data$lung_dang_work))/sd(data$lung_dang_work)
data$alcohol_per_year<-(data$alcohol_per_year-mean(data$alcohol_per_year))/sd(data$alcohol_per_year)
data$yrs_no_smoke<-(data$yrs_no_smoke-mean(data$yrs_no_smoke))/sd(data$yrs_no_smoke)


#-----------------------------------------------------
#CREATING DUMMIES
#-----------------------------------------------------

names(data)[names(data) == "rndgroup"] <- "rnd_group"
names(data)[names(data) == "marital"] <- "marital_status"
names(data)[names(data) == "de_stag_7thed"] <- "stage"
names(data)[names(data) == "smokelive"] <- "live_with_smoker"
names(data)[names(data) == "smokework"] <- "work_with_smoker"
names(data)[names(data) == "pkyr"] <- "package_year"
names(data)[names(data) == "pipe"] <- "smoke_pipe"
names(data)[names(data) == "cigar"] <- "smoke_sigar"
names(data)[names(data) == "diagasbe"] <- "Asbestosis"
names(data)[names(data) == "diagadas"] <- "Astma_adult"
names(data)[names(data) == "diagchas"] <- "Astma_childhood"
names(data)[names(data) == "diagbron"] <- "Bronchitis" 
names(data)[names(data) == "diagchro"] <- "Bronchitis_chronical" 
names(data)[names(data) == "diagcopd"] <- "COPD"   
names(data)[names(data) == "diagdiab"] <- "Diabetes"   
names(data)[names(data) == "diagemph"] <- "Emphysema"   
names(data)[names(data) == "diaghear"] <- "Heart_disease"   
names(data)[names(data) == "diagpneu"] <- "Pneumonia"   
names(data)[names(data) == "diagsarc"] <- "Sarcoidosis"   
names(data)[names(data) == "diagtube"] <- "Tuberculosis"   
names(data)[names(data) == "diaghype"] <- "Hypertension"      
names(data)[names(data) == "diagstro"] <- "Stroke"

#time in days from diagnosis to death or censoring
t<-data$death_days-data$diag_date
tmax <- round(max(t[!is.na(t)]))+50
t[which(is.na(t))]<-tmax

#write.csv(data,"data_for_freq.csv")


res.cox <- coxph(Surv(t, status) ~ factor(race)+ factor(gender)+factor(rnd_group)+age+package_year+lesionsize
                 +factor(live_with_smoker)+factor(work_with_smoker)+factor(marital_status)+factor(smoke_pipe)+factor(smoke_sigar)
                 +factor(Asbestosis)+factor(Astma_adult)+factor(Astma_childhood)+factor(Bronchitis)+factor(Bronchitis_chronical)
                 +factor(COPD)+factor(Diabetes)+factor(Emphysema)+factor(Heart_disease)+factor(Pneumonia)+factor(Sarcoidosis)+factor(Tuberculosis)
                 +factor(Hypertension)+factor(Stroke)+progression_num+BMI+yrs_no_smoke+lung_dang_work
                 +alcohol_per_year+factor(stage), data =  data)
summary(res.cox)

#                                    coef  exp(coef)   se(coef)      z Pr(>|z|)    
#factor(race)2                  3.865e-01  1.472e+00  8.264e-01  0.468 0.639986    
#factor(race)3                 -1.832e+01  1.104e-08  3.229e+03 -0.006 0.995473    
#factor(race)4                 -1.086e+00  3.376e-01  1.503e+00 -0.723 0.469938    
#factor(gender)2               -4.612e-01  6.305e-01  3.886e-01 -1.187 0.235314    
#factor(rnd_group)2             9.441e-01  2.570e+00  3.328e-01  2.837 0.004559 ** 
#age                           -5.114e-02  9.501e-01  1.646e-01 -0.311 0.756044    
#package_year                   5.441e-03  1.005e+00  1.776e-01  0.031 0.975560    
#lesionsize                     1.678e-01  1.183e+00  1.629e-01  1.030 0.302806    
#factor(live_with_smoker)1     -1.090e+00  3.363e-01  5.548e-01 -1.964 0.049515 *  
#factor(work_with_smoker)1      9.008e-01  2.461e+00  6.239e-01  1.444 0.148791    
#factor(marital_status)1       -1.232e-01  8.841e-01  8.835e-01 -0.139 0.889128    
#factor(marital_status)2       -2.123e-01  8.087e-01  9.251e-01 -0.229 0.818498    
#factor(smoke_pipe)1            1.152e-01  1.122e+00  4.945e-01  0.233 0.815795    
#factor(smoke_sigar)1           1.192e-01  1.127e+00  4.378e-01  0.272 0.785406    
#factor(Asbestosis)1            1.852e+00  6.371e+00  1.125e+00  1.647 0.099625 .  
#factor(Astma_adult)1           6.985e-03  1.007e+00  8.066e-01  0.009 0.993090    
#factor(Astma_childhood)1       1.353e+00  3.870e+00  8.199e-01  1.651 0.098838 .  
#factor(Bronchitis)1            4.902e-02  1.050e+00  1.294e+00  0.038 0.969783    
#factor(Bronchitis_chronical)1 -1.582e+00  2.056e-01  7.308e-01 -2.165 0.030400 *  
#factor(COPD)1                 -1.547e-01  8.567e-01  7.802e-01 -0.198 0.842880    
#factor(Diabetes)1              7.936e-01  2.211e+00  6.441e-01  1.232 0.217951    
#factor(Emphysema)1             1.596e+00  4.932e+00  5.750e-01  2.775 0.005517 ** 
#factor(Heart_disease)1        -1.384e+00  2.507e-01  5.724e-01 -2.417 0.015642 *  
#factor(Pneumonia)1             1.398e-01  1.150e+00  4.297e-01  0.325 0.744869    
#factor(Sarcoidosis)1           2.471e+00  1.183e+01  1.114e+00  2.217 0.026617 *  
#factor(Tuberculosis)1          3.895e+00  4.913e+01  1.575e+00  2.473 0.013388 *  
#factor(Hypertension)1          5.558e-01  1.743e+00  3.271e-01  1.699 0.089308 .  
#factor(Stroke)1                2.754e+00  1.571e+01  8.245e-01  3.341 0.000836 ***
#progression_num                6.287e-01  1.875e+00  1.167e-01  5.389  7.1e-08 ***
#BMI                            1.926e-01  1.212e+00  1.807e-01  1.066 0.286594    
#yrs_no_smoke                  -2.368e-02  9.766e-01  1.644e-01 -0.144 0.885470    
#lung_dang_work                 1.626e-01  1.177e+00  1.705e-01  0.954 0.340054    
#alcohol_per_year              -9.481e-03  9.906e-01  1.804e-01 -0.053 0.958082    
#factor(stage)2                -6.948e-01  4.992e-01  6.268e-01 -1.108 0.267661    
#factor(stage)3                 5.018e-01  1.652e+00  5.228e-01  0.960 0.337185    
#factor(stage)4                 7.575e-03  1.008e+00  7.987e-01  0.009 0.992433    
#factor(stage)5                 1.390e+00  4.013e+00  4.068e-01  3.416 0.000635 ***
#factor(stage)6                 2.245e+00  9.445e+00  7.957e-01  2.822 0.004771 ** 
#factor(stage)7                 2.030e+00  7.611e+00  6.353e-01  3.195 0.001399 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#                              exp(coef) exp(-coef) lower .95 upper .95
#factor(race)2                 1.472e+00  6.794e-01   0.29135    7.4359
#factor(race)3                 1.104e-08  9.056e+07   0.00000       Inf
#factor(race)4                 3.376e-01  2.962e+00   0.01774    6.4221
#factor(gender)2               6.305e-01  1.586e+00   0.29440    1.3505
#factor(rnd_group)2            2.570e+00  3.890e-01   1.33879    4.9351
#age                           9.501e-01  1.052e+00   0.68813    1.3119
#package_year                  1.005e+00  9.946e-01   0.70988    1.4241
#lesionsize                    1.183e+00  8.455e-01   0.85950    1.6275
#factor(live_with_smoker)1     3.363e-01  2.973e+00   0.11338    0.9977
#factor(work_with_smoker)1     2.461e+00  4.063e-01   0.72469    8.3605
#factor(marital_status)1       8.841e-01  1.131e+00   0.15648    4.9953
#factor(marital_status)2       8.087e-01  1.237e+00   0.13193    4.9574
#factor(smoke_pipe)1           1.122e+00  8.912e-01   0.42573    2.9575
#factor(smoke_sigar)1          1.127e+00  8.876e-01   0.47767    2.6571
#factor(Asbestosis)1           6.371e+00  1.570e-01   0.70305   57.7395
#factor(Astma_adult)1          1.007e+00  9.930e-01   0.20721    4.8938
#factor(Astma_childhood)1      3.870e+00  2.584e-01   0.77591   19.3019
#factor(Bronchitis)1           1.050e+00  9.522e-01   0.08314   13.2675
#factor(Bronchitis_chronical)1 2.056e-01  4.865e+00   0.04908    0.8609
#factor(COPD)1                 8.567e-01  1.167e+00   0.18565    3.9535
#factor(Diabetes)1             2.211e+00  4.522e-01   0.62569    7.8147
#factor(Emphysema)1            4.932e+00  2.028e-01   1.59794   15.2194
#factor(Heart_disease)1        2.507e-01  3.989e+00   0.08165    0.7698
#factor(Pneumonia)1            1.150e+00  8.695e-01   0.49543    2.6697
#factor(Sarcoidosis)1          1.183e+01  8.451e-02   1.33184  105.1208
#factor(Tuberculosis)1         4.913e+01  2.035e-02   2.24410 1075.7900
#factor(Hypertension)1         1.743e+00  5.736e-01   0.91819    3.3102
#factor(Stroke)1               1.571e+01  6.365e-02   3.12179   79.0666
#progression_num               1.875e+00  5.333e-01   1.49189    2.3570
#BMI                           1.212e+00  8.248e-01   0.85077    1.7276
#yrs_no_smoke                  9.766e-01  1.024e+00   0.70758    1.3479
#lung_dang_work                1.177e+00  8.499e-01   0.84243    1.6433
#alcohol_per_year              9.906e-01  1.010e+00   0.69556    1.4107
#factor(stage)2                4.992e-01  2.003e+00   0.14612    1.7053
#factor(stage)3                1.652e+00  6.055e-01   0.59278    4.6019
#factor(stage)4                1.008e+00  9.925e-01   0.21061    4.8206
#factor(stage)5                4.013e+00  2.492e-01   1.80830    8.9069
#factor(stage)6                9.445e+00  1.059e-01   1.98571   44.9230
#factor(stage)7                7.611e+00  1.314e-01   2.19130   26.4372

#Concordance= 0.837  (se = 0.026 )
#Likelihood ratio test= 120.4  on 39 df,   p=3e-10
#Wald test            = 102.2  on 39 df,   p=1e-07
#Score (logrank) test = 157.6  on 39 df,   p=3e-16


#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

data<-read.csv("merged_clin_params.csv")

names(data)[names(data) == "age.x"] <- "age"
names(data)[names(data) == "de_stag_7thed"] <- "stage"
names(data)[names(data) == "rndgroup"] <- "rnd_group"
names(data)[names(data) == "de_stag_7thed"] <- "stage"
names(data)[names(data) == "diagasbe"] <- "Asbestosis"
names(data)[names(data) == "diagchro"] <- "Bronchitis_chronical" 
names(data)[names(data) == "diagdiab"] <- "Diabetes"   
names(data)[names(data) == "diagemph"] <- "Emphysema"   
names(data)[names(data) == "diaghear"] <- "Heart_disease"   
names(data)[names(data) == "diagpneu"] <- "Pneumonia"   
names(data)[names(data) == "diagsarc"] <- "Sarcoidosis"   
names(data)[names(data) == "diagtube"] <- "Tuberculosis"   
names(data)[names(data) == "diaghype"] <- "Hypertension"      
names(data)[names(data) == "diagstro"] <- "Stroke"
names(data)[names(data) == "smokelive"] <- "live_with_smoker"


data$status<-data$death_days
data$status[which(!is.na(data$status))]<-1
data$status[which(is.na(data$status))]<-0

t<-data$death_days-data$diag_date
tmax <- round(max(t[!is.na(t)]))+50
t[which(is.na(t))]<-tmax



data$age<-(data$age-mean(data$age))/sd(data$age)

data$str.given.lym<-(data$str.given.lym-mean(data$str.given.lym))/sd(data$str.given.lym)
data$tum.given.lym<-(data$tum.given.lym-mean(data$tum.given.lym))/sd(data$tum.given.lym)
data$lym.given.str<-(data$lym.given.str-mean(data$lym.given.str))/sd(data$lym.given.str)

data$tum.given.str<-(data$tum.given.str-mean(data$tum.given.str))/sd(data$tum.given.str)
data$lym.given.tum<-(data$lym.given.tum-mean(data$lym.given.tum))/sd(data$lym.given.tum)
data$str.given.tum<-(data$str.given.tum-mean(data$str.given.tum))/sd(data$str.given.tum)


data$pi_str<-(data$pi_str-mean(data$pi_str))/sd(data$pi_str)
data$pi_lym<-(data$pi_lym-mean(data$pi_lym))/sd(data$pi_lym)
data$lambda<-(data$lambda-mean(data$lambda))/sd(data$lambda)

data$stage[which(data$stage==2)]<-1
data$stage[which(data$stage==3)]<-2
data$stage[which(data$stage==4)]<-2
data$stage[which(data$stage==5)]<-3
data$stage[which(data$stage==6)]<-3
data$stage[which(data$stage==7)]<-4







#As in paper:
fit <- coxph(Surv(t, status) ~ factor(stage)+age+gender+cigsmok+  str.given.lym + tum.given.lym + lym.given.str + tum.given.str + lym.given.tum + str.given.tum + lambda + pi_lym + pi_str + cluster(pid), data = data)
print(summary(fit))
AIC(fit) #7 levels of stage:5988.665 , 4 levels of stage: 6014.315
extractAIC(fit) #5988.665
BIC(fit) #7 levels of stage: 6062.227, 4 levels of stage: 6075.617


#//Output for the model that is reported in the paper
#//the results of p-value significance are far from being the same, BUT:
#
#                   coef exp(coef) se(coef) robust se      z Pr(>|z|)    
#factor(stage)2  0.92551   2.52315  0.13283   0.43514  2.127   0.0334 *  
#factor(stage)3  1.46141   4.31202  0.12024   0.36457  4.009 6.11e-05 ***
#factor(stage)4  1.76959   5.86845  0.19015   0.70056  2.526   0.0115 *  
#age             0.09140   1.09571  0.04943   0.15954  0.573   0.5667    
#gender         -0.23066   0.79401  0.10167   0.32703 -0.705   0.4806    
#cigsmok        -0.10349   0.90168  0.09873   0.31080 -0.333   0.7391    
#str.given.lym   0.10498   1.11068  0.06047   0.09879  1.063   0.2880    
#tum.given.lym   0.11741   1.12458  0.06616   0.11503  1.021   0.3074    
#lym.given.str   0.13129   1.14029  0.11290   0.20840  0.630   0.5287    
#tum.given.str   0.32433   1.38311  0.11033   0.21947  1.478   0.1395    
#lym.given.tum  -0.10547   0.89990  0.09141   0.15273 -0.691   0.4898    
#str.given.tum  -0.21784   0.80425  0.08107   0.17224 -1.265   0.2060    
#lambda         -0.29003   0.74824  0.06143   0.15093 -1.922   0.0546 .  
#pi_lym          0.13836   1.14839  0.07855   0.16065  0.861   0.3891    
#pi_str         -0.70647   0.49338  0.17288   0.35973 -1.964   0.0495 *  

#               exp(coef) exp(-coef) lower .95 upper .95
#factor(stage)2    2.5232     0.3963    1.0754    5.9201
#factor(stage)3    4.3120     0.2319    2.1104    8.8106
#factor(stage)4    5.8685     0.1704    1.4866   23.1656
#age               1.0957     0.9126    0.8015    1.4979
#gender            0.7940     1.2594    0.4183    1.5073
#cigsmok           0.9017     1.1090    0.4903    1.6581
#str.given.lym     1.1107     0.9003    0.9152    1.3480
#tum.given.lym     1.1246     0.8892    0.8976    1.4090
#lym.given.str     1.1403     0.8770    0.7579    1.7156
#tum.given.str     1.3831     0.7230    0.8996    2.1265
#lym.given.tum     0.8999     1.1112    0.6671    1.2140
#str.given.tum     0.8043     1.2434    0.5738    1.1272
#lambda            0.7482     1.3365    0.5566    1.0058
#pi_lym            1.1484     0.8708    0.8382    1.5734
#pi_str            0.4934     2.0268    0.2438    0.9986

#Concordance= 0.716  (se = 0.037 )
#Likelihood ratio test= 293.7  on 15 df,   p=<2e-16
#Wald test            = 69.04  on 15 df,   p=7e-09
#Score (logrank) test = 347.6  on 15 df,   p=<2e-16,   Robust = 24.81  p=0.05

#----------------------------------------------------------------------------
#model with "stage" subdivided not on 4 levels, but on 7:
#----------------------------------------------------------------------------

#coef exp(coef) se(coef) robust se      z Pr(>|z|)    
#factor(stage)2  0.01032   1.01037  0.20951   0.64186  0.016 0.987177    
#factor(stage)3  0.61981   1.85857  0.15728   0.52642  1.177 0.239035    
#factor(stage)4  1.80425   6.07539  0.20882   0.59222  3.047 0.002315 ** 
#factor(stage)5  1.40265   4.06597  0.12942   0.39309  3.568 0.000359 ***
#factor(stage)6  2.28577   9.83323  0.32515   0.67559  3.383 0.000716 ***
#factor(stage)7  1.72262   5.59917  0.19394   0.71942  2.394 0.016645 *  
#age             0.09753   1.10245  0.05016   0.15516  0.629 0.529618    
#gender         -0.36068   0.69720  0.10798   0.34332 -1.051 0.293453    
#cigsmok        -0.11183   0.89420  0.10385   0.33078 -0.338 0.735300    
#str.given.lym   0.10161   1.10695  0.06363   0.10738  0.946 0.344022    
#tum.given.lym   0.10579   1.11159  0.06739   0.12134  0.872 0.383312    
#lym.given.str   0.09044   1.09466  0.11378   0.21560  0.420 0.674843    
#tum.given.str   0.34155   1.40713  0.10871   0.22026  1.551 0.120977    
#lym.given.tum  -0.06368   0.93830  0.09111   0.15063 -0.423 0.672475    
#str.given.tum  -0.22063   0.80201  0.08160   0.17472 -1.263 0.206675    
#lambda         -0.26213   0.76941  0.06341   0.16176 -1.621 0.105119    
#pi_lym          0.19974   1.22109  0.07842   0.15637  1.277 0.201466    
#pi_str         -0.65844   0.51766  0.17121   0.35618 -1.849 0.064511 .  

#exp(coef) exp(-coef) lower .95 upper .95
#factor(stage)2    1.0104     0.9897    0.2872     3.555
#factor(stage)3    1.8586     0.5380    0.6624     5.215
#factor(stage)4    6.0754     0.1646    1.9032    19.394
#factor(stage)5    4.0660     0.2459    1.8818     8.785
#factor(stage)6    9.8332     0.1017    2.6160    36.962
#factor(stage)7    5.5992     0.1786    1.3669    22.935
#age               1.1024     0.9071    0.8134     1.494
#gender            0.6972     1.4343    0.3557     1.366
#cigsmok           0.8942     1.1183    0.4676     1.710
#str.given.lym     1.1069     0.9034    0.8969     1.366
#tum.given.lym     1.1116     0.8996    0.8763     1.410
#lym.given.str     1.0947     0.9135    0.7174     1.670
#tum.given.str     1.4071     0.7107    0.9138     2.167
#lym.given.tum     0.9383     1.0658    0.6984     1.261
#str.given.tum     0.8020     1.2469    0.5695     1.130
#lambda            0.7694     1.2997    0.5604     1.056
#pi_lym            1.2211     0.8189    0.8988     1.659
#pi_str            0.5177     1.9318    0.2575     1.040

#Concordance= 0.721  (se = 0.037 )
#Likelihood ratio test= 325.4  on 18 df,   p=<2e-16
#Wald test            = 71.66  on 18 df,   p=2e-08
#Score (logrank) test = 447.2  on 18 df,   p=<2e-16,   Robust = 25.72  p=0.1


fit <- coxph(Surv(t, status) ~Stroke+Tuberculosis+Heart_disease+Emphysema+rnd_group+Bronchitis_chronical+ live_with_smoker+progression_num+factor(stage)+str.given.lym + tum.given.lym + lym.given.str + tum.given.str + lym.given.tum + str.given.tum + lambda + pi_lym + pi_str+ cluster(pid), data = data)
print(summary(fit))
AIC(fit) #5745.871
BIC(fit) #5839.867
#When comparing models fitted by maximum likelihood to the same data, the smaller the AIC or BIC, the better the fit.
#------------------------------------
#4 stages:
#------------------------------------
#                          coef exp(coef)  se(coef) robust se      z Pr(>|z|)    
#Stroke                0.739778  2.095471  0.336489  0.985521  0.751 0.452865    
#Tuberculosis          1.097625  2.997040  0.343593  0.548229  2.002 0.045271 *  
#Heart_disease        -0.786000  0.455664  0.157624  0.552462 -1.423 0.154817    
#Emphysema             1.646275  5.187617  0.183774  0.480930  3.423 0.000619 ***
#rnd_group             0.919044  2.506892  0.107723  0.371581  2.473 0.013386 *  
#Bronchitis_chronical -1.160641  0.313285  0.210776  0.684700 -1.695 0.090055 .  
#live_with_smoker      1.157111  3.180732  0.459537  1.148114  1.008 0.313533    
#progression_num       0.483801  1.622229  0.043426  0.135661  3.566 0.000362 ***
#factor(stage)2        1.155322  3.175045  0.135313  0.410493  2.814 0.004886 ** 
#factor(stage)3        1.455571  4.286932  0.122363  0.410317  3.547 0.000389 ***
#factor(stage)4        1.494277  4.456116  0.187715  0.730691  2.045 0.040853 *  
#str.given.lym         0.076368  1.079360  0.071522  0.113055  0.675 0.499359    
#tum.given.lym         0.073002  1.075732  0.078152  0.130081  0.561 0.574662    
#lym.given.str         0.196648  1.217316  0.122370  0.200155  0.982 0.325863    
#tum.given.str         0.187387  1.206094  0.120982  0.214522  0.874 0.382386    
#lym.given.tum        -0.104789  0.900514  0.103184  0.161238 -0.650 0.515755    
#str.given.tum        -0.088461  0.915339  0.088164  0.171076 -0.517 0.605098    
#lambda               -0.038096  0.962620  0.069888  0.162509 -0.234 0.814654    
#pi_lym               -0.008149  0.991884  0.096162  0.151665 -0.054 0.957152    
#pi_str               -0.567546  0.566915  0.190255  0.362358 -1.566 0.117288    

#                     exp(coef) exp(-coef) lower .95 upper .95
#Stroke                  2.0955     0.4772   0.30367    14.460
#Tuberculosis            2.9970     0.3337   1.02339     8.777
#Heart_disease           0.4557     2.1946   0.15431     1.346
#Emphysema               5.1876     0.1928   2.02116    13.315
#rnd_group               2.5069     0.3989   1.21017     5.193
#Bronchitis_chronical    0.3133     3.1920   0.08187     1.199
#live_with_smoker        3.1807     0.3144   0.33516    30.186
#progression_num         1.6222     0.6164   1.24348     2.116
#factor(stage)2          3.1750     0.3150   1.42016     7.098
#factor(stage)3          4.2869     0.2333   1.91816     9.581
#factor(stage)4          4.4561     0.2244   1.06412    18.660
#str.given.lym           1.0794     0.9265   0.86484     1.347
#tum.given.lym           1.0757     0.9296   0.83364     1.388
#lym.given.str           1.2173     0.8215   0.82230     1.802
#tum.given.str           1.2061     0.8291   0.79210     1.836
#lym.given.tum           0.9005     1.1105   0.65651     1.235
#str.given.tum           0.9153     1.0925   0.65458     1.280
#lambda                  0.9626     1.0388   0.70005     1.324
#pi_lym                  0.9919     1.0082   0.73682     1.335
#pi_str                  0.5669     1.7639   0.27866     1.153

#------------------------------------
#7 stages:
#------------------------------------

#                          coef exp(coef)  se(coef) robust se      z Pr(>|z|)    
#Stroke                0.531799  1.701991  0.358397  1.100890  0.483 0.629051    
#Tuberculosis          0.968043  2.632787  0.367589  0.512939  1.887 0.059127 .  
#Heart_disease        -0.771967  0.462103  0.159188  0.563247 -1.371 0.170510    
#Emphysema             1.648075  5.196967  0.186879  0.503514  3.273 0.001064 ** 
#rnd_group             0.853700  2.348319  0.112169  0.397716  2.147 0.031833 *  
#Bronchitis_chronical -1.147254  0.317507  0.216540  0.731217 -1.569 0.116656    
#live_with_smoker      1.101817  3.009630  0.462177  1.165290  0.946 0.344388    
#progression_num       0.506903  1.660142  0.044645  0.142275  3.563 0.000367 ***
#factor(stage)2        0.175613  1.191977  0.215636  0.682472  0.257 0.796932    
#factor(stage)3        1.079957  2.944554  0.160523  0.545964  1.978 0.047920 *  
#factor(stage)4        1.432138  4.187645  0.232606  0.569874  2.513 0.011968 *  
#factor(stage)5        1.430939  4.182623  0.130579  0.436762  3.276 0.001052 ** 
#factor(stage)6        2.397100 10.991254  0.313642  0.674001  3.557 0.000376 ***
#factor(stage)7        1.508757  4.521107  0.191499  0.724507  2.082 0.037300 *  
#str.given.lym         0.053344  1.054792  0.072916  0.114150  0.467 0.640276    
#tum.given.lym         0.093826  1.098369  0.079181  0.136645  0.687 0.492308    
#lym.given.str         0.215251  1.240173  0.122895  0.209798  1.026 0.304897    
#tum.given.str         0.191147  1.210637  0.120078  0.213474  0.895 0.370568    
#lym.given.tum        -0.086013  0.917582  0.103345  0.164320 -0.523 0.600662    
#str.given.tum        -0.087386  0.916323  0.088383  0.174999 -0.499 0.617532    
#lambda               -0.005545  0.994470  0.070668  0.174816 -0.032 0.974695    
#pi_lym                0.001428  1.001429  0.096306  0.155887  0.009 0.992693    
#pi_str               -0.566096  0.567738  0.188781  0.361357 -1.567 0.117213    

#                     exp(coef) exp(-coef) lower .95 upper .95
#Stroke                  1.7020    0.58755   0.19673    14.724
#Tuberculosis            2.6328    0.37983   0.96339     7.195
#Heart_disease           0.4621    2.16402   0.15322     1.394
#Emphysema               5.1970    0.19242   1.93713    13.943
#rnd_group               2.3483    0.42584   1.07701     5.120
#Bronchitis_chronical    0.3175    3.14953   0.07574     1.331
#live_with_smoker        3.0096    0.33227   0.30663    29.540
#progression_num         1.6601    0.60236   1.25615     2.194
#factor(stage)2          1.1920    0.83894   0.31286     4.541
#factor(stage)3          2.9446    0.33961   1.00994     8.585
#factor(stage)4          4.1876    0.23880   1.37054    12.795
#factor(stage)5          4.1826    0.23908   1.77696     9.845
#factor(stage)6         10.9913    0.09098   2.93316    41.187
#factor(stage)7          4.5211    0.22118   1.09281    18.704
#str.given.lym           1.0548    0.94805   0.84334     1.319
#tum.given.lym           1.0984    0.91044   0.84030     1.436
#lym.given.str           1.2402    0.80634   0.82206     1.871
#tum.given.str           1.2106    0.82601   0.79672     1.840
#lym.given.tum           0.9176    1.08982   0.66493     1.266
#str.given.tum           0.9163    1.09132   0.65026     1.291
#lambda                  0.9945    1.00556   0.70597     1.401
#pi_lym                  1.0014    0.99857   0.73778     1.359
#pi_str                  0.5677    1.76138   0.27961     1.153

#Concordance= 0.806  (se = 0.029 )
#Likelihood ratio test= 578.2  on 23 df,   p=<2e-16
#Wald test            = 100.1  on 23 df,   p=1e-11
#Score (logrank) test = 717.1  on 23 df,   p=<2e-16,   Robust = 38.7  p=0.02


fit <- coxph(Surv(t, status) ~Tuberculosis+Heart_disease+Emphysema+rnd_group+Bronchitis_chronical+ live_with_smoker+progression_num+factor(stage) + tum.given.lym + lym.given.str  + pi_str+ cluster(pid), data = data)
print(summary(fit))
AIC(fit) #5747.833
BIC(fit) # 5805.048

#Statistical significance. The column marked “z” gives the Wald statistic value. 
#It corresponds to the ratio of each regression coefficient to its standard error (z = coef/se(coef)).
#The wald statistic evaluates, whether the beta (β) coefficient of a given variable is statistically 
#significantly different from 0. From the output above, we can conclude that the variable sex have highly
#statistically significant coefficients.

#The regression coefficients. The second feature to note in the Cox model results is the the sign of the
#regression coefficients (coef). A positive sign means that the hazard (risk of death) is higher, and thus
#the prognosis worse, for subjects with higher values of that variable. The variable sex is encoded as a 
#numeric vector. 1: male, 2: female. The R summary for the Cox model gives the hazard ratio (HR) for the 
#second group relative to the first group, that is, female versus male. The beta coefficient for sex = -0.53 
#indicates that females have lower risk of death (lower survival rates) than males, in these data.

#Hazard ratios. The exponentiated coefficients (exp(coef) = exp(-0.53) = 0.59), also known as hazard ratios,
#give the effect size of covariates. For example, being female (sex=2) reduces the hazard by a factor of 0.59, 
#or 41%. Being female is associated with good prognostic.




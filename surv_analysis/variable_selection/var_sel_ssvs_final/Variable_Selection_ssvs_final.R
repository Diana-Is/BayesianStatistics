library(rjags)
library(coda)

# for plots
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
require(gplots)
require(ggpubr)

setwd("C:/Users/Diana/Desktop/Bayes_project/var_sel_ssvs_final")

data <- read.csv("total_merge_pure_params.csv", header=T)

pids<-unique(data$pid)
filter1<-data[which(data$pid==pids[1]),]
rand_n <- sample(1:dim(filter1)[1], 1)
filter1<-filter1[rand_n,]
for(i in 2:length(pids)) {
  temp<-data[which(data$pid==pids[i]),]
  rand_n <- sample(1:dim(temp)[1], 1)
  temp<-temp[rand_n,]
  filter1<-rbind(filter1,temp)
  
}

data<-filter1

names(data)[names(data) == "age.x"] <- "age"
names(data)[names(data) == "rndgroup"] <- "rnd_group"
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
names(data)[names(data) == "diaghear"] <- "Heart disease"   
names(data)[names(data) == "diagpneu"] <- "Pneumonia"   
names(data)[names(data) == "diagsarc"] <- "Sarcoidosis"   
names(data)[names(data) == "diagtube"] <- "Tuberculosis"   
names(data)[names(data) == "diaghype"] <- "Hypertension"      
names(data)[names(data) == "diagstro"] <- "Stroke"

names(data)[names(data) == "lym"] <- "omega_lym"
names(data)[names(data) == "str"] <- "omega_str"
names(data)[names(data) == "lym_lym"] <- "theta_lym_lym"
names(data)[names(data) == "lym_str"] <- "theta_lym_str"
names(data)[names(data) == "lym_tum"] <- "theta_lym_tum"
names(data)[names(data) == "str_str"] <- "theta_str_str"
names(data)[names(data) == "str_tum"] <- "theta_str_tum"
#-------------------------------------------------------------------
#Standartization of non-categorical variables
#-------------------------------------------------------------------
data$package_year<-(data$package_year-mean(data$package_year))/sd(data$package_year)
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
data$rnd_group[which(data$rnd_group==1)]<-0 #Computer Tomography
data$rnd_group[which(data$rnd_group==2)]<-1 #X-Ray
data$gender[which(data$gender==1)]<-0 #male
data$gender[which(data$gender==2)]<-1 #female

data$race_asian<-data$race
data$race_asian[which(data$race_asian!=3)]<-0
data$race_asian[which(data$race_asian==3)]<-1

data$race_minor<-data$race
data$race_minor[which(data$race_minor!=4)]<-0
data$race_minor[which(data$race_minor==4)]<-1
data$race[which(data$race!=2)]<-0
data$race[which(data$race==2)]<-1


data$stage_4<-data$de_stag_7thed
data$stage_4[which(data$stage_4!=7)]<-0
data$stage_4[which(data$stage_4==7)]<-1

data$stage_3_6<-data$de_stag_7thed
data$stage_3_6[which(data$stage_3_6!=6)]<-0
data$stage_3_6[which(data$stage_3_6==6)]<-1

data$stage_3_5<-data$de_stag_7thed
data$stage_3_5[which(data$stage_3_5!=5)]<-0
data$stage_3_5[which(data$stage_3_5==5)]<-1

data$stage_3<-data$stage_3_5+data$stage_3_6
data[c("stage_3_5","stage_3_6")]<-list(NULL)

data$stage_2_4<-data$de_stag_7thed
data$stage_2_4[which(data$stage_2_4!=4)]<-0
data$stage_2_4[which(data$stage_2_4==4)]<-1

data$stage_2_3<-data$de_stag_7thed
data$stage_2_3[which(data$stage_2_3!=3)]<-0
data$stage_2_3[which(data$stage_2_3==3)]<-1

data$stage_2<-data$stage_2_4+data$stage_2_3
data[c("stage_2_4","stage_2_3","de_stag_7thed")]<-list(NULL)



#time in days from diagnosis to death or censoring
t<-data$death_days-data$diag_date
cens<-data$fup_days-data$diag_date
tmax <- round(max(t[!is.na(t)]))+1
cens[!is.na(t)] <- tmax
#removing irrelevant variables
data[c("pid","cigsmok","fup_days","death_days","diag_date","X","event","marital")]<-list(NULL)
data[c("omega_lym","omega_str","theta_lym_lym","theta_lym_str","theta_lym_tum","theta_str_str","theta_str_tum","lambda")]<-list(NULL)
#reordering the variables
#To change the order as in the above question do df2[,c(1,3,2,4)]
#write.csv(data,"var_3.csv")
data<-data[,c(1,30,31,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,34,33,32)]

X <- as.matrix(data)

N <- dim(X)[1]
p <- dim(X)[2]


########################################################################
#################### Spike and Slab Weibull SSVS########################
########################################################################

###Parameters of the  spike slab prior Set 1
c_ss <- 100
intersect <- 0.05
tau_ss <- intersect / sqrt(2 * log(c_ss) * c_ss^2/(c_ss^2 - 1))

## With this choice of hyperameter, c_ss and intersection 
## the variance of the quasi-spike prior is
tau_ss^2

##While the variance of the slab is 
(tau_ss*c_ss)^2

## data to pass to JAGS (see the code in SSVS_probit.bug)
data_JAGS_SSVS <-list(N = N, p = p, t = t, cens=cens, X = as.matrix(X), 
                   tau_ss = tau_ss, c_ss = c_ss)

## A list of initial value for the MCMC algorithm 
# that WinBUGS will implement
inits = function() {
  list(beta0 = 0.0, beta = rep(0,p), g = rep(0,p),
       .RNG.seed = 321, .RNG.name = 'base::Wichmann-Hill') }

model=jags.model("SSVS_Weibull.bug",
                 data = data_JAGS_SSVS,
                 n.adapt = 10000,
                 inits = inits,
                 n.chains = 2) 

# if we want to perform a larger burn in with not adaptation.
update(model,n.iter=50000)

  # Posterior parameters WinBUGS  has to save NO "sigma2" here
param <- c("beta0", "beta", "g", "mdl")

# number of iterations & thinning
nit <- 70000
thin <-10

##The command coda.samle() calls jags from R passing the data and initial value just defined
output <- coda.samples(model = model,
                       variable.names = param,
                       n.iter = nit,
                       thin = thin)
#write.csv(data,"for_var_sel2.csv")
save(output,file='SpSl_weib_ssvs_final.dat')#save the chain
load('SpSl_weib_ssvs_final.dat')


###We give a look at the trace plot and at the density summary
### of the posterior chains for each of the parameters
#x11()
#plot(output,ask=T)
#dev.off()

# summary
#summary(output)

### To work with the posterior chains it is better 
### to cast the output to be an array object of R

output2 <- as.matrix(output)

###Some variable selection thecniques:
# The median probability moodel (MPM)
# pick variables with estimated posterior inclusion probabilities 
# higher than 0.5
# Notice that the estimated posterior inclusion probabilities are the
# posterior means of the gamma variables (in the code we called g)
head(output2)

##We save the posterior chain of the inclusion variable in post_g
post_g <-as.matrix(output2[,36:69])
post_mean_g <- apply(post_g,2,"mean") 


p1 <- data.frame(value = post_mean_g, var = colnames(X)) %>%
  ggplot(aes(y = value, x = var, fill = var)) + 
  geom_bar(color="blue", fill="#B2B2FE",stat="identity") + 
  geom_hline(mapping = aes(yintercept = .5), col = 2, lwd = 1.1) +
  ylim(0,1)+
  coord_flip() + 
  theme_minimal() + 
 theme(legend.position="none",axis.text = element_text(size = 13),axis.title = element_text(size = 20)) + 
  ylab("Posterior inclusion probabilities") + 
  xlab("")
x11()
p1




# ANOTHER WAY TO CHOOSE THE MODEL:
# Highest posterior density model (HPD)
# pick a model with the highest estimated posterior probability 


## HPD model 
# getting the unique profiles
# and sort the results
unique_model <- unique(post_g, MARGIN  = 1)
freq <- apply(unique_model, 1, function(b) sum(apply(post_g, MARGIN = 1, function(a) all(a == b))))
cbind(unique_model[order(freq,decreasing = T),], sort(freq,decreasing = T))

# the HPD model is 
colnames(X)[as.logical(unique_model[which.max(freq),])]

#[1] "race_asian"      "gender"          "rnd_group"      
#[4] "Asbestosis"      "progression_num" "stage_4"


########################################################################
#################### Spike and Slab lognorm#############################
########################################################################
###Parameters of the  spike slab prior Set 1
c_ss <- 100
intersect <- 0.05
tau_ss <- intersect / sqrt(2 * log(c_ss) * c_ss^2/(c_ss^2 - 1))

## With this choice of hyperameter, c_ss and intersection 
## the variance of the quasi-spike prior is
tau_ss^2

##While the variance of the slab is 
(tau_ss*c_ss)^2

## data to pass to JAGS (see the code in SSVS_probit.bug)
data_JAGS_SSVS <-list(N = N, p = p, t = t, cens=cens, X = as.matrix(X), 
                      tau_ss = tau_ss, c_ss = c_ss)

## A list of initial value for the MCMC algorithm 
# that WinBUGS will implement
inits = function() {
  list(beta0 = 0.0, beta = rep(0,p), g = rep(0,p),
       .RNG.seed = 321, .RNG.name = 'base::Wichmann-Hill') }


model=jags.model("SSVS_lognorm.bug",
                 data = data_JAGS_SSVS,
                 n.adapt = 10000,
                 inits = inits,
                 n.chains = 2) 

# if we want to perform a larger burn in with not adaptation.
update(model,n.iter=50000)

# Posterior parameters WinBUGS  has to save NO "sigma2" here
param <- c("beta0", "beta", "g", "mdl")

# number of iterations & thinning
nit <- 70000
thin <-10

##The command coda.samle() calls jags from R passing the data and initial value just defined
output_lognorm <- coda.samples(model = model,
                       variable.names = param,
                       n.iter = nit,
                       thin = thin)
save(output_lognorm,file='SpSl_lognorm_SSVS_final.dat')#save the chain
load('SpSl_lognorm_SSVS_final.dat')


### To work with the posterior chains it is better 
### to cast the output to be an array object of R

output_lognorm2 <- as.matrix(output_lognorm)

###Some variable selection thecniques:
# The median probability moodel (MPM)
# pick variables with estimated posterior inclusion probabilities 
# higher than 0.5
# Notice that the estimated posterior inclusion probabilities are the
# posterior means of the gamma variables (in the code we called g)

##We save the posterior chain of the inclusion variable in post_g
post_g_norm <-as.matrix(output_lognorm2[,36:69])
post_mean_g_norm <- apply(post_g_norm,2,"mean") 


p1 <- data.frame(value = post_mean_g_norm, var = colnames(X)) %>%
  ggplot(aes(y = value, x = var, fill = var)) + 
  geom_bar(color="blue", fill="#B2B2FE",stat="identity") + 
  geom_hline(mapping = aes(yintercept = .5), col = 2, lwd = 1.1) +
  ylim(0,1)+
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="none",axis.text = element_text(size = 13),axis.title = element_text(size = 20)) + 
  ylab("Posterior inclusion probabilities") + 
  xlab("")
x11()
p1


## HPD model 
# getting the unique profiles
# and sort the results
unique_model <- unique(post_g_norm, MARGIN  = 1)
freq <- apply(unique_model, 1, function(b) sum(apply(post_g_norm, MARGIN = 1, function(a) all(a == b))))
cbind(unique_model[order(freq,decreasing = T),], sort(freq,decreasing = T))
colnames(X)[as.logical(unique_model[which.max(freq),])]

#[1] "race_asian"       "race_minor"       "gender"          
#[4] "rnd_group"        "work_with_smoker" "Asbestosis"      
#[7] "Tuberculosis"     "progression_num"  "BMI"             
#[10] "stage_4"         


########################################################################
#################### Spike and Slab loglogist###########################
########################################################################
c_ss <- 100
intersect <- 0.05
tau_ss <- intersect / sqrt(2 * log(c_ss) * c_ss^2/(c_ss^2 - 1))

## With this choice of hyperameter, c_ss and intersection 
## the variance of the quasi-spike prior is
tau_ss^2

##While the variance of the slab is 
(tau_ss*c_ss)^2

logt<-log(t)
logc<-log(cens)

## data to pass to JAGS (see the code in SSVS_probit.bug)
data_JAGS_SSVS <-list(N = N, p = p, t = logt, cens=logc, X = as.matrix(X), 
                      tau_ss = tau_ss, c_ss = c_ss)

## A list of initial value for the MCMC algorithm 
# that WinBUGS will implement
inits = function() {
  list(beta0 = 0.0, beta = rep(0,p), g = rep(0,p),
       .RNG.seed = 312, .RNG.name = 'base::Wichmann-Hill') }


model=jags.model("loglog_ssvs.bug",
                 data = data_JAGS_SSVS,
                 n.adapt = 10000,
                 inits = inits,
                 n.chains = 2) 

# if we want to perform a larger burn in with not adaptation.
update(model,n.iter=50000)

# Posterior parameters WinBUGS  has to save NO "sigma2" here
param <- c("beta0", "beta", "g", "mdl")

# number of iterations & thinning
nit <- 70000
thin <-10

##The command coda.samle() calls jags from R passing the data and initial value just defined
output_loglog <- coda.samples(model = model,
                               variable.names = param,
                               n.iter = nit,
                               thin = thin)
save(output_loglog,file='SpSl_loglog_ssvs_final.dat')#save the chain
load('SpSl_loglog_ssvs_final.dat')


### To work with the posterior chains it is better 
### to cast the output to be an array object of R

output_loglog2 <- as.matrix(output_loglog)

###Some variable selection thecniques:
# The median probability moodel (MPM)
# pick variables with estimated posterior inclusion probabilities 
# higher than 0.5
# Notice that the estimated posterior inclusion probabilities are the
# posterior means of the gamma variables (in the code we called g)


##We save the posterior chain of the inclusion variable in post_g
post_g_log <-as.matrix(output_loglog2[,36:69])
post_mean_g_log <- apply(post_g_log,2,"mean") 



p1 <- data.frame(value = post_mean_g_log, var = colnames(X)) %>%
  ggplot(aes(y = value, x = var, fill = var)) + 
  geom_bar(color="blue", fill="#B2B2FE",stat="identity") + 
  geom_hline(mapping = aes(yintercept = .5), col = 2, lwd = 1.1) +
  ylim(0,1)+
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="none",axis.text = element_text(size = 13),axis.title = element_text(size = 20)) + 
  ylab("Posterior inclusion probabilities") + 
  xlab("")
x11()
p1

unique_model <- unique(post_g_log, MARGIN  = 1)
freq <- apply(unique_model, 1, function(b) sum(apply(post_g_log, MARGIN = 1, function(a) all(a == b))))
cbind(unique_model[order(freq,decreasing = T),], sort(freq,decreasing = T))
colnames(X)[as.logical(unique_model[which.max(freq),])]

#[1] "race_asian"       "gender"           "rnd_group"       
#[4] "work_with_smoker" "Asbestosis"       "Sarcoidosis"     
#[7] "Tuberculosis"     "progression_num"  "stage_4"     

r<-data.frame(value = post_mean_g, var = colnames(X))
r$base<-"Weibull"

r3<-data.frame(value = post_mean_g_norm, var = colnames(X))
r3$base<-"LogNorm"

r2<-data.frame(value = post_mean_g_log, var = colnames(X))
r2$base<-"LogLog"

rr<-rbind(r,r3,r2)


a<-ggplot(rr, aes(x=var,y=value,fill=base)) + 
  geom_bar(position="dodge",stat="identity") +
  coord_flip()+
  ylim(0,1)+
  ylab("Posterior inclusion probabilities") + 
  xlab("")+
  theme_minimal()+
theme(axis.text = element_text(size = 13),axis.title = element_text(size = 20))
x11()
a
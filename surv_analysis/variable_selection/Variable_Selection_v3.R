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

setwd("C:/Users/Diana/Desktop/Bayes_project/")

data <- read.csv("data_cleaned_add.csv", header=T)

#-------------------------------------------------------------------
#Standartization of non-categorical variables
#-------------------------------------------------------------------
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
data$rndgroup[which(data$rndgroup==1)]<-0 #Computer Tomography
data$rndgroup[which(data$rndgroup==2)]<-1 #X-Ray
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

data$mar_sep_wid<-data$marital
data$mar_sep_wid[which(data$mar_sep_wid!=2)]<-0
data$mar_sep_wid[which(data$mar_sep_wid==2)]<-1

data$marital[which(data$marital==2)]<-0

data$stage_7<-data$de_stag_7thed
data$stage_7[which(data$stage_7!=7)]<-0
data$stage_7[which(data$stage_7==7)]<-1

data$stage_6<-data$de_stag_7thed
data$stage_6[which(data$stage_6!=6)]<-0
data$stage_6[which(data$stage_6==6)]<-1

data$stage_5<-data$de_stag_7thed
data$stage_5[which(data$stage_5!=5)]<-0
data$stage_5[which(data$stage_5==5)]<-1

data$stage_4<-data$de_stag_7thed
data$stage_4[which(data$stage_4!=4)]<-0
data$stage_4[which(data$stage_4==4)]<-1

data$stage_3<-data$de_stag_7thed
data$stage_3[which(data$stage_3!=3)]<-0
data$stage_3[which(data$stage_3==3)]<-1

data$de_stag_7thed[which(data$de_stag_7thed!=2)]<-0
data$de_stag_7thed[which(data$de_stag_7thed==2)]<-1


names(data)[names(data) == "rndgroup"] <- "rnd_group"
names(data)[names(data) == "marital"] <- "marital_status"
names(data)[names(data) == "de_stag_7thed"] <- "stage"
names(data)[names(data) == "smokelive"] <- "live_with_smoker"
names(data)[names(data) == "smokework"] <- "work_with_smoker"
names(data)[names(data) == "pkyr"] <- "package-year"
names(data)[names(data) == "pipe"] <- "smoke_pipe"
names(data)[names(data) == "cigar"] <- "smoke_sigar"
names(data)[names(data) == "diagasbe"] <- "Asbestosis"
names(data)[names(data) == "diagadas"] <- "Astma_adult"
names(data)[names(data) == "diagchas"] <- "Astma_childhood"
names(data)[names(data) == "diagbron"] <- "Bronchitis" 
names(data)[names(data) == "diagchro"] <- "Bronchitis chronical" 
names(data)[names(data) == "diagcopd"] <- "COPD"   
names(data)[names(data) == "diagdiab"] <- "Diabetes"   
names(data)[names(data) == "diagemph"] <- "Emphysema"   
names(data)[names(data) == "diaghear"] <- "Heart disease"   
names(data)[names(data) == "diagpneu"] <- "Pneumonia"   
names(data)[names(data) == "diagsarc"] <- "Sarcoidosis"   
names(data)[names(data) == "diagtube"] <- "Tuberculosis"   
names(data)[names(data) == "diaghype"] <- "Hypertension"      
names(data)[names(data) == "diagstro"] <- "Stroke"

#time in days from diagnosis to death or censoring
t<-data$death_days-data$diag_date
cens<-data$fup_days-data$diag_date
tmax <- round(max(t[!is.na(t)]))+1
cens[!is.na(t)] <- tmax
#removing irrelevant variables
data[c("pid","cigsmok","fup_days","death_days","diag_date","X","event")]<-list(NULL)

#reordering the variables
#To change the order as in the above question do df2[,c(1,3,2,4)]
#write.csv(data,"var_3.csv")
data<-data[,c(1,32,33,2,3,4,5,6,7,8,9,34,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,39,38,37,36,35,27,28,29,30,31)]

X <- as.matrix(data)

N <- dim(X)[1]
p <- dim(X)[2]



########################################################################
#################### Spike and Slab Weibull#############################
########################################################################

data_JAGS_SpSl <-list(N = N, p = p, X = as.matrix(X), 
                      var_beta = rep(1, p),alpha=1,cens=cens,t=t)

## A list of initial value for the MCMC algorithm 
# that WinBUGS will implement
inits = function() {
  list(beta0 = 0.0, beta_temp = rep(0,p), g = rep(0,p), theta = rep(0.5, p),
       .RNG.seed = 48, .RNG.name = 'base::Wichmann-Hill') 
}

model=jags.model("SpSl_Weibull.bug",
                 data = data_JAGS_SpSl,
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
save(output,file='SpSl_weib_new3.dat')#save the chain
load('SpSl_weib_new3.dat')


# the output is an mcmc object of the library coda
str(output)

###We give a look at the trace plot and at the density summary
### of the posterior chains for each of the parameters
x11()
plot(output,ask=T)
dev.off()

# summary
summary(output)

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
post_g <-as.matrix(output2[,41:79])
apply(post_g,2,"mean")
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

# we will compare later the model with other methods
mp_SpSl <- as.vector(which(post_mean_g >= 0.5))
post_mean_g[mp_SpSl]



# ANOTHER WAY TO CHOOSE THE MODEL:
# Highest posterior density model (HPD)
# pick a model with the highest estimated posterior probability 

# Recall that we have represented the model index
# using a binary coding as 
# mdl=1+2^g1+...+2^gp
# the visited model are saved in the chain

plot(output[,"mdl"], pch = 20)

# for example at iteration 10 the chain explored the 
# model 

output[10,"mdl"]

### We start to analyze how many models have been visited 
## by the posterior chain:
length(unique( output2[,"mdl"]))

# model visited out of 4096
## Now we compute the posterior frequency of the visited  models
visited_models<-table(output2[,"mdl"])
visited_models

## HPD model 
# getting the unique profiles
# and sort the results
unique_model <- unique(post_g, MARGIN  = 1)
freq <- apply(unique_model, 1, function(b) sum(apply(post_g, MARGIN = 1, function(a) all(a == b))))
cbind(unique_model[order(freq,decreasing = T),], sort(freq,decreasing = T))

# the HPD model is 
colnames(X)[as.logical(unique_model[which.max(freq),])]
HDP_SpSl <- c(1:21)[as.logical(unique_model[which.max(freq),])]



u<-summary(output)
for(l in 1:p){
  if(u$quantiles[l,1]<0 && u$quantiles[l,5]>0)
  {
    # cat("*** variable ", colnames(X)[l], " excluded \n")
  }
  else
  {
    cat("*** variable ", colnames(X)[l], " included \n")
  }
}

########################################################################
#################### Spike and Slab lognorm#############################
########################################################################

data_JAGS_SpSl <-list(N = N, p = p, X = as.matrix(X), 
                      var_beta = rep(1, p),cens=cens,t=t)

## A list of initial value for the MCMC algorithm 
# that WinBUGS will implement
inits = function() {
  list(beta0 = 0.0, beta_temp = rep(0,p), g = rep(0,p), theta = rep(0.5, p),
       .RNG.seed = 48, .RNG.name = 'base::Wichmann-Hill') 
}

model=jags.model("SpSl_lognorm.bug",
                 data = data_JAGS_SpSl,
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
save(output_lognorm,file='SpSl_lognorm_new2.dat')#save the chain
load('SpSl_lognorm_new2.dat')


x11()
plot(output_lognorm,ask=T)
dev.off()

# summary
summary(output_lognorm)

### To work with the posterior chains it is better 
### to cast the output to be an array object of R

output_lognorm2 <- as.matrix(output_lognorm)

###Some variable selection thecniques:
# The median probability moodel (MPM)
# pick variables with estimated posterior inclusion probabilities 
# higher than 0.5
# Notice that the estimated posterior inclusion probabilities are the
# posterior means of the gamma variables (in the code we called g)
head(output_lognorm2)

##We save the posterior chain of the inclusion variable in post_g
post_g <-as.matrix(output_lognorm2[,41:79])
apply(post_g,2,"mean")
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


## HPD model 
# getting the unique profiles
# and sort the results
unique_model <- unique(post_g, MARGIN  = 1)
freq <- apply(unique_model, 1, function(b) sum(apply(post_g, MARGIN = 1, function(a) all(a == b))))
cbind(unique_model[order(freq,decreasing = T),], sort(freq,decreasing = T))
colnames(X)[as.logical(unique_model[which.max(freq),])]

#[1] "race_minor"    "package-year"  "lesionsize"   
#[4] "mar_sep_wid"   "Astma_adult"   "Emphysema"    
#[7] "Heart disease" "Sarcoidosis"   "Hypertension" 
#[10] "Stroke"        "BMI"  

########################################################################
#################### Spike and Slab loglogist###########################
########################################################################

logt<-log(t)
logc<-log(cens)

data_JAGS_SpSl <-list(N = N, p = p, X = as.matrix(X), 
                      var_beta = rep(1, p),cens=logc,t=logt)

## A list of initial value for the MCMC algorithm 
# that WinBUGS will implement
inits = function() {
  list(beta0 = 0.0, beta_temp = rep(0,p), g = rep(0,p), theta = rep(0.5, p),
       .RNG.seed = 48, .RNG.name = 'base::Wichmann-Hill') 
}

model=jags.model("SpSl_loglog.bug",
                 data = data_JAGS_SpSl,
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
save(output_loglog,file='SpSl_loglog_new2.dat')#save the chain
load('SpSl_loglog_new2.dat')


x11()
plot(output_loglog,ask=T)
dev.off()

# summary
summary(output_loglog)

### To work with the posterior chains it is better 
### to cast the output to be an array object of R

output_loglog2 <- as.matrix(output_loglog)

###Some variable selection thecniques:
# The median probability moodel (MPM)
# pick variables with estimated posterior inclusion probabilities 
# higher than 0.5
# Notice that the estimated posterior inclusion probabilities are the
# posterior means of the gamma variables (in the code we called g)
head(output_loglog2)

##We save the posterior chain of the inclusion variable in post_g
post_g <-as.matrix(output_loglog2[,41:79])
apply(post_g,2,"mean")
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

unique_model <- unique(post_g, MARGIN  = 1)
freq <- apply(unique_model, 1, function(b) sum(apply(post_g, MARGIN = 1, function(a) all(a == b))))
cbind(unique_model[order(freq,decreasing = T),], sort(freq,decreasing = T))
colnames(X)[as.logical(unique_model[which.max(freq),])]


#[1] "race_asian"      "race_minor"      "mar_sep_wid"    
#[4] "Astma_adult"     "Bronchitis"      "COPD"           
#[7] "Emphysema"       "Heart disease"   "Sarcoidosis"    
#[10] "Hypertension"    "stage_5"         "progression_num"


# R & STAN are friends!
library(rjags)
library(coda)
library(survival)
#McGilchrist and Aisbett (1991) kidney
# for plots
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
require(gplots)
require(ggpubr)

###################################################################
#SURV RANDOM EFF WEIBULL
###################################################################
setwd("C:/Users/Diana/Desktop/Bayes_project/random_full")
data <- read.csv("total_merge_pure_params.csv", header=T)

myvars <- c("pid","death_days","diag_date","fup_days", 
            "gender", "race","smokework","rndgroup","diagasbe",
            "diagsarc","de_stag_7thed","progression_num","lym_lym",
            "lym_str","lym_tum","str_str","str_tum","lym","str","lambda")

data <- data[myvars]

names(data)[names(data) == "rndgroup"] <- "rnd_group"
names(data)[names(data) == "de_stag_7thed"] <- "stage"
names(data)[names(data) == "diagasbe"] <- "Asbestosis"
names(data)[names(data) == "diagsarc"] <- "Sarcoidosis"
names(data)[names(data) == "smokework"] <- "work_with_smoker"
names(data)[names(data) == "lym"] <- "omega_lym"
names(data)[names(data) == "str"] <- "omega_str"
names(data)[names(data) == "lym_lym"] <- "theta_lym_lym"
names(data)[names(data) == "lym_str"] <- "theta_lym_str"
names(data)[names(data) == "lym_tum"] <- "theta_lym_tum"
names(data)[names(data) == "str_str"] <- "theta_str_str"
names(data)[names(data) == "str_tum"] <- "theta_str_tum"

data$rnd_group[which(data$rnd_group==1)]<-0 #Computer Tomography
data$rnd_group[which(data$rnd_group==2)]<-1 #X-Ray

data$gender[which(data$gender==1)]<-0 #male 
data$gender[which(data$gender==2)]<-1 #female

data$race_asian<-data$race
data$race_asian[which(data$race_asian!=3)]<-0 
data$race_asian[which(data$race_asian==3)]<-1

data$stage[which(data$stage==2)]<-1
data$stage[which(data$stage==3)]<-2
data$stage[which(data$stage==4)]<-2
data$stage[which(data$stage==5)]<-3
data$stage[which(data$stage==6)]<-3
data$stage[which(data$stage==7)]<-4


pids<-unique(data$pid)
filter5<-data[which(data$pid==pids[1]),]
filter5<-filter5[1:5,]

for(i in 2:length(pids)) {
  temp<-data[which(data$pid==pids[i]),]
  if (dim(temp)[1]>5) temp<-temp[1:5,]
  filter5<-rbind(filter5,temp)
  
}

data<-filter5

#-------------------------------------------------------------------
#Standartization of non-categorical variables
#-------------------------------------------------------------------
t<-data$death_days-data$diag_date
cens<-data$fup_days-data$diag_date
tmax <- round(max(t[!is.na(t)]))+50
cens[!is.na(t)] <- tmax
#removing irrelevant variables
data[c("fup_days","death_days","diag_date","race")]<-list(NULL)


data$theta_lym_lym<-(data$theta_lym_lym-mean(data$theta_lym_lym))/sd(data$theta_lym_lym)
data$theta_lym_str<-(data$theta_lym_str-mean(data$theta_lym_str))/sd(data$theta_lym_str)
data$theta_lym_tum<-(data$theta_lym_tum-mean(data$theta_lym_tum))/sd(data$theta_lym_tum)
data$theta_str_str<-(data$theta_str_str-mean(data$theta_str_str))/sd(data$theta_str_str)
data$theta_str_tum<-(data$theta_str_tum-mean(data$theta_str_tum))/sd(data$theta_str_tum)

data$omega_str<-(data$omega_str-mean(data$omega_str))/sd(data$omega_str)
data$omega_lym<-(data$omega_lym-mean(data$omega_lym))/sd(data$omega_lym)
data$lambda<-(data$lambda-mean(data$lambda))/sd(data$lambda)


N <- dim(data)[1]/5
M<-5
t<-t(matrix(t,nrow = 5,ncol = 174))
cens<-t(matrix(cens,nrow = 5,ncol = 174))
s_theta_lym_lym<-t(matrix(data$theta_lym_lym,nrow = 5,ncol = 174))
s_theta_lym_str<-t(matrix(data$theta_lym_str,nrow = 5,ncol = 174))
s_theta_lym_tum<-t(matrix(data$theta_lym_tum,nrow = 5,ncol = 174))
s_theta_str_str<-t(matrix(data$theta_str_str,nrow = 5,ncol = 174))
s_theta_str_tum<-t(matrix(data$theta_str_tum,nrow = 5,ncol = 174))
s_omega_str<-t(matrix(data$omega_str,nrow = 5,ncol = 174))
s_omega_lym<-t(matrix(data$omega_lym,nrow = 5,ncol = 174))
lambda<-t(matrix(data$lambda,nrow = 5,ncol = 174))

work_with_smoker<-t(matrix(data$work_with_smoker,nrow = 5,ncol = 174))[,1]
stage<-t(matrix(data$stage,nrow = 5,ncol = 174))[,1]
rnd_group<-t(matrix(data$rnd_group,nrow = 5,ncol = 174))[,1]
progression_num<-t(matrix(data$progression_num,nrow = 5,ncol = 174))[,1]
Asbestosis<-t(matrix(data$Asbestosis,nrow = 5,ncol = 174))[,1]
Sarcoidosis<-t(matrix(data$Sarcoidosis,nrow = 5,ncol = 174))[,1]
gender<-t(matrix(data$gender,nrow = 5,ncol = 174))[,1]
race_asian<-t(matrix(data$race_asian,nrow = 5,ncol = 174))[,1]

progression_num<-(progression_num-mean(progression_num))/sd(progression_num)


dat <- list(cens=cens,
            t=t,
            N=N,
            M=M,
            stage=stage,
            rnd_group=rnd_group,
            progression_num=progression_num,
            Asbestosis=Asbestosis,
            Sarcoidosis=Sarcoidosis,
            gender=gender,
            work_with_smoker=work_with_smoker,
            race_asian=race_asian,
            lym_lym=s_theta_lym_lym,
            lym_str=s_theta_lym_str,
            lym_tum=s_theta_lym_tum,
            str_str=s_theta_str_str,
            str_tum=s_theta_str_tum,
            omega_lym=s_omega_lym,
            omega_str=s_omega_str,
            lambda=lambda
            
            
)

inits = function() {list(beta.intercept=0,
                         beta.stage2=0,
                         beta.stage3=0,
                         beta.stage4=0,
                         beta.rnd=0,
                         beta.sarcoid=0,
                         beta.gender=0,
                         beta.race_asian=0,
                         beta.progr=0,
                         beta.work_smoke=0,
                         beta.asbest=0,
                         beta.omega_str=0,
                         beta.omega_lym=0, 
                         beta.lym_lym=0, 
                         beta.lym_str=0, 
                         beta.lym_tum=0,
                         beta.str_str=0, 
                         beta.str_tum=0, 
                         beta.lambda=0,
                         tau=0.3,
                         tau2=1,
                         .RNG.seed = 418, .RNG.name = 'base::Wichmann-Hill') 
}


modelAFT <- jags.model("surv_random_eff_weib.bug",data=dat,inits=inits,n.adapt=10000,n.chains=2)
update(modelAFT,50000)

  variable.names=c("beta.intercept",
                 "beta.stage2",
                 "beta.stage3",
                 "beta.stage4",
                 "beta.rnd",
                 "beta.progr",
                 "beta.sarcoid",
                 "beta.asbest",
                 "beta.work_smoke",
                 "beta.omega_str",
                 "beta.omega_lym", 
                 "beta.lym_lym", 
                 "beta.lym_str", 
                 "beta.lym_tum",
                 "beta.str_str", 
                 "beta.str_tum", 
                 "beta.lambda",
                 "sigma"
                 ) # monitoring
#



n.iter=50000
thin=10

outputAFT <- coda.samples(model=modelAFT,variable.names=variable.names,n.iter=n.iter,thin=thin)
save(outputAFT,file='AFT_Weibull_output_rand.Rdata')   # we save the chain
load('AFT_Weibull_output_rand.Rdata' )


data.out <- as.matrix(outputAFT)
data.out <- data.frame(data.out)

beta.intercept<-as.data.frame(data.out$beta.intercept)
beta.intercept$cathegory<-"beta.intercept"
names(beta.intercept)[names(beta.intercept) == "data.out$beta.intercept"] <- "value"


beta.asbest<-as.data.frame(data.out$beta.asbest)
beta.asbest$cathegory<-"beta.asbest"
names(beta.asbest)[names(beta.asbest) == "data.out$beta.asbest"] <- "value"

beta.lambda<-as.data.frame(data.out$beta.lambda)
beta.lambda$cathegory<-"beta.lambda"
names(beta.lambda)[names(beta.lambda) == "data.out$beta.lambda"] <- "value"

beta.lym_lym<-as.data.frame(data.out$beta.lym_lym)
beta.lym_lym$cathegory<-"beta.lym_lym"
names(beta.lym_lym)[names(beta.lym_lym) == "data.out$beta.lym_lym"] <- "value"

beta.lym_str<-as.data.frame(data.out$beta.lym_str)
beta.lym_str$cathegory<-"beta.lym_str"
names(beta.lym_str)[names(beta.lym_str) == "data.out$beta.lym_str"] <- "value"

beta.lym_tum<-as.data.frame(data.out$beta.lym_tum)
beta.lym_tum$cathegory<-"beta.lym_tum"
names(beta.lym_tum)[names(beta.lym_tum) == "data.out$beta.lym_tum"] <- "value"


beta.omega_lym<-as.data.frame(data.out$beta.omega_lym)
beta.omega_lym$cathegory<-"beta.omega_lym"
names(beta.omega_lym)[names(beta.omega_lym) == "data.out$beta.omega_lym"] <- "value"

beta.omega_str<-as.data.frame(data.out$beta.omega_str)
beta.omega_str$cathegory<-"beta.omega_str"
names(beta.omega_str)[names(beta.omega_str) == "data.out$beta.omega_str"] <- "value"

beta.progr<-as.data.frame(data.out$beta.progr)
beta.progr$cathegory<-"beta.progr"
names(beta.progr)[names(beta.progr) == "data.out$beta.progr"] <- "value"

beta.rnd<-as.data.frame(data.out$beta.rnd)
beta.rnd$cathegory<-"beta.rnd"
names(beta.rnd)[names(beta.rnd) == "data.out$beta.rnd"] <- "value"

beta.sarcoid<-as.data.frame(data.out$beta.sarcoid)
beta.sarcoid$cathegory<-"beta.sarcoid"
names(beta.sarcoid)[names(beta.sarcoid) == "data.out$beta.sarcoid"] <- "value"

beta.stage2<-as.data.frame(data.out$beta.stage2)
beta.stage2$cathegory<-"beta.stage2"
names(beta.stage2)[names(beta.stage2) == "data.out$beta.stage2"] <- "value"

beta.stage3<-as.data.frame(data.out$beta.stage3)
beta.stage3$cathegory<-"beta.stage3"
names(beta.stage3)[names(beta.stage3) == "data.out$beta.stage3"] <- "value"

beta.stage4<-as.data.frame(data.out$beta.stage4)
beta.stage4$cathegory<-"beta.stage4"
names(beta.stage4)[names(beta.stage4) == "data.out$beta.stage4"] <- "value"

beta.str_str<-as.data.frame(data.out$beta.str_str)
beta.str_str$cathegory<-"beta.str_str"
names(beta.str_str)[names(beta.str_str) == "data.out$beta.str_str"] <- "value"

beta.str_tum<-as.data.frame(data.out$beta.str_tum)
beta.str_tum$cathegory<-"beta.str_tum"
names(beta.str_tum)[names(beta.str_tum) == "data.out$beta.str_tum"] <- "value"


beta.work_smoke<-as.data.frame(data.out$beta.work_smoke)
beta.work_smoke$cathegory<-"beta.work_smoke"
names(beta.work_smoke)[names(beta.work_smoke) == "data.out$beta.work_smoke"] <- "value"

sigma<-as.data.frame(data.out$sigma)
sigma$cathegory<-"sigma"
names(sigma)[names(sigma) == "data.out$sigma"] <- "value"

data_for_plot<-rbind(beta.asbest,
                     beta.intercept,
                     beta.lambda,
                     beta.lym_lym,
                     beta.lym_str,
                     beta.lym_tum,
                     beta.omega_lym,
                     beta.omega_str,
                     beta.progr,
                     beta.rnd,
                     beta.stage2,
                     beta.stage3,
                     beta.stage4,
                     beta.str_str,
                     beta.str_tum,
                     beta.work_smoke,
                     sigma)
library(ggridges)

p<-ggplot(data_for_plot,aes(x = value, y = cathegory, fill = cathegory)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
x11()
p









library(coda)

x11()
plot(outputAFT,ask=T)
dev.off()


data.out <- as.matrix(outputAFT)
data.out <- data.frame(data.out)
attach(data.out)
n.chain <- dim(data.out)[1]
n.chain

summary(data.out)


my_WAIC <- function(fit, param){
  
  llik   <- rstan::extract(fit, param)[[1]]
  p_WAIC <- sum(apply(llik, 2, var))
  lppd   <- sum(apply(llik, 2, function(x) log(mean(exp(x)))))
  WAIC   <- - 2 * lppd + 2 * p_WAIC
  return(WAIC)
  
}


my_WAIC(model, "log_lik")
my_WAIC(model_fix, "log_lik")




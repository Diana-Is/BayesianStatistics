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


setwd("C:/Users/Diana/Desktop/Bayes_project/")
data <- read.csv("merged_clin_params.csv", header=T)

data[c("X","lym","tum","str","lym_lym","lym_str","lym_tum","str_str","str_tum","tum_tum")]<-list(NULL)

myvars <- c("fup_days","death_days","diag_date","pid","age.x","cigsmok","gender","de_stag_7thed","lym.given.lym",
            "str.given.lym","tum.given.lym","lym.given.str",
            "str.given.str","tum.given.str","lym.given.tum",
            "str.given.tum","tum.given.tum","pi_lym",
            "pi_str","pi_tum","lambda")

data <- data[myvars]

names(data)[names(data) == "age.x"] <- "age"
names(data)[names(data) == "cigsmok"] <- "smoke_stat"
names(data)[names(data) == "de_stag_7thed"] <- "stage"


data$gender[which(data$gender==1)]<-0 #male
data$gender[which(data$gender==2)]<-1 #female

pids<-unique(data$pid)
filter5<-data[which(data$pid==pids[1]),]
filter5<-filter5[1:5,]

for(i in 2:length(pids)) {
  temp<-data[which(data$pid==pids[i]),]
  if (dim(temp)[1]>5) temp<-temp[1:5,]
  filter5<-rbind(filter5,temp)
  
}

#-------------------------------------------------------------------
#Standartization of non-categorical variables
#-------------------------------------------------------------------

t<-data$death_days-data$diag_date
cens<-data$fup_days-data$diag_date
tmax <- round(max(t[!is.na(t)]))+50
cens[!is.na(t)] <- tmax
#removing irrelevant variables
data[c("fup_days","death_days","diag_date")]<-list(NULL)

data$str.given.lym<-(data$str.given.lym-mean(data$str.given.lym))/sd(data$str.given.lym)
data$tum.given.lym<-(data$tum.given.lym-mean(data$tum.given.lym))/sd(data$tum.given.lym)
data$lym.given.str<-(data$lym.given.str-mean(data$lym.given.str))/sd(data$lym.given.str)

data$tum.given.str<-(data$tum.given.str-mean(data$tum.given.str))/sd(data$tum.given.str)
data$lym.given.tum<-(data$lym.given.tum-mean(data$lym.given.tum))/sd(data$lym.given.tum)
data$str.given.tum<-(data$str.given.tum-mean(data$str.given.tum))/sd(data$str.given.tum)


data$pi_str<-(data$pi_str-mean(data$pi_str))/sd(data$pi_str)
data$pi_lym<-(data$pi_lym-mean(data$pi_lym))/sd(data$pi_lym)
data$lambda<-(data$lambda-mean(data$lambda))/sd(data$lambda)


N <- dim(data)[1]/5
M<-5
t<-t(matrix(t,nrow = 5,ncol = 296))
cens<-t(matrix(cens,nrow = 5,ncol = 296))
s_str_lym<-t(matrix(data$str.given.lym,nrow = 5,ncol = 296))
s_tum_lym<-t(matrix(data$tum.given.lym,nrow = 5,ncol = 296))
s_lym_str<-t(matrix(data$lym.given.str,nrow = 5,ncol = 296))

s_tum_str<-t(matrix(data$tum.given.str,nrow = 5,ncol = 296))
s_lym_tum<-t(matrix(data$lym.given.tum,nrow = 5,ncol = 296))
s_str_tum<-t(matrix(data$str.given.tum,nrow = 5,ncol = 296))

s_pi_str<-t(matrix(data$pi_str,nrow = 5,ncol = 296))
s_pi_lym<-t(matrix(data$pi_lym,nrow = 5,ncol = 296))
lambda<-t(matrix(data$lambda,nrow = 5,ncol = 296))

smoke_stat<-t(matrix(data$smoke_stat,nrow = 5,ncol = 296))[,1]
age<-t(matrix(data$age,nrow = 5,ncol = 296))[,1]
age<-(age-mean(age))/sd(age)
gender<-t(matrix(data$gender,nrow = 5,ncol = 296))[,1]
stage<-t(matrix(data$stage,nrow = 5,ncol = 296))[,1]

dat <- list(cens=cens,
            t=t,
            N=N,
            M=M,
            stage=stage,
            gender=gender,
            age=age,
            smoke_stat=smoke_stat,
            str_lym=s_str_lym,
            tum_lym=s_tum_lym,
            lym_str=s_lym_str,
            tum_str=s_tum_str,
            lym_tum=s_lym_tum,
            str_tum=s_str_tum,
            pi_lym=s_pi_lym,
            pi_str=s_pi_str,
            lambda=lambda
)
#--------------------------------------------------
#HPD metrics output:
#--------------------------------------------------

inits = function() {list(beta.intercept=0,
                         beta.stage2=0,
                         beta.stage3=0,
                         beta.stage4=0,
                         beta.stage5=0,
                         beta.stage6=0,
                         beta.stage7=0,
                         beta.gender=0,
                         beta.age=0,
                         beta.smoke_stat=0,
                         beta.pi_str=0,
                         beta.pi_lym=0, 
                         beta.strlym=0, 
                         beta.tumlym=0, 
                         beta.lymstr=0,
                         beta.tumstr=0, 
                         beta.lymtum=0, 
                         beta.strtum=0,
                         beta.lambda=0,
                         tau=0.3,
                         tau2=1,
                         .RNG.seed = 48, .RNG.name = 'base::Wichmann-Hill') 
}

modelAFT <- jags.model("surv_random_base_weib.bug",data=dat,inits=inits,n.adapt=5000,n.chains=2)

update(modelAFT,30000)

variable.names=c("beta.intercept",
                 "beta.stage2",
                 "beta.stage3",
                 "beta.stage4",
                 "beta.stage5",
                 "beta.stage6",
                 "beta.stage7",
                 "beta.gender",
                 "beta.age",
                 "beta.smoke_stat",
                 "beta.pi_str",
                 "beta.pi_lym", 
                 "beta.strlym", 
                 "beta.tumlym", 
                 "beta.lymstr",
                 "beta.tumstr", 
                 "beta.lymtum", 
                 "beta.strtum",
                 "beta.lambda",
                 "sigma"
) # monitoring
#

n.iter=30000
thin=10

outputAFT <- coda.samples(model=modelAFT,variable.names=variable.names,n.iter=n.iter,thin=thin)
save(outputAFT,file='AFT_output_rand_base.Rdata')   # we save the chain

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

x11()
par(mfrow=c(2,3))
acf(data.out[,'beta.intercept'],lwd=3,col="red3",main="autocorrelation of beta.intercept")
acf(data.out[,'beta.stage2'],lwd=3,col="red3",main="autocorrelation of beta.stage2")
acf(data.out[,'beta.stage3'],lwd=3,col="red3",main="autocorrelation of beta.stage3")
acf(data.out[,'beta.stage4'],lwd=3,col="red3",main="autocorrelation of beta.stage4")
acf(data.out[,'beta.stage5'],lwd=3,col="red3",main="autocorrelation of beta.stage5")
acf(data.out[,'beta.stage6'],lwd=3,col="red3",main="autocorrelation of beta.stage6")

x11()
par(mfrow=c(2,3))
acf(data.out[,'beta.stage7'],lwd=3,col="red3",main="autocorrelation of beta.stage7")
acf(data.out[,'beta.gender'],lwd=3,col="red3",main="autocorrelation of beta.gender")
acf(data.out[,'beta.rnd'],lwd=3,col="red3",main="autocorrelation of beta.rnd")


x11()
par(mfrow=c(2,3))
plot(ts(data.out[,'beta.1.']),xlab="t",ylab="beta1")
plot(ts(data.out[,'beta.2.']),xlab="t",ylab="beta2")
plot(ts(data.out[,'beta.3.']),xlab="t",ylab="beta3")
plot(ts(data.out[,'beta.4.']),xlab="t",ylab="beta4")
plot(ts(data.out[,'beta.5.']),xlab="t",ylab="beta5")
plot(ts(data.out[,'beta.6.']),xlab="t",ylab="beta6")


x11()
par(mfrow=c(2,3))
plot(density(data.out[,'beta.intercept'],adj=2),  xlab=expression("beta.intercept"),main="")
abline(v=quantile(data.out[,'beta.intercept'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.stage2'],adj=2),  xlab=expression("beta.stage2"),main="")
abline(v=quantile(data.out[,'beta.stage2'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.stage3'],adj=2),  xlab=expression("beta.stage3"),main="")
abline(v=quantile(data.out[,'beta.stage3'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.stage4'],adj=2),  xlab=expression("beta.stage4"),main="")
abline(v=quantile(data.out[,'beta.stage4'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.stage5'],adj=2),  xlab=expression("beta.stage5"),main="")
abline(v=quantile(data.out[,'beta.stage5'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.stage6'],adj=2),  xlab=expression("beta.stage6"),main="")
abline(v=quantile(data.out[,'beta.stage6'],prob=c(.025,.975)),lwd=2,col="red")

x11()
par(mfrow=c(2,3))
plot(density(data.out[,'beta.stage7'],adj=2),  xlab=expression("beta.stage7"),main="")
abline(v=quantile(data.out[,'beta.stage7'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.gender'],adj=2),  xlab=expression("beta.gender"),main="")
abline(v=quantile(data.out[,'beta.gender'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.rnd'],adj=2),  xlab=expression("beta.rnd"),main="")
abline(v=quantile(data.out[,'beta.rnd'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.lesion'],adj=2),  xlab=expression("beta.lesion"),main="")
abline(v=quantile(data.out[,'beta.lesion'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.pipe'],adj=2),  xlab=expression("beta.pipe"),main="")
abline(v=quantile(data.out[,'beta.pipe'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.worksmoke'],adj=2),  xlab=expression("beta.worksmoke"),main="")
abline(v=quantile(data.out[,'beta.worksmoke'],prob=c(.025,.975)),lwd=2,col="red")


x11()
par(mfrow=c(2,3))
plot(density(data.out[,'beta.progression'],adj=2),  xlab=expression("beta.progression"),main="")
abline(v=quantile(data.out[,'beta.progression'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.lung_dang_work'],adj=2),  xlab=expression("beta.lung_dang_work"),main="")
abline(v=quantile(data.out[,'beta.lung_dang_work'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.asbest'],adj=2),  xlab=expression("beta.asbest"),main="")
abline(v=quantile(data.out[,'beta.asbest'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.bronch'],adj=2),  xlab=expression("beta.bronch"),main="")
abline(v=quantile(data.out[,'beta.bronch'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.diabet'],adj=2),  xlab=expression("beta.diabet"),main="")
abline(v=quantile(data.out[,'beta.diabet'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.emphys'],adj=2),  xlab=expression("beta.emphys"),main="")
abline(v=quantile(data.out[,'beta.emphys'],prob=c(.025,.975)),lwd=2,col="red")


x11()
par(mfrow=c(2,3))
plot(density(data.out[,'beta.heart'],adj=2),  xlab=expression("beta.heart"),main="")
abline(v=quantile(data.out[,'beta.heart'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.sarcoid'],adj=2),  xlab=expression("beta.sarcoid"),main="")
abline(v=quantile(data.out[,'beta.sarcoid'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.tub'],adj=2),  xlab=expression("beta.tub"),main="")
abline(v=quantile(data.out[,'beta.tub'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.stroke'],adj=2),  xlab=expression("beta.stroke"),main="")
abline(v=quantile(data.out[,'beta.stroke'],prob=c(.025,.975)),lwd=2,col="red")




my_WAIC <- function(fit, param){
  
  llik   <- rstan::extract(fit, param)[[1]]
  p_WAIC <- sum(apply(llik, 2, var))
  lppd   <- sum(apply(llik, 2, function(x) log(mean(exp(x)))))
  WAIC   <- - 2 * lppd + 2 * p_WAIC
  return(WAIC)
  
}


my_WAIC(model, "log_lik")
my_WAIC(model_fix, "log_lik")




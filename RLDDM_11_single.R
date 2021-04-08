

#rm(list=ls())  # remove all variables 

library(rstan)
library(dplyr)

# set seed
set.seed(08826)

# source HDIofMCMC.R to calculate HDI
source("HDIofMCMC.R") 

# read the data file
dat = read.table("TST_nspn.txt", header=T, sep="\t")
dat1 = filter(dat, measurment==1)
dat2 = filter(dat, measurment==2)
#generate 50 subjects from 543x
s <- c(186, 222,  67, 480,  51, 254, 401, 422, 105, 209, 538, 492, 110, 206, 329, 120,  61, 456,  9,  22, 112, 438, 399, 549, 264, 292, 548, 526, 156,
       395, 26, 159, 350, 117, 345, 114, 267, 311, 195,  59, 140,  87, 447, 542, 512, 139, 106, 331,  63, 337)
length(s)
dat2 = filter(dat2, subject %in% s)
dat2

allSubjs = unique(dat2$subject)  # all subject IDs

N = length(allSubjs)
N # number of subjects

T = max(table(dat2$subject))  # number of trials per subject
T

###maybe we should see what T is; is it different for each subject? If that's the case
Tsubj = as.numeric(as.list(table(dat2$subject)))
numIter = 100             # number of iterations to find global minimum values
numPars = 3               # number of parameters

t = split(dat2$transition, dat2$subject, drop = FALSE)
transition = t[1]
for (i in 2:N){
  transition = rbind(transition,t[i])
}
transition = lapply(transition,as.numeric)
transition = lapply(transition, 'length<-', T)
transition = lapply(transition, function(x) replace(x, is.na(x), 0))

s = split((dat2$second_stage_state+1), dat2$subject, drop = FALSE)
secstagestate = s[1]
for (i in 2:N){ # skip loop if single subject
  secstagestate = rbind(secstagestate,s[i])
}
secstagestate = lapply(secstagestate,as.numeric)
secstagestate = lapply(secstagestate, 'length<-', T)
secstagestate = lapply(secstagestate, function(x) replace(x, is.na(x), 2))

c1 = split(dat2$choice1, dat2$subject, drop = FALSE)
choice1 = c1[1]
for (i in 2:N){ # skip loop if single subject
  choice1 = rbind(choice1,c1[i])
}
choice1 = lapply(choice1,as.numeric)
choice1 = lapply(choice1, 'length<-', T)
choice1 = lapply(choice1, function(x) replace(x, is.na(x), 1))

c2 = split(dat2$choice2, dat2$subject, drop = FALSE)
choice2 = c2[1]
for (i in 2:N){ # skip loop if single subject
  choice2 = rbind(choice2,c2[i])
}
choice2 = lapply(choice2,as.numeric)
choice2 = lapply(choice2, 'length<-', T)
choice2 = lapply(choice2, function(x) replace(x, is.na(x), 1))

rt1 = split(dat2$RT1, dat2$subject, drop = FALSE)
RT1 = rt1[1]
for (i in 2:N){ # skip loop if single subject
  RT1 = rbind(RT1,rt1[i])
}
RT1 = lapply(RT1,as.numeric)
RT1 = lapply(RT1, 'length<-', T)
RT1 = lapply(RT1, function(x) replace(x, is.na(x), 0))

rt2 = split(dat2$RT2, dat2$subject, drop = FALSE)
RT2 = rt2[1]
for (i in 2:N){ # skip loop if single subject
  RT2 = rbind(RT2,rt2[i])
}
RT2 = lapply(RT2,as.numeric)
RT2 = lapply(RT2, 'length<-', T)
RT2 = lapply(RT2, function(x) replace(x, is.na(x), 0))

r = split(dat2$reward, dat2$subject, drop = FALSE)
reward = r[1]
for (i in 2:N){ # skip loop if single subject
  reward = rbind(reward,r[i])
}
reward = lapply(reward,as.numeric)
reward = lapply(reward, 'length<-', T)
reward = lapply(reward, function(x) replace(x, is.na(x), 0))

rt1 = split(dat2$RT1, dat2$subject, drop = FALSE)
minRT1 = min(rt1[[1]])
for (i in 2:N){ # skip loop if single subject
  minRT1 = rbind(minRT1,min(rt1[[i]]))
}
minRT1 = lapply(minRT1,as.numeric)
minRT1 = unlist(minRT1, use.names=FALSE)

rt2 = split(dat2$RT2, dat2$subject, drop = FALSE)
minRT2 = min(rt2[[1]])
for (i in 2:N){ # skip loop if single subject
  minRT2 = rbind(minRT2,min(rt2[[i]]))
}
minRT2 = lapply(minRT2,as.numeric)
minRT2 = unlist(minRT2, use.names=FALSE)

RTbound1 = rep(0.1,N)
RTbound2 = rep(0.1,N)


# use all data
dataList3 <- list(
  N                = N,
  T                = T,
  Tsubj            = Tsubj,
  transition       = transition,
  secstagestate    = secstagestate,
  choice1          = choice1,
  choice2          = choice2,
  RT1              = RT1,
  RT2              = RT2,
  reward           = reward,
  minRT1           = minRT1,
  minRT2           = minRT2,
  RTbound1         = RTbound1,
  RTbound2         = RTbound2
)

# run!
# model3 = stan("RLDDM_11_single.stan", data = dataList3, pars = c("alpha1", "alpha2", "p", "omega", "lambda", 
#                                                                   "a1", "a2", "b1", "b2", "tau1", "tau2",
#                                                                   "rt1","rt2","ch1","ch2"),
#                 iter = 4000, warmup=2000, chains=4, cores=4, control = list(adapt_delta = 0.99))
# 
# save.image(file='model3.RData')
# traceplot & posterior histogram - put in loop
load('model3.RData')
for (i in 1:N){
  png(file=paste("model3_",i,"_traceplot.png",sep=''),width=600, height=350)
  t<-traceplot(model3, pars=c(paste('alpha1[',i,']',sep=''),paste('alpha2[',i,']',sep=''),paste('p[',i,']',sep=''),paste('omega[',i,']',sep=''),
                              paste('lambda[',i,']',sep=''),paste('a1[',i,']',sep=''),paste('a2[',i,']',sep=''),
                              paste('b1[',i,']',sep=''),paste('b2[',i,']',sep=''),paste('tau1[',i,']',sep=''),paste('tau2[',i,']',sep='')), inc_warmup=FALSE)
  print(t)
  dev.off()
  
}

# print summary
print(model3, pars=c("alpha1"), include=TRUE)
print(model3, pars=c("alpha2"), include=TRUE)
print(model3, pars=c("p"), include=TRUE)
print(model3, pars=c("lambda"), include=TRUE)
print(model3, pars=c("omega"), include=TRUE)
print(model3, pars=c("a1"), include=TRUE)
print(model3, pars=c("a2"), include=TRUE)
print(model3, pars=c("b1"), include=TRUE)
print(model3, pars=c("b2"), include=TRUE)
print(model3, pars=c("tau1"), include=TRUE)
print(model3, pars=c("tau2"), include=TRUE)


# extract Stan fit object (parameters)
parameters <- rstan::extract(model3)

library(bayesplot)

png(file="model3_alpha1_posterior.png",width=2000, height=1200)
plot(model3, plotfun = "hist", pars = c("alpha1"), include = TRUE)
dev.off()
png(file="model3_alpha2_posterior.png",width=2000, height=1200)
plot(model3, plotfun = "hist", pars = c("alpha2"), include = TRUE)
dev.off()
png(file="model3_p_posterior.png",width=2000, height=1200)
plot(model3, plotfun = "hist", pars = c("p"), include = TRUE)
dev.off()
png(file="model3_omega_posterior.png",width=2000, height=1200)
plot(model3, plotfun = "hist", pars = c("omega"), include = TRUE)
dev.off()
png(file="model3_lambda_posterior.png",width=2000, height=1200)
plot(model3, plotfun = "hist", pars = c("lambda"), include = TRUE)
dev.off()
png(file="model3_a1_posterior.png",width=2000, height=1200)
plot(model3, plotfun = "hist", pars = c('a1'), include = TRUE)
dev.off()
png(file="model3_a2_posterior.png",width=2000, height=1200)
plot(model3, plotfun = "hist", pars = c('a2'), include = TRUE)
dev.off()
png(file="model3_b1_posterior.png",width=2000, height=1200)
plot(model3, plotfun = "hist", pars = c('b1'), include = TRUE)
dev.off()
png(file="model3_b2_posterior.png",width=2000, height=1200)
plot(model3, plotfun = "hist", pars = c('b2'), include = TRUE)
dev.off()
png(file="model3_tau1_posterior.png",width=2000, height=1200)
plot(model3, plotfun = "hist", pars = c('tau1'), include = TRUE)
dev.off()
png(file="model3_tau2_posterior.png",width=2000, height=1200)
plot(model3, plotfun = "hist", pars = c('tau2'), include = TRUE)
dev.off()

# 95% HDI of rho
for (i in 1:N){
  assign(paste('alpha1_95',i,sep=''),HDIofMCMC(parameters$alpha1[,i], credMass = 0.95))
  assign(paste('alpha2_95',i,sep=''),HDIofMCMC(parameters$alpha2[,i], credMass = 0.95))
  assign(paste('p_95',i,sep=''),HDIofMCMC(parameters$p[,i], credMass = 0.95))
  assign(paste('omega_95',i,sep=''),HDIofMCMC(parameters$omega[,i], credMass = 0.95))
  assign(paste('lambda_95',i,sep=''),HDIofMCMC(parameters$lambda[,i], credMass = 0.95))
  assign(paste('a1_95',i,sep=''),HDIofMCMC(parameters$a1[,i], credMass = 0.95))
  assign(paste('a2_95',i,sep=''),HDIofMCMC(parameters$a2[,i], credMass = 0.95))
  assign(paste('b1_95',i,sep=''),HDIofMCMC(parameters$b1[,i], credMass = 0.95))
  assign(paste('b2_95',i,sep=''),HDIofMCMC(parameters$b2[,i], credMass = 0.95))
  assign(paste('tau1_95',i,sep=''),HDIofMCMC(parameters$tau1[,i], credMass = 0.95))
  assign(paste('tau2_95',i,sep=''),HDIofMCMC(parameters$tau2[,i], credMass = 0.95))
}











####PPC####

alpha1 = parameters$alpha1
alpha1_mean = rep(NA, length=N)

alpha2 = parameters$alpha2
alpha2_mean = rep(NA, length=N)

p = parameters$p
p_mean = rep(NA, length=N)

omega = parameters$omega
omega_mean = rep(NA, length=N)

lambda = parameters$lambda
lambda_mean = rep(NA, length=N)

a1 = parameters$a1
a1_mean = rep(NA, length=N)

a2 = parameters$a2
a2_mean = rep(NA, length=N)

b1 = parameters$b1
b1_mean = rep(NA, length=N)

b2 = parameters$b2
b2_mean = rep(NA, length=N)

tau1 = parameters$tau1
tau1_mean = rep(NA, length=N)

tau2 = parameters$tau2
tau2_mean = rep(NA, length=N)

for (i in 1:N){
  alpha1_mean[i] = mean(alpha1[i,])
  alpha2_mean[i] = mean(alpha2[i,])
  p_mean[i] = mean(p[i,])
  omega_mean[i] = mean(omega[i,])
  lambda_mean[i] = mean(lambda[i,])
  a1_mean[i] = mean(a1[i,])
  a2_mean[i] = mean(a2[i,])
  b1_mean[i] = mean(b1[i,])
  b2_mean[i] = mean(b2[i,])
  tau1_mean[i] = mean(tau1[i,])
  tau2_mean[i] = mean(tau2[i,])
}

#plot posterior predictive check
library(bayesplot)

#for each subject

ch1_rep <- array(parameters$ch1,c(T,N,2000))

ch2_rep <- array(parameters$ch2,c(T,N,2000))

rt1_rep <- array(parameters$rt1,c(T,N,2000))

rt2_rep <- array(parameters$rt2,c(T,N,2000))

for (i in 1:N){
  png(file=paste("RLDDM_11_single_ch1_subject",i,"__ppc.png",sep=''),width=600, height=350)
  ch1<-ppc_bars(unlist(choice1[i]), matrix(ch1_rep[,i,],ncol=T), freq=FALSE)
  print(ch1)
  dev.off()
  png(file=paste("RLDDM_11_single_ch2_subject",i,"__ppc.png",sep=''),width=600, height=350)
  ch2<-ppc_bars(unlist(choice2[i]), matrix(ch2_rep[,i,],ncol=T), freq=FALSE)
  print(ch2)
  dev.off()
  png(file=paste("RLDDM_11_single_rt1_subject",i,"__ppc.png",sep=''),width=600, height=350)
  rt1<-ppc_stat(unlist(RT1[i]), matrix(rt1_rep[,i,],ncol=T)) 
  print(rt1)
  dev.off()
  png(file=paste("RLDDM_11_single_rt2_subject",i,"__ppc.png",sep=''),width=600, height=350)
  rt2<-ppc_stat(unlist(RT2[i]), matrix(rt2_rep[,i,],ncol=T))
  print(rt2)
  dev.off()
}












####parameter recovery####

##for simulating 1000 trials

dataList1_sim <- list(
  N                = N,
  T                = 1000,
  Tsubj            = rep(1000,N),
  alpha1           = alpha1_mean,
  alpha2           = alpha2_mean,
  p                = p_mean,
  omega            = omega_mean,
  lambda           = lambda_mean,
  a1               = a1_mean,
  a2               = a2_mean,
  b1               = b1_mean,
  b2               = b2_mean,
  tau1             = tau1_mean,
  tau2             = tau2_mean
)

model3_gen_1000 = stan("RLDDM_11_sim.stan", data = dataList1_sim, pars = c("rt1","rt2","ch1","ch2","arwrd","asecstagestate","trn"),
                       iter = 1, chains=1, cores=1, algorithm="Fixed_param")

save.image(file='model3.RData')
print(model3_gen_1000)

parameters_gen_1000 <- rstan::extract(model3_gen_1000)


#for each subject

ch1_rep_1000 <- parameters_gen_1000$ch1

ch2_rep_1000 <- parameters_gen_1000$ch2

rt1_rep_1000 <- parameters_gen_1000$rt1

rt2_rep_1000 <- parameters_gen_1000$rt2

#preprocess minRT data from generated data(from parameters, not generated data since minRT doesn't work well with that)

minRT1r = min(c(rt1_rep[,1,]))
for (i in 2:N){ # skip loop if single subject
  minRT1r = c(minRT1r,min(c(rt1_rep[,i,])))
}
minRT1r = lapply(minRT1r,as.numeric)
minRT1r = unlist(minRT1r, use.names=FALSE)
minRT1r

minRT2r = min(c(rt2_rep[,1,]))
for (i in 2:N){ # skip loop if single subject
  minRT2r = c(minRT2r,min(c(rt2_rep[,i,])))
}
minRT2r = lapply(minRT2r,as.numeric)
minRT2r = unlist(minRT2r, use.names=FALSE)
RTbound = rep(0.1,N)

##for 1000 simulated trials (try larger trials)

dataList1_recovery <- list(
  N                = N,
  T                = 1000,
  Tsubj            = rep(1000,N),
  transition       = parameters_gen_1000$trn[1,,],
  secstagestate    = parameters_gen_1000$asecstagestate[1,,]+1,
  choice1          = ch1_rep_1000[1,,],
  choice2          = ch2_rep_1000[1,,],
  RT1              = rt1_rep_1000[1,,],
  RT2              = rt2_rep_1000[1,,],
  reward           = parameters_gen_1000$arwrd[1,,],
  minRT1           = minRT1r,
  minRT2           = minRT2r,
  RTbound1         = RTbound,
  RTbound2         = RTbound
)

# Let's fit the simulated data
model3_recovery_1000_RL = stan("RLDDM_11_single.stan", data = dataList1_recovery, pars = c("alpha1", "alpha2", "p", "omega", "lambda"),
                       iter = 4000, warmup=2000, chains=4, cores=4)

model3_recovery_1000_DDM = stan("RLDDM_11_single.stan", data = dataList1_recovery, pars = c("a1", "a2", "b1", "b2", "tau1", "tau2"),
                       iter = 4000, warmup=2000, chains=4, cores=4)

parameters_recovery <- rstan::extract(model3_recovery)

save.image(file='model3.RData')

#plot recovery posteriors if you want..


#plot recovery plots

alpha1r = parameters_recovery$alpha1
alpha1_recovered_mean = rep(NA, length=N)

alpha2r = parameters_recovery$alpha2
alpha2_recovered_mean = rep(NA, length=N)

pr = parameters_recovery$p
p_recovered_mean = rep(NA, length=N)

omegar = parameters_recovery$omega
omega_recovered_mean = rep(NA, length=N)

lambdar = parameters_recovery$lambda
lambda_recovered_mean = rep(NA, length=N)

a1r = parameters_recovery$a1
a1_recovered_mean = rep(NA, length=N)

a2r = parameters_recovery$a2
a2_recovered_mean = rep(NA, length=N)

b1r = parameters_recovery$b1
b1_recovered_mean = rep(NA, length=N)

b2r = parameters_recovery$b2
b2_recovered_mean = rep(NA, length=N)

tau1r = parameters_recovery$tau1
tau1_recovered_mean = rep(NA, length=N)

tau2r = parameters_recovery$tau2
tau2_recovered_mean = rep(NA, length=N)

for (i in 1:N){
  alpha1_recovered_mean[i] = mean(alpha1r[i,])
  alpha2_recovered_mean[i] = mean(alpha2r[i,])
  p_recovered_mean[i] = mean(pr[i,])
  omega_recovered_mean[i] = mean(omegar[i,])
  lambda_recovered_mean[i] = mean(lambdar[i,])
  a1_recovered_mean[i] = mean(a1r[i,])
  a2_recovered_mean[i] = mean(a2r[i,])
  b1_recovered_mean[i] = mean(b1r[i,])
  b2_recovered_mean[i] = mean(b2r[i,])
  tau1_recovered_mean[i] = mean(tau1r[i,])
  tau2_recovered_mean[i] = mean(tau2r[i,])
}

par(mfrow=c(4,3))

library(ggplot2)

everyone = seq(1:N)
allparams = data.frame(subject = everyone, alpha1 = alpha1_mean, alpha2 = alpha2_mean, alpha1_r = alpha1_recovered_mean, alpha2_r = alpha2_recovered_mean,
                       p = p_mean, p_r = p_recovered_mean,
                       omega = omega_mean, omega_r = omega_recovered_mean, lambda = lambda_mean, lambda_r = lambda_recovered_mean,
                       a1 = a1_mean, a1_r = a1_recovered_mean, a2 = a2_mean, a2_r = a2_recovered_mean, b1 = b1_mean, b1_r = b1_recovered_mean, b2 = b2_mean, b2_r = b2_recovered_mean,
                       tau1 = tau1_mean, tau1_r = tau1_recovered_mean, tau2 = tau2_mean, tau2_r = tau2_recovered_mean)

alpha1_corr<-round(cor(allparams$alpha1, allparams$alpha1_r),digits=3)
alpha1_corr
png(file="RLDDM_11_single_alpha1_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=alpha1, y=alpha1_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('alpha1_recovery_1000') +
  geom_text(x=median(alpha1_mean),y=max(alpha1_recovered_mean),
            label = paste("r=",alpha1_corr))
dev.off()

alpha2_corr<-round(cor(allparams$alpha2, allparams$alpha2_r),digits=3)
alpha2_corr
png(file="RLDDM_11_single_alpha2_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=alpha2, y=alpha2_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('alpha2_recovery_1000') +
  geom_text(x=median(alpha2_mean),y=max(alpha2_recovered_mean),
            label = paste("r=",alpha2_corr))
dev.off()

p_corr<-round(cor(allparams$p, allparams$p_r),digits=3)
p_corr
png(file="RLDDM_11_single_p_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=p, y=p_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('p_recovery_1000') +
  geom_text(x=median(p_mean),y=max(p_recovered_mean),
            label = paste("r=",p_corr))
dev.off()

omega_corr<-round(cor(allparams$omega, allparams$omega_r),digits=3)
omega_corr
png(file="RLDDM_11_single_omega_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=omega, y=omega_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('omega_recovery_1000') +
  geom_text(x=median(omega_mean),y=max(omega_recovered_mean),
            label = paste("r=",omega_corr))
dev.off()

lambda_corr<-round(cor(allparams$lambda, allparams$lambda_r),digits=3)
lambda_corr
png(file="RLDDM_11_single_lambda_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=lambda, y=lambda_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('lambda_recovery_1000') +
  geom_text(x=median(lambda_mean),y=max(lambda_recovered_mean),
            label = paste("r=",lambda_corr))
dev.off()

a1_corr<-round(cor(allparams$a1, allparams$a1_r),digits=3)
a1_corr
png(file="RLDDM_11_single_a1_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=a1, y=a1_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('a1_recovery_1000') +
  geom_text(x=median(a1_mean),y=max(a1_recovered_mean),
            label = paste("r=",a1_corr))
dev.off()

a2_corr<-round(cor(allparams$a2, allparams$a2_r),digits=3)
a2_corr
png(file="RLDDM_11_single_a2_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=a2, y=a2_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('a2_recovery_1000') +
  geom_text(x=median(a2_mean),y=max(a2_recovered_mean),
            label = paste("r=",a2_corr))
dev.off()

b1_corr<-round(cor(allparams$b1, allparams$b1_r),digits=3)
b1_corr
png(file="RLDDM_11_single_b1_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=b1, y=b1_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('b1_recovery_1000') +
  geom_text(x=median(b1_mean),y=max(b1_recovered_mean),
            label = paste("r=",b1_corr))
dev.off()

b2_corr<-round(cor(allparams$b2, allparams$b2_r),digits=3)
b2_corr
png(file="RLDDM_11_single_b2_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=b2, y=b2_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('b2_recovery_1000') +
  geom_text(x=median(b2_mean),y=max(b2_recovered_mean),
            label = paste("r=",b2_corr))
dev.off()

tau1_corr<-round(cor(allparams$tau1, allparams$tau1_r),digits=3)
tau1_corr
png(file="RLDDM_11_single_tau1_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=tau1, y=tau1_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('tau1_recovery_1000') +
  geom_text(x=median(tau1_mean),y=max(tau1_recovered_mean),
            label = paste("r=",tau1_corr))
dev.off()

tau2_corr<-round(cor(allparams$tau2, allparams$tau2_r),digits=3)
tau2_corr
png(file="RLDDM_11_single_tau2_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=tau2, y=tau2_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('tau2_recovery_1000') +
  geom_text(x=median(tau2_mean),y=max(tau2_recovered_mean),
            label = paste("r=",tau2_corr))
dev.off()








####parameter recovery####

##for simulating 200 trials

dataList1_sim <- list(
  N                = N,
  T                = 200,
  Tsubj            = rep(200,N),
  alpha1           = alpha1_mean,
  alpha2           = alpha2_mean,
  p                = p_mean,
  omega            = omega_mean,
  lambda           = lambda_mean,
  a1               = a1_mean,
  a2               = a2_mean,
  b1               = b1_mean,
  b2               = b2_mean,
  tau1             = tau1_mean,
  tau2             = tau2_mean
)

model3_gen_200 = stan("RLDDM_11_sim.stan", data = dataList1_sim, pars = c("rt1","rt2","ch1","ch2","arwrd","asecstagestate","trn"),
                      iter = 1, chains=1, cores=1, algorithm="Fixed_param")

save.image(file='model3.RData')
print(model3_gen_800)

parameters_gen_800 <- rstan::extract(model3_gen_800)


#for each subject

ch1_rep_200 <- parameters_gen_200$ch1

ch2_rep_200 <- parameters_gen_200$ch2

rt1_rep_200 <- parameters_gen_200$rt1

rt2_rep_200 <- parameters_gen_200$rt2

#preprocess minRT data from generated data(from parameters, not generated data since minRT doesn't work well with that)

minRT1r = min(c(rt1_rep[,1,]))
for (i in 2:N){ # skip loop if single subject
  minRT1r = c(minRT1r,min(c(rt1_rep[,i,])))
}
minRT1r = lapply(minRT1r,as.numeric)
minRT1r = unlist(minRT1r, use.names=FALSE)
minRT1r

minRT2r = min(c(rt2_rep[,1,]))
for (i in 2:N){ # skip loop if single subject
  minRT2r = c(minRT2r,min(c(rt2_rep[,i,])))
}
minRT2r = lapply(minRT2r,as.numeric)
minRT2r = unlist(minRT2r, use.names=FALSE)
RTbound = rep(0.1,N)

##for 800 simulated trials (try larger trials)

dataList1_recovery <- list(
  N                = N,
  T                = 200,
  Tsubj            = rep(200,N),
  transition       = parameters_gen_200$trn[1,,],
  secstagestate    = parameters_gen_200$asecstagestate[1,,]+1,
  choice1          = ch1_rep_200[1,,],
  choice2          = ch2_rep_200[1,,],
  RT1              = rt1_rep_200[1,,],
  RT2              = rt2_rep_200[1,,],
  reward           = parameters_gen_200$arwrd[1,,],
  minRT1           = minRT1r,
  minRT2           = minRT2r,
  RTbound1         = RTbound,
  RTbound2         = RTbound
)

# Let's fit the simulated data
model3_recovery_200 = stan("RLDDM_11_single.stan", data = dataList1_recovery, pars = c("alpha1", "alpha2", "p", "omega", "lambda",
                                                                                          "a1", "a2", "b1", "b2", "tau1", "tau2"),
                               iter = 4000, warmup=2000, chains=4, cores=4)

parameters_recovery <- rstan::extract(model3_recovery_200)

save.image(file='model3.RData')

#plot recovery posteriors if you want..


#plot recovery plots

alpha1r = parameters_recovery$alpha1
alpha1_recovered_mean = rep(NA, length=N)

alpha2r = parameters_recovery$alpha2
alpha2_recovered_mean = rep(NA, length=N)

pr = parameters_recovery$p
p_recovered_mean = rep(NA, length=N)

omegar = parameters_recovery$omega
omega_recovered_mean = rep(NA, length=N)

lambdar = parameters_recovery$lambda
lambda_recovered_mean = rep(NA, length=N)

a1r = parameters_recovery$a1
a1_recovered_mean = rep(NA, length=N)

a2r = parameters_recovery$a2
a2_recovered_mean = rep(NA, length=N)

b1r = parameters_recovery$b1
b1_recovered_mean = rep(NA, length=N)

b2r = parameters_recovery$b2
b2_recovered_mean = rep(NA, length=N)

tau1r = parameters_recovery$tau1
tau1_recovered_mean = rep(NA, length=N)

tau2r = parameters_recovery$tau2
tau2_recovered_mean = rep(NA, length=N)

for (i in 1:N){
  alpha1_recovered_mean[i] = mean(alpha1r[i,])
  alpha2_recovered_mean[i] = mean(alpha2r[i,])
  p_recovered_mean[i] = mean(pr[i,])
  omega_recovered_mean[i] = mean(omegar[i,])
  lambda_recovered_mean[i] = mean(lambdar[i,])
  a1_recovered_mean[i] = mean(a1r[i,])
  a2_recovered_mean[i] = mean(a2r[i,])
  b1_recovered_mean[i] = mean(b1r[i,])
  b2_recovered_mean[i] = mean(b2r[i,])
  tau1_recovered_mean[i] = mean(tau1r[i,])
  tau2_recovered_mean[i] = mean(tau2r[i,])
}

par(mfrow=c(4,3))

library(ggplot2)

everyone = seq(1:N)
allparams = data.frame(subject = everyone, alpha1 = alpha1_mean, alpha2 = alpha2_mean, alpha1_r = alpha1_recovered_mean, alpha2_r = alpha2_recovered_mean,
                       p = p_mean, p_r = p_recovered_mean,
                       omega = omega_mean, omega_r = omega_recovered_mean, lambda = lambda_mean, lambda_r = lambda_recovered_mean,
                       a1 = a1_mean, a1_r = a1_recovered_mean, a2 = a2_mean, a2_r = a2_recovered_mean, b1 = b1_mean, b1_r = b1_recovered_mean, b2 = b2_mean, b2_r = b2_recovered_mean,
                       tau1 = tau1_mean, tau1_r = tau1_recovered_mean, tau2 = tau2_mean, tau2_r = tau2_recovered_mean)

alpha1_corr<-round(cor(allparams$alpha1, allparams$alpha1_r),digits=3)
alpha1_corr
png(file="RLDDM_11_single_alpha1_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=alpha1, y=alpha1_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('alpha1_recovery_800') +
  geom_text(x=median(alpha1_mean),y=max(alpha1_recovered_mean),
            label = paste("r=",alpha1_corr))
dev.off()

alpha2_corr<-round(cor(allparams$alpha2, allparams$alpha2_r),digits=3)
alpha2_corr
png(file="RLDDM_11_single_alpha2_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=alpha2, y=alpha2_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('alpha2_recovery_800') +
  geom_text(x=median(alpha2_mean),y=max(alpha2_recovered_mean),
            label = paste("r=",alpha2_corr))
dev.off()

p_corr<-round(cor(allparams$p, allparams$p_r),digits=3)
p_corr
png(file="RLDDM_11_single_p_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=p, y=p_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('p_recovery_800') +
  geom_text(x=median(p_mean),y=max(p_recovered_mean),
            label = paste("r=",p_corr))
dev.off()

omega_corr<-round(cor(allparams$omega, allparams$omega_r),digits=3)
omega_corr
png(file="RLDDM_11_single_omega_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=omega, y=omega_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('omega_recovery_800') +
  geom_text(x=median(omega_mean),y=max(omega_recovered_mean),
            label = paste("r=",omega_corr))
dev.off()

lambda_corr<-round(cor(allparams$lambda, allparams$lambda_r),digits=3)
lambda_corr
png(file="RLDDM_11_single_lambda_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=lambda, y=lambda_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('lambda_recovery_800') +
  geom_text(x=median(lambda_mean),y=max(lambda_recovered_mean),
            label = paste("r=",lambda_corr))
dev.off()

a1_corr<-round(cor(allparams$a1, allparams$a1_r),digits=3)
a1_corr
png(file="RLDDM_11_single_a1_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=a1, y=a1_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('a1_recovery_800') +
  geom_text(x=median(a1_mean),y=max(a1_recovered_mean),
            label = paste("r=",a1_corr))
dev.off()

a2_corr<-round(cor(allparams$a2, allparams$a2_r),digits=3)
a2_corr
png(file="RLDDM_11_single_a2_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=a2, y=a2_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('a2_recovery_800') +
  geom_text(x=median(a2_mean),y=max(a2_recovered_mean),
            label = paste("r=",a2_corr))
dev.off()

b1_corr<-round(cor(allparams$b1, allparams$b1_r),digits=3)
b1_corr
png(file="RLDDM_11_single_b1_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=b1, y=b1_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('b1_recovery_800') +
  geom_text(x=median(b1_mean),y=max(b1_recovered_mean),
            label = paste("r=",b1_corr))
dev.off()

b2_corr<-round(cor(allparams$b2, allparams$b2_r),digits=3)
b2_corr
png(file="RLDDM_11_single_b2_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=b2, y=b2_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('b2_recovery_800') +
  geom_text(x=median(b2_mean),y=max(b2_recovered_mean),
            label = paste("r=",b2_corr))
dev.off()

tau1_corr<-round(cor(allparams$tau1, allparams$tau1_r),digits=3)
tau1_corr
png(file="RLDDM_11_single_tau1_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=tau1, y=tau1_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('tau1_recovery_800') +
  geom_text(x=median(tau1_mean),y=max(tau1_recovered_mean),
            label = paste("r=",tau1_corr))
dev.off()

tau2_corr<-round(cor(allparams$tau2, allparams$tau2_r),digits=3)
tau2_corr
png(file="RLDDM_11_single_tau2_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=tau2, y=tau2_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('tau2_recovery_800') +
  geom_text(x=median(tau2_mean),y=max(tau2_recovered_mean),
            label = paste("r=",tau2_corr))
dev.off()


save.image(file='model3.RData')


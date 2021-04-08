
#rm(list=ls())  # remove all variables 

library(rstan)
library(dplyr)

set.seed(08826)

# source HDIofMCMC.R to calculate HDI
source("HDIofMCMC.R") 

# read the data file
dat = read.table("TST_nspn.txt", header=T, sep="\t")
dat1 = filter(dat, measurment==1)
dat2 = filter(dat, measurment==2)
#generate 50 subjects from N
allSubjs = unique(dat2$subject)  # all subject IDs

N = length(allSubjs)
N # number of subjects
s <- c(186, 222,  67, 480,  51)
s
dat2 = filter(dat2, subject %in% s)
dat2

allSubjs = unique(dat2$subject)  # all subject IDs

N = length(allSubjs)
N # number of subjects
T = max(table(dat2$subject))  # max number of trials
T

###maybe we should see what T is; is it different for each subject? If that's the case
Tsubj = as.numeric(as.list(table(dat2$subject)))
numIter = 500             # number of iterations to find global minimum values
numPars = 7               # number of parameters

t = split(dat2$transition, dat2$subject, drop = FALSE)
transition = t[1]
for (i in 2:N){
  transition = rbind(transition,t[i])
}
transition = lapply(transition,as.numeric)
transition = lapply(transition, 'length<-', max(Tsubj))
transition = lapply(transition, function(x) replace(x, is.na(x), 0))

s = split((dat2$second_stage_state+1), dat2$subject, drop = FALSE)
secstagestate = s[1]
for (i in 2:N){ # skip loop if hierarchical subject
  secstagestate = rbind(secstagestate,s[i])
}
secstagestate = lapply(secstagestate,as.numeric)
secstagestate = lapply(secstagestate, 'length<-', max(Tsubj))
secstagestate = lapply(secstagestate, function(x) replace(x, is.na(x), 2))

c1 = split(dat2$choice1, dat2$subject, drop = FALSE)
choice1 = c1[1]
for (i in 2:N){ # skip loop if hierarchical subject
  choice1 = rbind(choice1,c1[i])
}
choice1 = lapply(choice1,as.numeric)
choice1 = lapply(choice1, 'length<-', max(Tsubj))
choice1 = lapply(choice1, function(x) replace(x, is.na(x), 1))

c2 = split(dat2$choice2, dat2$subject, drop = FALSE)
choice2 = c2[1]
for (i in 2:N){ # skip loop if hierarchical subject
  choice2 = rbind(choice2,c2[i])
}
choice2 = lapply(choice2,as.numeric)
choice2 = lapply(choice2, 'length<-', max(Tsubj))
choice2 = lapply(choice2, function(x) replace(x, is.na(x), 1))

rt1 = split(dat2$RT1, dat2$subject, drop = FALSE)
RT1 = rt1[1]
for (i in 2:N){ # skip loop if hierarchical subject
  RT1 = rbind(RT1,rt1[i])
}
RT1 = lapply(RT1,as.numeric)
RT1 = lapply(RT1, 'length<-', max(Tsubj))
RT1 = lapply(RT1, function(x) replace(x, is.na(x), 0))

rt2 = split(dat2$RT2, dat2$subject, drop = FALSE)
RT2 = rt2[1]
for (i in 2:N){ # skip loop if single subject
  RT2 = rbind(RT2,rt2[i])
}
RT2 = lapply(RT2,as.numeric)
RT2 = lapply(RT2, 'length<-', max(Tsubj))
RT2 = lapply(RT2, function(x) replace(x, is.na(x), 0))

r = split(dat2$reward, dat2$subject, drop = FALSE)
reward = r[1]
for (i in 2:N){ # skip loop if single subject
  reward = rbind(reward,r[i])
}
reward = lapply(reward,as.numeric)
reward = lapply(reward, 'length<-', max(Tsubj))
reward = lapply(reward, function(x) replace(x, is.na(x), 0))

rt1 = split(c(dat2$RT1,dat2$RT2), dat2$subject, drop = FALSE)
minRT = min(rt1[[1]])
for (i in 2:N){ # skip loop if single subject
  minRT = rbind(minRT,min(rt1[[i]]))
}
minRT = lapply(minRT,as.numeric)
minRT = unlist(minRT, use.names=FALSE)

RTbound = 0.1


# use all data
dataList2 <- list(
  N                = N,
  T                = T,
  Tsubj            = Tsubj,
  transition       = transition,
  secstagestate    = secstagestate,   # absolute value
  choice1          = choice1,
  choice2          = choice2,
  RT1              = RT1,
  RT2              = RT2,
  reward           = reward,
  minRT            = minRT,
  RTbound          = RTbound
)


#run!
model2 = stan("RLDDM_7_hierarchical.stan", data = dataList2, pars = c("alpha", "p", "omega", "lambda", 
                                                                "a", "b", "tau","mu_alpha","mu_pp","mu_omega","mu_lambda",
                                                                "mu_a","mu_b","mu_tau","sigma",
                                                                "rt1","rt2","ch1","ch2","arwrd","asecstagestate","trn"),
              iter = 4000, warmup=2000, chains=4, cores=4)

save.image(file='model2.RData')

# traceplot
png(file=paste("model2_group_mean_traceplot.png",sep=''),width=600, height=350)
traceplot(model2, pars=c('mu_alpha','mu_pp','mu_omega',
                         'mu_lambda','mu_a','mu_b','mu_tau'), inc_warmup=TRUE)
dev.off()
png(file=paste("model2_group_sigma_traceplot.png",sep=''),width=600, height=350)
traceplot(model2, pars=c('sigma'), inc_warmup=TRUE)
dev.off()

# traceplot & posterior histogram - put in loop

for (i in 1:N){
  png(file=paste("model2_",i,"_traceplot.png",sep=''),width=600, height=350)
  t<-traceplot(model2, pars=c(paste('alpha[',i,']',sep=''),paste('p[',i,']',sep=''),paste('omega[',i,']',sep=''),
                              paste('lambda[',i,']',sep=''),paste('a[',i,']',sep=''),paste('b[',i,']',sep=''),paste('tau[',i,']',sep='')), inc_warmup=FALSE)
  print(t)
  dev.off()
  
}

# print summary
print(model2, pars=c("alpha"), include=TRUE)
print(model2, pars=c("p"), include=TRUE)
print(model2, pars=c("lambda"), include=TRUE)
print(model2, pars=c("omega"), include=TRUE)
print(model2, pars=c("a"), include=TRUE)
print(model2, pars=c("b"), include=TRUE)
print(model2, pars=c("tau"), include=TRUE)
print(model2, pars=c("mu_alpha","mu_pp","mu_omega","mu_lambda",
                     "mu_a","mu_b","mu_tau","sigma"), include=TRUE)

# extract Stan fit object (parameters)
parameters <- rstan::extract(model2)

library(bayesplot)

png(file="model2_alpha_posterior.png",width=2000, height=1200)
plot(model2, plotfun = "hist", pars = c("alpha"), include = TRUE)
dev.off()
png(file="model2_p_posterior.png",width=2000, height=1200)
plot(model2, plotfun = "hist", pars = c("p"), include = TRUE)
dev.off()
png(file="model2_omega_posterior.png",width=2000, height=1200)
plot(model2, plotfun = "hist", pars = c("omega"), include = TRUE)
dev.off()
png(file="model2_lambda_posterior.png",width=2000, height=1200)
plot(model2, plotfun = "hist", pars = c("lambda"), include = TRUE)
dev.off()
png(file="model2_a_posterior.png",width=2000, height=1200)
plot(model2, plotfun = "hist", pars = c('a'), include = TRUE)
dev.off()
png(file="model2_b_posterior.png",width=1000, height=600)
plot(model2, plotfun = "hist", pars = c('b'), include = TRUE)
dev.off()
png(file="model2_tau_posterior.png",width=1000, height=600)
plot(model2, plotfun = "hist", pars = c('tau'), include = TRUE)
dev.off()
png(file="model2_mu_posterior.png",width=1000, height=600)
plot(model2, plotfun = "hist", pars = c('mu_alpha','mu_pp','mu_omega','mu_lambda',
                                        'mu_a','mu_b','mu_tau'), include = TRUE)
dev.off()
png(file="model2_sigma_posterior.png",width=1000, height=600)
plot(model2, plotfun = "hist", pars = c('sigma'), include = TRUE)
dev.off()

# 95% HDI of rho
for (i in 1:N){
  assign(paste('alpha_95',i,sep=''),HDIofMCMC(parameters$alpha[,i], credMass = 0.95))
  assign(paste('p_95',i,sep=''),HDIofMCMC(parameters$p[,i], credMass = 0.95))
  assign(paste('omega_95',i,sep=''),HDIofMCMC(parameters$omega[,i], credMass = 0.95))
  assign(paste('lambda_95',i,sep=''),HDIofMCMC(parameters$lambda[,i], credMass = 0.95))
  assign(paste('a_95',i,sep=''),HDIofMCMC(parameters$a[,i], credMass = 0.95))
  assign(paste('b_95',i,sep=''),HDIofMCMC(parameters$b[,i], credMass = 0.95))
  assign(paste('tau_95',i,sep=''),HDIofMCMC(parameters$tau[,i], credMass = 0.95))
}
assign(paste('sigma_95',i,sep=''),HDIofMCMC(parameters$sigma, credMass = 0.95))
assign(paste('mu_alpha_95',i,sep=''),HDIofMCMC(parameters$mu_alpha, credMass = 0.95))
assign(paste('mu_p_95',i,sep=''),HDIofMCMC(parameters$mu_pp, credMass = 0.95))
assign(paste('mu_omega_95',i,sep=''),HDIofMCMC(parameters$mu_omega, credMass = 0.95))
assign(paste('mu_lambda_95',i,sep=''),HDIofMCMC(parameters$mu_lambda, credMass = 0.95))
assign(paste('mu_a_95',i,sep=''),HDIofMCMC(parameters$mu_a, credMass = 0.95))
assign(paste('mu_b_95',i,sep=''),HDIofMCMC(parameters$mu_b, credMass = 0.95))
assign(paste('mu_tau_95',i,sep=''),HDIofMCMC(parameters$mu_tau, credMass = 0.95))











##PPC (200 trial simulation)

#first we need to define datalist for simulation

alpha = parameters$alpha
alpha_mean = rep(NA, length=N)

p = parameters$p
p_mean = rep(NA, length=N)

omega = parameters$omega
omega_mean = rep(NA, length=N)

lambda = parameters$lambda
lambda_mean = rep(NA, length=N)

a = parameters$a
a_mean = rep(NA, length=N)

b = parameters$b
b_mean = rep(NA, length=N)


tau = parameters$tau
tau_mean = rep(NA, length=N)


for (i in 1:N){
  alpha_mean[i] = mean(alpha[i,])
  p_mean[i] = mean(p[i,])
  omega_mean[i] = mean(omega[i,])
  lambda_mean[i] = mean(lambda[i,])
  a_mean[i] = mean(a[i,])
  b_mean[i] = mean(b[i,])
  tau_mean[i] = mean(tau[i,])
}


#plot posterior predictive check
library(bayesplot)

#for each subject

ch1_rep <- array(parameters$ch1,c(T,N,2000))

ch2_rep <- array(parameters$ch2,c(T,N,2000))

rt1_rep <- array(parameters$rt1,c(T,N,2000))

rt2_rep <- array(parameters$rt2,c(T,N,2000))

for (i in 1:N){
  png(file=paste("RLDDM_7_hierarchical_ch1_subject",i,"__ppc.png",sep=''),width=600, height=350)
  ch1<-ppc_bars(unlist(choice1[i]), matrix(ch1_rep[,i,],ncol=T), freq=FALSE)
  print(ch1)
  dev.off()
  png(file=paste("RLDDM_7_hierarchical_ch2_subject",i,"__ppc.png",sep=''),width=600, height=350)
  ch2<-ppc_bars(unlist(choice2[i]), matrix(ch2_rep[,i,],ncol=T), freq=FALSE)
  print(ch2)
  dev.off()
  png(file=paste("RLDDM_7_hierarchical_rt1_subject",i,"__ppc.png",sep=''),width=600, height=350)
  rt1<-ppc_stat(unlist(RT1[i]), matrix(rt1_rep[,i,],ncol=T)) 
  print(rt1)
  dev.off()
  png(file=paste("RLDDM_7_hierarchical_rt2_subject",i,"__ppc.png",sep=''),width=600, height=350)
  rt2<-ppc_stat(unlist(RT2[i]), matrix(rt2_rep[,i,],ncol=T))
  print(rt2)
  dev.off()
}











#parameter recovery

##for simulating 1000 trials

dataList1_sim <- list(
  N                = N,
  T                = 1000,
  Tsubj            = rep(1000,N),
  alpha            = alpha_mean,
  p                = p_mean,
  omega            = omega_mean,
  lambda           = lambda_mean,
  a                = a_mean,
  b                = b_mean,
  tau              = tau_mean
)

model2_gen_1000 = stan("RLDDM_7_sim.stan", data = dataList1_sim, pars = c("rt1","rt2","ch1","ch2","arwrd","asecstagestate","trn"),
                       iter = 1, chains=1, cores=1, algorithm="Fixed_param")

print(model2_gen_1000)

parameters_gen_1000 <- rstan::extract(model2_gen_1000)

save.image(file='model2.RData')

#for each subject

ch1_rep_1000 <- parameters_gen_1000$ch1

ch2_rep_1000 <- parameters_gen_1000$ch2

rt1_rep_1000 <- parameters_gen_1000$rt1

rt2_rep_1000 <- parameters_gen_1000$rt2

#preprocess minRT data from generated data(from parameters, not generated data since minRT doesn't work well with that)

minRTr = min(c(rt1_rep[,1,],rt2_rep[,1,]))
for (i in 2:N){ # skip loop if hierarchical subject
  minRTr = c(minRTr,min(c(rt1_rep[,i,],rt2_rep[,i,])))
}
minRTr = lapply(minRTr,as.numeric)
minRTr = unlist(minRTr, use.names=FALSE)
minRTr
RTbound = 0.1

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
  minRT            = minRTr,
  RTbound          = RTbound
)

# Let's fit the simulated data
model2_recovery = stan("RLDDM_7_hierarchical.stan", data = dataList1_recovery, pars = c("alpha", "p", "omega", "lambda", 
                                                                                  "a", "b", "tau"),
                       iter = 4000, warmup=2000, chains=4, cores=4)

parameters_recovery <- rstan::extract(model2_recovery)
save.image(file='model2.RData')

#plot recovery posteriors if you want..


#plot recovery plots

alphar = parameters_recovery$alpha
alpha_recovered_mean = rep(NA, length=N)

pr = parameters_recovery$p
p_recovered_mean = rep(NA, length=N)

omegar = parameters_recovery$omega
omega_recovered_mean = rep(NA, length=N)

lambdar = parameters_recovery$lambda
lambda_recovered_mean = rep(NA, length=N)

ar = parameters_recovery$a
a_recovered_mean = rep(NA, length=N)

br = parameters_recovery$b
b_recovered_mean = rep(NA, length=N)

taur = parameters_recovery$tau
tau_recovered_mean = rep(NA, length=N)


for (i in 1:N){
  alpha_recovered_mean[i] = mean(alphar[i,])
  p_recovered_mean[i] = mean(pr[i,])
  omega_recovered_mean[i] = mean(omegar[i,])
  lambda_recovered_mean[i] = mean(lambdar[i,])
  a_recovered_mean[i] = mean(ar[i,])
  b_recovered_mean[i] = mean(br[i,])
  tau_recovered_mean[i] = mean(taur[i,])
}

par(mfrow=c(4,3))

library(ggplot2)

everyone = seq(1:N)
allparams = data.frame(subject = everyone, alpha = alpha_mean, alpha_r = alpha_recovered_mean,
                       p = p_mean, p_r = p_recovered_mean,
                       omega = omega_mean, omega_r = omega_recovered_mean, lambda = lambda_mean, lambda_r = lambda_recovered_mean,
                       a = a_mean, a_r = a_recovered_mean, b = b_mean, b_r = b_recovered_mean,
                       tau = tau_mean, tau_r = tau_recovered_mean)

alpha_corr<-round(cor(allparams$alpha, allparams$alpha_r),digits=3)
alpha_corr
png(file="RLDDM_7_hierarchical_alpha_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=alpha, y=alpha_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('alpha_recovery_1000') +
  geom_text(x=median(alpha_mean),y=max(alpha_recovered_mean),
            label = paste("r=",alpha_corr))
dev.off()

p_corr<-round(cor(allparams$p, allparams$p_r),digits=3)
p_corr
png(file="RLDDM_7_hierarchical_p_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=p, y=p_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('p_recovery_1000') +
  geom_text(x=median(p_mean),y=max(p_recovered_mean),
            label = paste("r=",p_corr))
dev.off()

omega_corr<-round(cor(allparams$omega, allparams$omega_r),digits=3)
omega_corr
png(file="RLDDM_7_hierarchical_omega_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=omega, y=omega_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('omega_recovery_1000') +
  geom_text(x=median(omega_mean),y=max(omega_recovered_mean),
            label = paste("r=",omega_corr))
dev.off()

lambda_corr<-round(cor(allparams$lambda, allparams$lambda_r),digits=3)
lambda_corr
png(file="RLDDM_7_hierarchical_lambda_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=lambda, y=lambda_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('lambda_recovery_1000') +
  geom_text(x=median(lambda_mean),y=max(lambda_recovered_mean),
            label = paste("r=",lambda_corr))
dev.off()

a_corr<-round(cor(allparams$a, allparams$a_r),digits=3)
a_corr
png(file="RLDDM_7_hierarchical_a_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=a, y=a_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('a_recovery_1000') +
  geom_text(x=median(a_mean),y=max(a_recovered_mean),
            label = paste("r=",a_corr))
dev.off()

b_corr<-round(cor(allparams$b, allparams$b_r),digits=3)
b_corr
png(file="RLDDM_7_hierarchical_b_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=b, y=b_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('b_recovery_1000') +
  geom_text(x=median(b_mean),y=max(b_recovered_mean),
            label = paste("r=",b_corr))
dev.off()


tau_corr<-round(cor(allparams$tau, allparams$tau_r),digits=3)
tau_corr
png(file="RLDDM_7_hierarchical_tau_recovery_1000.png",width=600, height=350)
ggplot(allparams, aes(x=tau, y=tau_r)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle('tau_recovery_1000') +
  geom_text(x=median(tau_mean),y=max(tau_recovered_mean),
            label = paste("r=",tau_corr))
dev.off()










##for simulating 800 trials

dataList1_sim <- list(
  N                = N,
  T                = 800,
  Tsubj            = rep(800,N),
  alpha            = alpha_mean,
  p                = p_mean,
  omega            = omega_mean,
  lambda           = lambda_mean,
  a                = a_mean,
  b                = b_mean,
  tau              = tau_mean
)

model2_gen_800 = stan("RLDDM_7_sim.stan", data = dataList1_sim, pars = c("rt1","rt2","ch1","ch2","arwrd","asecstagestate","trn"),
                      iter = 1, chains=1, cores=1, algorithm="Fixed_param")

print(model2_gen_800)

parameters_gen_800 <- rstan::extract(model2_gen_800)


#for each subject

ch1_rep_800 <- parameters_gen_800$ch1

ch2_rep_800 <- parameters_gen_800$ch2

rt1_rep_800 <- parameters_gen_800$rt1

rt2_rep_800 <- parameters_gen_800$rt2

#preprocess minRT data from generated data(from parameters, not generated data since minRT doesn't work well with that)

minRTr = min(c(rt1_rep[,1,],rt2_rep[,1,]))
for (i in 2:N){ # skip loop if single subject
  minRTr = c(minRTr,min(c(rt1_rep[,i,],rt2_rep[,i,])))
}
minRTr = lapply(minRTr,as.numeric)
minRTr = unlist(minRTr, use.names=FALSE)
minRTr
RTbound = 0.1

##for 800 simulated trials (try larger trials)

dataList1_recovery_800 <- list(
  N                = N,
  T                = 800,
  Tsubj            = rep(800,N),
  transition       = parameters_gen_800$trn[1,,],
  secstagestate    = parameters_gen_800$asecstagestate[1,,]+1,
  choice1          = ch1_rep_800[1,,],
  choice2          = ch2_rep_800[1,,],
  RT1              = rt1_rep_800[1,,],
  RT2              = rt2_rep_800[1,,],
  reward           = parameters_gen_800$arwrd[1,,],
  minRT            = minRTr,
  RTbound          = RTbound
)

# Let's fit the simulated data
model2_recovery_800 = stan("RLDDM_7_hierarchical.stan", data = dataList1_recovery_800, pars = c("alpha", "p", "omega", "lambda", 
                                                                                          "a", "b", "tau"),
                           iter = 4000, warmup=2000, chains=4, cores=4)

parameters_recovery_800 <- rstan::extract(model2_recovery_800)

#plot recovery posteriors if you want



#plot recovery plots

alphar = parameters_recovery$alpha
alpha_recovered_mean = rep(NA, length=N)

pr = parameters_recovery$p
p_recovered_mean = rep(NA, length=N)

omegar = parameters_recovery$omega
omega_recovered_mean = rep(NA, length=N)

lambdar = parameters_recovery$lambda
lambda_recovered_mean = rep(NA, length=N)

ar = parameters_recovery$a
a_recovered_mean = rep(NA, length=N)

br = parameters_recovery$b
b_recovered_mean = rep(NA, length=N)

taur = parameters_recovery$tau
tau_recovered_mean = rep(NA, length=N)

for (i in 1:N){
  alpha_recovered_mean[i] = mean(alphar[i,])
  p_recovered_mean[i] = mean(pr[i,])
  omega_recovered_mean[i] = mean(omegar[i,])
  lambda_recovered_mean[i] = mean(lambdar[i,])
  a_recovered_mean[i] = mean(ar[i,])
  b_recovered_mean[i] = mean(br[i,])
  tau_recovered_mean[i] = mean(taur[i,])
}

par(mfrow=c(4,3))

library(ggplot2)

everyone = as.list(seq(1:N))
allparams_800 = data.frame(subject = everyone, alpha = alpha_mean, alpha_r = alpha_recovered_mean,
                           p = p_mean, p_r = p_recovered_mean,
                           omega = omega_mean, omega_r = omega_recovered_mean, lambda = lambda_mean, lambda_r = lambda_recovered_mean,
                           a = a_mean, a_r = a_recovered_mean, b = b_mean, b_r = b_recovered_mean,
                           tau = tau_mean, tau_r = tau_recovered_mean)

png(file="RLDDM_7_hierarchical_alpha_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=alpha, y=alpha_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_p_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=p, y=p_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_omega_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=omega, y=omega_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_lambda_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=lambda, y=lambda_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_a_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=a, y=a_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_b_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=b, y=b_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_tau_recovery_800.png",width=600, height=350)
ggplot(allparams, aes(x=tau, y=tau_r)) + geom_point()
dev.off()










##for simulating 600 trials

dataList1_sim <- list(
  N                = N,
  T                = 600,
  Tsubj            = rep(600,N),
  alpha            = alpha_mean,
  p                = p_mean,
  omega            = omega_mean,
  lambda           = lambda_mean,
  a                = a_mean,
  b                = b_mean,
  tau              = tau_mean
)

model2_gen_600 = stan("RLDDM_7_sim.stan", data = dataList1_sim, pars = c("rt1","rt2","ch1","ch2","arwrd","asecstagestate","trn"),
                      iter = 1, chains=1, cores=1, algorithm="Fixed_param")

print(model2_gen_600)

parameters_gen_600 <- rstan::extract(model2_gen_600)


#for each subject

ch1_rep_600 <- parameters_gen_600$ch1

ch2_rep_600 <- parameters_gen_600$ch2

rt1_rep_600 <- parameters_gen_600$rt1

rt2_rep_600 <- parameters_gen_600$rt2

#preprocess minRT data from generated data(from parameters, not generated data since minRT doesn't work well with that)

minRTr = min(c(rt1_rep[,1,],rt2_rep[,1,]))
for (i in 2:N){ # skip loop if single subject
  minRTr = c(minRTr,min(c(rt1_rep[,i,],rt2_rep[,i,])))
}
minRTr = lapply(minRTr,as.numeric)
minRTr = unlist(minRTr, use.names=FALSE)
minRTr
RTbound = 0.1

##for 800 simulated trials (try larger trials)

dataList1_recovery_600 <- list(
  N                = N,
  T                = 600,
  Tsubj            = rep(600,N),
  transition       = parameters_gen_600$trn[1,,],
  secstagestate    = parameters_gen_600$asecstagestate[1,,]+1,
  choice1          = ch1_rep_600[1,,],
  choice2          = ch2_rep_600[1,,],
  RT1              = rt1_rep_600[1,,],
  RT2              = rt2_rep_600[1,,],
  reward           = parameters_gen_600$arwrd[1,,],
  minRT            = minRTr,
  RTbound          = RTbound
)

# Let's fit the simulated data
model2_recovery_600 = stan("RLDDM_7_hierarchical.stan", data = dataList1_recovery_600, pars = c("alpha", "p", "omega", "lambda", 
                                                                                          "a", "b", "tau"),
                           iter = 4000, warmup=2000, chains=4, cores=4)

parameters_recovery_600 <- rstan::extract(model2_recovery_600)

#plot recovery posteriors if you want



#plot recovery plots

alphar = parameters_recovery$alpha
alpha_recovered_mean = rep(NA, length=N)

pr = parameters_recovery$p
p_recovered_mean = rep(NA, length=N)

omegar = parameters_recovery$omega
omega_recovered_mean = rep(NA, length=N)

lambdar = parameters_recovery$lambda
lambda_recovered_mean = rep(NA, length=N)

ar = parameters_recovery$a
a_recovered_mean = rep(NA, length=N)

br = parameters_recovery$b
b_recovered_mean = rep(NA, length=N)

taur = parameters_recovery$tau
tau_recovered_mean = rep(NA, length=N)

for (i in 1:N){
  alpha_recovered_mean[i] = mean(alphar[i,])
  p_recovered_mean[i] = mean(pr[i,])
  omega_recovered_mean[i] = mean(omegar[i,])
  lambda_recovered_mean[i] = mean(lambdar[i,])
  a_recovered_mean[i] = mean(ar[i,])
  b_recovered_mean[i] = mean(br[i,])
  tau_recovered_mean[i] = mean(taur[i,])
}

par(mfrow=c(4,3))

library(ggplot2)

everyone = as.list(seq(1:N))
allparams = data.frame(subject = everyone, alpha = alpha_mean, alpha_r = alpha_recovered_mean,
                       p = p_mean, p_r = p_recovered_mean,
                       omega = omega_mean, omega_r = omega_recovered_mean, lambda = lambda_mean, lambda_r = lambda_recovered_mean,
                       a = a_mean, a_r = a_recovered_mean, b = b_mean, b_r = b_recovered_mean,
                       tau = tau_mean, tau_r = tau_recovered_mean)

png(file="RLDDM_7_hierarchical_alpha_recovery_600.png",width=600, height=350)
ggplot(allparams, aes(x=alpha, y=alpha_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_p_recovery_600.png",width=600, height=350)
ggplot(allparams, aes(x=p, y=p_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_omega_recovery_600.png",width=600, height=350)
ggplot(allparams, aes(x=omega, y=omega_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_lambda_recovery_600.png",width=600, height=350)
ggplot(allparams, aes(x=lambda, y=lambda_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_a_recovery_600.png",width=600, height=350)
ggplot(allparams, aes(x=a, y=a_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_b_recovery_600.png",width=600, height=350)
ggplot(allparams, aes(x=b, y=b_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_tau_recovery_600.png",width=600, height=350)
ggplot(allparams, aes(x=tau, y=tau_r)) + geom_point()
dev.off()










##for simulating 400 trials

dataList1_sim <- list(
  N                = N,
  T                = 400,
  Tsubj            = rep(400,N),
  alpha            = alpha_mean,
  p                = p_mean,
  omega            = omega_mean,
  lambda           = lambda_mean,
  a                = a_mean,
  b                = b_mean,
  tau              = tau_mean
)

model2_gen_400 = stan("RLDDM_7_sim.stan", data = dataList1_sim, pars = c("rt1","rt2","ch1","ch2","arwrd","asecstagestate","trn"),
                      iter = 1, chains=1, cores=1, algorithm="Fixed_param")

print(model2_gen_400)

parameters_gen_400 <- rstan::extract(model2_gen_400)


#for each subject

ch1_rep_400 <- parameters_gen_400$ch1

ch2_rep_400 <- parameters_gen_400$ch2

rt1_rep_400 <- parameters_gen_400$rt1

rt2_rep_400 <- parameters_gen_400$rt2

#preprocess minRT data from generated data(from parameters, not generated data since minRT doesn't work well with that)

minRTr = min(c(rt1_rep[,1,],rt2_rep[,1,]))
for (i in 2:N){ # skip loop if hierarchical subject
  minRTr = c(minRTr,min(c(rt1_rep[,i,],rt2_rep[,i,])))
}
minRTr = lapply(minRTr,as.numeric)
minRTr = unlist(minRTr, use.names=FALSE)
minRTr
RTbound = 0.1

##for 400 simulated trials (try larger trials)

dataList1_recovery_400 <- list(
  N                = N,
  T                = 400,
  Tsubj            = rep(400,N),
  transition       = parameters_gen_400$trn[1,,],
  secstagestate    = parameters_gen_400$asecstagestate[1,,]+1,
  choice1          = ch1_rep_400[1,,],
  choice2          = ch2_rep_400[1,,],
  RT1              = rt1_rep_400[1,,],
  RT2              = rt2_rep_400[1,,],
  reward           = parameters_gen_400$arwrd[1,,],
  minRT            = minRTr,
  RTbound          = RTbound
)

# Let's fit the simulated data
model2_recovery_400 = stan("RLDDM_7_hierarchical.stan", data = dataList1_recovery_400, pars = c("alpha", "p", "omega", "lambda", 
                                                                                          "a", "b", "tau"),
                           iter = 4000, warmup=2000, chains=4, cores=4)

parameters_recovery_400 <- rstan::extract(model2_recovery_400)

#plot recovery posteriors if you want



#plot recovery plots

alphar = parameters_recovery$alpha
alpha_recovered_mean = rep(NA, length=N)

pr = parameters_recovery$p
p_recovered_mean = rep(NA, length=N)

omegar = parameters_recovery$omega
omega_recovered_mean = rep(NA, length=N)

lambdar = parameters_recovery$lambda
lambda_recovered_mean = rep(NA, length=N)

ar = parameters_recovery$a
a_recovered_mean = rep(NA, length=N)

br = parameters_recovery$b
b_recovered_mean = rep(NA, length=N)

taur = parameters_recovery$tau
tau_recovered_mean = rep(NA, length=N)

for (i in 1:N){
  alpha_recovered_mean[i] = mean(alphar[i,])
  p_recovered_mean[i] = mean(pr[i,])
  omega_recovered_mean[i] = mean(omegar[i,])
  lambda_recovered_mean[i] = mean(lambdar[i,])
  a_recovered_mean[i] = mean(ar[i,])
  b_recovered_mean[i] = mean(br[i,])
  tau_recovered_mean[i] = mean(taur[i,])
}

par(mfrow=c(4,3))

library(ggplot2)

everyone = as.list(seq(1:N))
allparams = data.frame(subject = everyone, alpha = alpha_mean, alpha_r = alpha_recovered_mean,
                       p = p_mean, p_r = p_recovered_mean,
                       omega = omega_mean, omega_r = omega_recovered_mean, lambda = lambda_mean, lambda_r = lambda_recovered_mean,
                       a = a_mean, a_r = a_recovered_mean, b = b_mean, b_r = b_recovered_mean,
                       tau = tau_mean, tau_r = tau_recovered_mean)

png(file="RLDDM_7_hierarchical_alpha_recovery_400.png",width=600, height=350)
ggplot(allparams, aes(x=alpha, y=alpha_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_p_recovery_400.png",width=600, height=350)
ggplot(allparams, aes(x=p, y=p_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_omega_recovery_400.png",width=600, height=350)
ggplot(allparams, aes(x=omega, y=omega_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_lambda_recovery_400.png",width=600, height=350)
ggplot(allparams, aes(x=lambda, y=lambda_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_a_recovery_400.png",width=600, height=350)
ggplot(allparams, aes(x=a, y=a_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_b_recovery_400.png",width=600, height=350)
ggplot(allparams, aes(x=b, y=b_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_tau_recovery_400.png",width=600, height=350)
ggplot(allparams, aes(x=tau, y=tau_r)) + geom_point()
dev.off()










##for simulating 200 trials

dataList1_sim <- list(
  N                = N,
  T                = 200,
  Tsubj            = rep(200,N),
  alpha            = alpha_mean,
  p                = p_mean,
  omega            = omega_mean,
  lambda           = lambda_mean,
  a                = a_mean,
  b                = b_mean,
  tau              = tau_mean
)

model2_gen_200 = stan("RLDDM_7_sim.stan", data = dataList1_sim, pars = c("rt1","rt2","ch1","ch2","arwrd","asecstagestate","trn"),
                      iter = 1, chains=1, cores=1, algorithm="Fixed_param")

print(model2_gen_200)

parameters_gen_200 <- rstan::extract(model2_gen_200)


#for each subject

ch1_rep_200 <- parameters_gen_200$ch1

ch2_rep_200 <- parameters_gen_200$ch2

rt1_rep_200 <- parameters_gen_200$rt1

rt2_rep_200 <- parameters_gen_200$rt2

#preprocess minRT data from generated data(from parameters, not generated data since minRT doesn't work well with that)

minRTr = min(c(rt1_rep[,1,],rt2_rep[,1,]))
for (i in 2:N){ # skip loop if single subject
  minRTr = c(minRTr,min(c(rt1_rep[,i,],rt2_rep[,i,])))
}
minRTr = lapply(minRTr,as.numeric)
minRTr = unlist(minRTr, use.names=FALSE)
minRTr
RTbound = 0.1

##for 200 simulated trials (try larger trials)

dataList1_recovery_200 <- list(
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
  minRT            = minRTr,
  RTbound          = RTbound
)

# Let's fit the simulated data
model2_recovery_200 = stan("RLDDM_7_hierarchical.stan", data = dataList1_recovery_200, pars = c("alpha", "p", "omega", "lambda", 
                                                                                          "a", "b", "tau"),
                           iter = 4000, warmup=2000, chains=4, cores=4)

parameters_recovery_200 <- rstan::extract(model2_recovery_200)

#plot recovery posteriors if you want



#plot recovery plots

alphar = parameters_recovery$alpha
alpha_recovered_mean = rep(NA, length=N)

pr = parameters_recovery$p
p_recovered_mean = rep(NA, length=N)

omegar = parameters_recovery$omega
omega_recovered_mean = rep(NA, length=N)

lambdar = parameters_recovery$lambda
lambda_recovered_mean = rep(NA, length=N)

ar = parameters_recovery$a
a_recovered_mean = rep(NA, length=N)

br = parameters_recovery$b
b_recovered_mean = rep(NA, length=N)

taur = parameters_recovery$tau
tau_recovered_mean = rep(NA, length=N)

for (i in 1:N){
  alpha_recovered_mean[i] = mean(alphar[i,])
  p_recovered_mean[i] = mean(pr[i,])
  omega_recovered_mean[i] = mean(omegar[i,])
  lambda_recovered_mean[i] = mean(lambdar[i,])
  a_recovered_mean[i] = mean(ar[i,])
  b_recovered_mean[i] = mean(br[i,])
  tau_recovered_mean[i] = mean(taur[i,])
}

par(mfrow=c(4,3))

library(ggplot2)

everyone = as.list(seq(1:N))
allparams = data.frame(subject = everyone, alpha = alpha_mean, alpha_r = alpha_recovered_mean,
                       p = p_mean, p_r = p_recovered_mean,
                       omega = omega_mean, omega_r = omega_recovered_mean, lambda = lambda_mean, lambda_r = lambda_recovered_mean,
                       a = a_mean, a_r = a_recovered_mean, b = b_mean, b_r = b_recovered_mean,
                       tau = tau_mean, tau_r = tau_recovered_mean)

png(file="RLDDM_7_hierarchical_alpha_recovery_200.png",width=600, height=350)
ggplot(allparams, aes(x=alpha, y=alpha_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_p_recovery_200.png",width=600, height=350)
ggplot(allparams, aes(x=p, y=p_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_omega_recovery_200.png",width=600, height=350)
ggplot(allparams, aes(x=omega, y=omega_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_lambda_recovery_200.png",width=600, height=350)
ggplot(allparams, aes(x=lambda, y=lambda_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_a_recovery_200.png",width=600, height=350)
ggplot(allparams, aes(x=a, y=a_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_b_recovery_200.png",width=600, height=350)
ggplot(allparams, aes(x=b, y=b_r)) + geom_point()
dev.off()

png(file="RLDDM_7_hierarchical_tau_recovery_200.png",width=600, height=350)
ggplot(allparams, aes(x=tau, y=tau_r)) + geom_point()
dev.off()









save.image(file='model2.RData')

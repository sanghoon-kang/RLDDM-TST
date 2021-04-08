######For model comparison######

###I'll only fit model parameters and not output predicted values###

loo1<-loo(model1, cores=getOption("mc.cores", 2))
loo2<-loo(model2, cores=getOption("mc.cores", 2))
loo3<-loo(model3, cores=getOption("mc.cores", 2))
loo4<-loo(model4, cores=getOption("mc.cores", 2))

loo_campare(loo1, loo2, loo3, loo4, criterion=c("loo","kfold",
                                                "waic"),detail=TRUE)

################################################################
## This code is written by Valery Fuh Ngwa                    ##
## Here we train the survival models on the time              ##
## to first relapse and compute association parameters thereof##
################################################################

library(JMbayes)
library(matrixStats)
library(gamlss)  #for non-linear regression splines
#library(JM)     # do not uncomment
library(lme4)
###Set working directory
#setwd("C:/Users/Owner/Documents/Survival_Analysis/DynPred_JointModels/DynPreds_JM/DynPred_JointModels_1")
setwd("C:/Users/Owner/Documents/Survival_Analysis/DynPred_JointModels/DynPreds_JM/PhD_Study_1/Article_2/Analyses")

## Read in the AusLong training data data sets. 
##This contains the genetic and clinical prognostic index developed before

xLong<-read.delim("xLong_Train.txt")  # AusLong training data


##Detrend the DMTD
xLong$vitD<-1
xLong$vitD[xLong$vitd=="Yes"]=-1
fit0 <- lmer(DMTD ~ Obstime + (Obstime| auslongid), data = xLong)
xLong$DMTDT<- abs(xLong$DMTD- fitted(fit0))

xSurv<-xLong[!duplicated(xLong$auslongid),]  

##Create Formulars for bayesain mixed models#
#############################################

(qt1=quantile(xLong$Obstime, c(0.33, 0.66)));range(xLong$Obstime)     
w0<-c("ns(Obstime, 3)", "(ns(Obstime,3)|auslongid)" )                                           ##Null model
w1<-c("DMTClass" , "DMTClass*log(DMTDT+0.5)"  , "ns(Obstime, 3)", "(ns(Obstime,3)|auslongid)" ) ## DMT model             
w2<-c("vitD*ns(Obstime, 3)", "(ns(Obstime,3)|auslongid)" )                                      ## VitD Model
w3<-c("vitD*DMTClass","vitD*log(DMTDT+0.5)","vitD*ns(Obstime, 3)", "(ns(Obstime,3)|auslongid)" )       ## VitD and DMT models

###RRE Formulars
r0<-c("ns(Obstime, 3)", "(ns(Obstime,3)|auslongid)" )                                                  ##Null model
r1<-c("DMTClass" , "DMTClass*log(DMTDT+0.5)"  , "ns(Obstime, 3)", "(ns(Obstime,3)|auslongid)" )        ## DMT model             
r2<-c("vitd*ns(Obstime, 3)", "(ns(Obstime,3)|auslongid)" )                                             ## VitD Model
r3<-c("vitd*DMTClass","vitd*log(DMTDT+0.5)","vitd*ns(Obstime, 3)", "(ns(Obstime,3)|auslongid)" )       ## VitD and DMT models

##Response variables
####################

resvar1<-"biomrkW"   ## Worsening specific genetic prognostic index (WS-GPI)
resvar2<-"biomrkR"   ## Relapse-specific genetic prognostic index (RS-GPI)

## Formula estimating the causal effect of WoD on RRE
formW0<-as.formula(paste(resvar1, paste(w0, collapse=" + "), sep=" ~")) 
formW1<-as.formula(paste(resvar1, paste(w1, collapse=" + "), sep=" ~"))
formW2<-as.formula(paste(resvar1, paste(w2, collapse=" + "), sep=" ~")) 
formW3<-as.formula(paste(resvar1, paste(w3, collapse=" + "), sep=" ~"))


formR0<-as.formula(paste(resvar2, paste(r0, collapse=" + "), sep=" ~")) 
formR1<-as.formula(paste(resvar2, paste(r1, collapse=" + "), sep=" ~"))
formR2<-as.formula(paste(resvar2, paste(r2, collapse=" + "), sep=" ~")) 
formR3<-as.formula(paste(resvar2, paste(r3, collapse=" + "), sep=" ~"))


###Set option for parallel processing
#####################################
options(mc.cores = parallel::detectCores()-4)

## Step 1: Fit a Linear Mixed model to describe each evolution
##############################################################
##WoD Model

lmeFitW0  <- mvglmer(list(formW0),data = xLong, families = list(gaussian),engine = "JAGS") 
lmeFitW1  <- mvglmer(list(formW1),data = xLong, families = list(gaussian),engine = "JAGS") 
lmeFitW2  <- mvglmer(list(formW2),data = xLong, families = list(gaussian),engine = "JAGS") 
lmeFitW3  <- mvglmer(list(formW3),data = xLong, families = list(gaussian),engine = "JAGS")

##Extract Performance statistics
#################################
mods<-list(lmeFitW0, lmeFitW1, lmeFitW2, lmeFitW3)
dic<-c(rep(NA, 4))
pD<-c(rep(NA, 4))
dbar<-c(rep(NA, 4))
dtheta<-c(rep(NA, 4)) 
loglik<-c(rep(NA, 4))
for(i in 1:4){
  dic[i]<-mods[[i]]$DIC
  pD[i]<-mods[[i]]$pD
  dbar[i]<-dic[i]-pD[i]
  dtheta[i]<- -pD[i]+dbar[i]
  loglik[i]<--2*(-dic[i]-2*pD[i])
}
dicW<-matrix(c(dbar, dtheta, loglik, pD, dic), nrow = 4, ncol =5 )
rownames(dicW)<-c("nullModel", "DMTModel", "vitDModel", "DMTxvitDModel")
colnames(dicW)<-c("dbar", "dtheta", "loglik", "pD", "dic")

##Fit RRE Models
################

lmeFitR0  <- mvglmer(list(formR0),data = xLong, families = list(gaussian),engine = "JAGS") 
lmeFitR1  <- mvglmer(list(formR1),data = xLong, families = list(gaussian),engine = "JAGS") 
lmeFitR2  <- mvglmer(list(formR2),data = xLong, families = list(gaussian),engine = "JAGS") 
lmeFitR3  <- mvglmer(list(formR3),data = xLong, families = list(gaussian),engine = "JAGS")

##Extract Performance statistics
#################################
mods<-list(lmeFitR0, lmeFitR1, lmeFitR2, lmeFitR3)
dic<-c(rep(NA, 4))
pD<-c(rep(NA, 4))
dbar<-c(rep(NA, 4))
dtheta<-c(rep(NA, 4)) 
loglik<-c(rep(NA, 4))
for(i in 1:4){
  dic[i]<-mods[[i]]$DIC
  pD[i]<-mods[[i]]$pD
  dbar[i]<-dic[i]-pD[i]
  dtheta[i]<- -pD[i]+dbar[i]
  loglik[i]<--2*(-dic[i]-2*pD[i])
}
dicR<-matrix(c(dbar, dtheta, loglik, pD, dic), nrow = 4, ncol =5 )
rownames(dicR)<-c("nullModel", "DMTModel", "vitDModel", "DMTxvitDModel")
colnames(dicR)<-c("dbar", "dtheta", "loglik", "pD", "dic")
dicR
########################################################
####Stage 2: Train Survival submodels on baseline data##
########################################################

##Proportional baseline Cox Model
##################################
xSurv$wstatus=tdCox.pbc
CoxREL <- coxph(Surv(timeR, RRE)  ~ Sex+I(age)+I(log(age))+RelCont+LatExp+rT2LScore+WoD, data = xSurv, model = TRUE, x=TRUE)
CoxWDS <- coxph(Surv(timeW, WoD)  ~ Sex+I(age)+I(log(age))+RelCont+LatExp+wT2LScore+RRE, data = xSurv, model = TRUE, x=TRUE)
summary(CoxREL);
summary(CoxWDS)

##PH Models with 0.3 months lag-effect
########################################

## What happens 3months before time to first events
CoxLag1 <- coxph(Surv(pmax(timeR-0.25,0), RRE, type = "right") ~ Sex+I(age)+I(log(age))+RelCont+LatExp+rT2LScore, data = xSurv, model = TRUE, x=TRUE)
CoxLag2 <- coxph(Surv(pmax(timeW-0.25,0), WoD, type = "right") ~ Sex+I(age)+I(log(age))+RelCont+LatExp+rT2LScore, data = xSurv, model = TRUE, x=TRUE)
summary(CoxLag1);
summary(CoxLag2)
####Final Model is Time-Dependent Cox Model
############################################
### Relapse
rLong      <- tmerge(xSurv, xSurv, id=auslongid, rstatus=event(timeR, RRE), wstatus=tdc(timeW))
rLong      <- rLong[order(rLong$id, rLong$tstart, rLong$tstop),]
CoxREL.tdc <- coxph(Surv(tstart, tstop, rstatus, type = "counting") ~ Sex+I(age)+I(log(age))+RelCont+LatExp+wT2LScore+wstatus +cluster(id),data=rLong, model = TRUE, x=TRUE)
CoxREL.tdc0 <- coxph(Surv(tstart,tstop, rstatus, type = "counting") ~ zbetaW+wT2LScore+wstatus+cluster(id), data=rLong, model = TRUE, x=TRUE)
###Worsening
wLong <- tmerge(xSurv, xSurv, id=auslongid, wstatus=event(timeW, WoD),rstatus=tdc(timeR))
wLong <- wLong[order(wLong$id, wLong$tstart, wLong$tstop),]
CoxWDS.tdc  <- coxph(Surv(tstart,tstop, wstatus, type = "counting") ~ Sex+I(age)+I(log(age))+RelCont+LatExp+wT2LScore+rstatus+cluster(id), data=wLong, model = TRUE, x=TRUE)
CoxWDS.tdc1 <- coxph(Surv(tstart,tstop, wstatus, type = "counting") ~ zbetaR+wT2LScore+rstatus+cluster(id), data=wLong, model = TRUE, x=TRUE)

##########################################################################
##Stage 3: Joint Models: Connect the Mixed Models with Survival models  ###
###########################################################################

#############################################
##Model 1: Constant Coefficient Joint Models#
#############################################

##Proportional Hazards
#######################

JMFitW0 <- mvJointModelBayes(lmeFitW0, CoxREL,  timeVar = "Obstime") #Cox PH
JMFitW1 <- mvJointModelBayes(lmeFitW1, CoxREL,  timeVar = "Obstime")  #Cox PH
JMFitW2 <- mvJointModelBayes(lmeFitW2, CoxREL,  timeVar = "Obstime")  #Cox PH
JMFitW3 <- mvJointModelBayes(lmeFitW3, CoxREL,  timeVar = "Obstime")  #Cox PH

JMFitR0 <- mvJointModelBayes(lmeFitR0, CoxWDS,  timeVar = "Obstime") #Cox PH
JMFitR1 <- mvJointModelBayes(lmeFitR1, CoxWDS,  timeVar = "Obstime")  #Cox PH
JMFitR2 <- mvJointModelBayes(lmeFitR2, CoxWDS,  timeVar = "Obstime")  #Cox PH
JMFitR3 <- mvJointModelBayes(lmeFitR3, CoxWDS,  timeVar = "Obstime")  #Cox PH

#Time-Dependent Cox Models 
#############################

##Predict Relapse
JMFitW0.clin <- mvJointModelBayes(lmeFitW0, CoxREL.tdc,  timeVar = "Obstime")  # Original Model
JMFitW0.cpi  <- mvJointModelBayes(lmeFitW0, CoxREL.tdc0, timeVar = "Obstime")  # Model with CPI as confounder
JMFitW1.clin <- mvJointModelBayes(lmeFitW1, CoxREL.tdc,  timeVar = "Obstime")  # Original Model
JMFitW1.cpi  <- mvJointModelBayes(lmeFitW1, CoxREL.tdc0, timeVar = "Obstime")  # Model with CPI as confounder
JMFitW2.clin <- mvJointModelBayes(lmeFitW2, CoxREL.tdc,  timeVar = "Obstime")  # Original Model
JMFitW2.cpi  <- mvJointModelBayes(lmeFitW2, CoxREL.tdc0, timeVar = "Obstime")  # Model with CPI as confounder
JMFitW3.clin <- mvJointModelBayes(lmeFitW3, CoxREL.tdc,  timeVar = "Obstime")  # Original Model
JMFitW3.cpi  <- mvJointModelBayes(lmeFitW3, CoxREL.tdc0, timeVar = "Obstime")  # Model with CPI as confounder

summary(JMFitW0.clin);summary(JMFitW0.cpi) ## Null Models
summary(JMFitW1.clin);summary(JMFitW1.cpi) ## DMT Models
summary(JMFitW2.clin);summary(JMFitW2.cpi) ## VitD Models
summary(JMFitW3.clin);summary(JMFitW3.cpi) ## DMT X VitD Models

##Predict Worsening
###################

JMFitR0.clin <- mvJointModelBayes(lmeFitR0, CoxWDS.tdc,  timeVar = "Obstime")
JMFitR0.cpi  <- mvJointModelBayes(lmeFitR0, CoxWDS.tdc1, timeVar = "Obstime")
JMFitR1.clin <- mvJointModelBayes(lmeFitR1, CoxWDS.tdc,  timeVar = "Obstime")
JMFitR1.cpi  <- mvJointModelBayes(lmeFitR1, CoxWDS.tdc1, timeVar = "Obstime")
JMFitR2.clin <- mvJointModelBayes(lmeFitR2, CoxWDS.tdc,  timeVar = "Obstime")
JMFitR2.cpi  <- mvJointModelBayes(lmeFitR2, CoxWDS.tdc1, timeVar = "Obstime")
JMFitR3.clin <- mvJointModelBayes(lmeFitR3, CoxWDS.tdc,  timeVar = "Obstime")
JMFitR3.cpi  <- mvJointModelBayes(lmeFitR3, CoxWDS.tdc1, timeVar = "Obstime")

summary(JMFitR0.clin);summary(JMFitR0.cpi) ## Null Models
summary(JMFitR1.clin);summary(JMFitR1.cpi) ## DMT Models
summary(JMFitR2.clin);summary(JMFitR2.cpi) ## VitD Models
summary(JMFitR3.clin);summary(JMFitR3.cpi) ## DMT X VitD Models

############################## END CCJM ####################################

#####################################
## Varying Coefficient Joint Models##
#####################################

### Create TVE
tve1 <-list("biomrkW_value" = ~ 0 + tve(timeR,   df = 8))
tve2 <-list("biomrkR_value" = ~ 0 + tve(timeW,   df = 8))
tve1.tdc <-list("biomrkW_value" = ~ 0 + tve(tstart,   df = 8))
tve2.tdc <-list("biomrkR_value" = ~ 0 + tve(tstart,   df = 8))

##Update Proportional Hazard Models
JMFitW0.tve.ph   <- update(JMFitW0, timeVar = "Obstime",    Interactions = tve1)
JMFitW1.tve.ph   <- update(JMFitW1, timeVar = "Obstime",    Interactions = tve1)
JMFitW2.tve.ph   <- update(JMFitW2, timeVar = "Obstime",    Interactions = tve1)
JMFitW3.tve.ph   <- update(JMFitW3, timeVar = "Obstime",    Interactions = tve1)

#JMFitW0.tve.ph   <- update(JMFitW0, timeVar = "Time2",    Interactions = tve1)
#JMFitW1.tve.ph   <- update(JMFitW1, timeVar = "Time2",    Interactions = tve1)
#JMFitW2.tve.ph   <- update(JMFitW2, timeVar = "Time2",    Interactions = tve1)
#JMFitW3.tve.ph   <- update(JMFitW3, timeVar = "Time2",    Interactions = tve1)

JMFitR0.tve.ph   <- update(JMFitR0, timeVar = "Obstime",    Interactions = tve2)
JMFitR1.tve.ph   <- update(JMFitR1, timeVar = "Obstime",    Interactions = tve2)
JMFitR2.tve.ph   <- update(JMFitR2, timeVar = "Obstime",    Interactions = tve2)
JMFitR3.tve.ph   <- update(JMFitR3, timeVar = "Obstime",    Interactions = tve2)

#JMFitR0.tve.ph   <- update(JMFitR0, timeVar = "Time1",    Interactions = tve2)
#JMFitR1.tve.ph   <- update(JMFitR1, timeVar = "Time1",    Interactions = tve2)
#JMFitR2.tve.ph   <- update(JMFitR2, timeVar = "Time1",    Interactions = tve2)
#JMFitR3.tve.ph   <- update(JMFitR3, timeVar = "Time1",    Interactions = tve2)

summary(JMFitR3.tve.cpi);summary(JMFitR.tve)

############################################
##Update Time-dependent Cox for WoD Outcome
############################################

JMFitW0.tve<- update(JMFitW0.clin, timeVar = "Time2",  Interactions = tve1.tdc) ## allow effects to vary with the 
JMFitW1.tve<- update(JMFitW1.clin, timeVar = "Time2",  Interactions = tve1.tdc) ## of the inetr-attack inetrvals
JMFitW2.tve<- update(JMFitW2.clin, timeVar = "Time2",  Interactions = tve1.tdc) ## i.e event stop - start
JMFitW3.tve<- update(JMFitW3.clin, timeVar = "Time2",  Interactions = tve1.tdc)

#JMFitW0.tve<- update(JMFitW0.clin, timeVar = "Obstime",  Interactions = tve1.tdc) ## Vary with longitudinal followup period
#JMFitW1.tve<- update(JMFitW1.clin, timeVar = "Obstime",  Interactions = tve1.tdc)
#JMFitW2.tve<- update(JMFitW2.clin, timeVar = "Obstime",  Interactions = tve1.tdc)
#JMFitW3.tve<- update(JMFitW3.clin, timeVar = "Obstime",  Interactions = tve1.tdc)

##models with CPI adjustment
JMFitW0.tve.cpi<- update(JMFitW0.cpi, timeVar = "Time2",  Interactions = tve1.tdc) ## allow effects to vary with the 
JMFitW1.tve.cpi<- update(JMFitW1.cpi, timeVar = "Time2",  Interactions = tve1.tdc) ## of the inetr-attack inetrvals
JMFitW2.tve.cpi<- update(JMFitW2.cpi, timeVar = "Time2",  Interactions = tve1.tdc) ## i.e event stop - start
JMFitW3.tve.cpi<- update(JMFitW3.cpi, timeVar = "Time2",  Interactions = tve1.tdc)

#JMFitW0.tve.cpi<- update(JMFitW0.cpi, timeVar ="Obstime",Interactions = tve1.tdc) ## Vary with longitudinal followup period
#JMFitW1.tve.cpi<- update(JMFitW1.cpi, timeVar ="Obstime",Interactions = tve1.tdc) 
#JMFitW2.tve.cpi<- update(JMFitW2.cpi, timeVar ="Obstime",Interactions = tve1.tdc) 
#JMFitW3.tve.cpi<- update(JMFitW3.cpi, timeVar ="Obstime",Interactions = tve1.tdc)

############################################
##Update Time-dependent Cox for RRE Outcome
############################################

JMFitR0.tve<- update(JMFitR0.clin, timeVar = "Time1",  Interactions = tve2.tdc) ## allow effects to vary with the 
JMFitR1.tve<- update(JMFitR1.clin, timeVar = "Time1",  Interactions = tve2.tdc) ## of the inetr-attack inetrvals
JMFitR2.tve<- update(JMFitR2.clin, timeVar = "Time1",  Interactions = tve2.tdc) ## i.e event stop - start
JMFitR3.tve<- update(JMFitR3.clin, timeVar = "Time1",  Interactions = tve2.tdc)

#JMFitR0.tve<- update(JMFitR0.clin, timeVar = "Obstime",  Interactions = tve2.tdc) ## Vary with longitudinal followup period
#JMFitR1.tve<- update(JMFitR1.clin, timeVar = "Obstime",  Interactions = tve2.tdc)
#JMFitR2.tve<- update(JMFitR2.clin, timeVar = "Obstime",  Interactions = tve2.tdc)
#JMFitR3.tve<- update(JMFitR3.clin, timeVar = "Obstime",  Interactions = tve2.tdc)

##models with CPI adjustment
JMFitR0.tve.cpi<- update(JMFitR0.cpi, timeVar = "Time1",  Interactions = tve2.tdc) ## allow effects to vary with the 
JMFitR1.tve.cpi<- update(JMFitR1.cpi, timeVar = "Time1",  Interactions = tve2.tdc) ## of the inetr-attack inetrvals
JMFitR2.tve.cpi<- update(JMFitR2.cpi, timeVar = "Time1",  Interactions = tve2.tdc) ## i.e event stop - start
JMFitR3.tve.cpi<- update(JMFitR3.cpi, timeVar = "Time1",  Interactions = tve2.tdc)

#JMFitR0.tve.cpi<- update(JMFitR0.cpi, timeVar ="Obstime",Interactions = tve2.tdc) ## Vary with longitudinal followup period
#JMFitR1.tve.cpi<- update(JMFitR1.cpi, timeVar ="Obstime",Interactions = tve2.tdc) 
#JMFitR2.tve.cpi<- update(JMFitR2.cpi, timeVar ="Obstime",Interactions = tve2.tdc) 
#JMFitR3.tve.cpi<- update(JMFitR3.cpi, timeVar ="Obstime",Interactions = tve2.tdc)

############################## END VCJM ####################################

save(lmeFitR0, lmeFitR1, lmeFitR2, lmeFitR3,
     lmeFitW0, lmeFitW1, lmeFitW2, lmeFitW3,
     CoxREL, CoxREL.tdc, CoxREL.tdc0, CoxWDS, 
     CoxWDS.tdc, CoxWDS.tdc1,dicW, dicR, 
     JMFitR0,JMFitR1,JMFitR2,JMFitR3,
     JMFitW0,JMFitW1,JMFitW2,JMFitW3,
     JMFitR0.clin, JMFitR1.clin,JMFitR2.clin,JMFitR3.clin, 
     JMFitW0.clin, JMFitW1.clin, JMFitW2.clin,JMFitW3.clin, 
     JMFitR0.cpi, JMFitR1.cpi, JMFitR2.cpi,JMFitR3.cpi, 
     JMFitW0.cpi, JMFitW1.cpi, JMFitW2.cpi,JMFitW3.cpi, 
     JMFitW0.tve, JMFitW1.tve,JMFitW2.tve, JMFitW3.tve,
     JMFitR0.tve,JMFitR1.tve, JMFitR2.tve,JMFitR3.tve,
     JMFitW0.tve.ph, JMFitW1.tve.ph,JMFitW2.tve.ph, JMFitW3.tve.ph,
     JMFitR0.tve.ph,JMFitR1.tve.ph, JMFitR2.tve.ph,JMFitR3.tve.ph,
     JMFitR1.tve, JMFitR2.tve, JMFitR3.tve, JMFitW0.tve,
     JMFitW1.tve,JMFitW2.tve,JMFitW3.tve, JMFitR0.tve.cpi,
     JMFitR1.tve.cpi, JMFitR2.tve.cpi, JMFitR3.tve.cpi, JMFitW0.tve.cpi,
     JMFitW1.tve.cpi,JMFitW2.tve.cpi,JMFitW3.tve.cpi,file = "JointModels.RData")

####################################
## Plot Time-Varying Effects     ###
####################################


##Call the Plotting Function

source("plot_tveffectsR.R")
source("plot_tveffectsW.R")


###UnAdjusted Model
##################

fitW<-JMFitW0 ##Cox PH with CCJM
fitR<-JMFitR0 ##Cox PH with CCJM

##replace the statistic of the object
fitW$statistics<-JMFitW0.clin$statistics ##td-cox with CCJM
fitR$statistics<-JMFitR0.clin$statistics ##td-cox with CCJM

##Create a new object to enable plotting of VCJM with td-cox
fitW.tve<-JMFitW0.tve.ph ##td-cox with VCJM
fitR.tve<-JMFitR0.tve.ph ##td-cox with VCJM

fitW.tve$statistics<-JMFitW0.tve$statistics ##td-cox with VCJM
fitR.tve$statistics<-JMFitR0.tve$statistics

oldpar <- par(no.readonly=TRUE) # save graphical parameters
vdist  <- hdist <- 1.2
layout(matrix(1:2, 1, 2, byrow=TRUE),widths=c(10,10),heights=c(10,10))
vdist<-4.5; hdist<-0.5
par(mar= c(vdist, 5, 4, hdist))
plot_tveffectsW(fitW.tve ,JMFitW.tve, fitW, include = "both", scale = "logHR", ccol=8, tvcol=1, conf.int.const=FALSE)
plot_tveffectsR(fitR.tve, fitR, include = "both", scale = "logHR", ccol=8, tvcol=1, conf.int.const=FALSE)
#plot_tveslope(JMFit1.tve.cvs , JMFit1, JMFit1.tve,  include = "both", scale = "slope", ccol=8, tvcol=1, conf.int.const=FALSE)
layout(matrix(1, 1, 1))
par(oldpar) #
box(which = "outer")
####################################################

### DMT Adjusted Model
######################

fitW<-JMFitW1 ##Cox PH with CCJM
fitR<-JMFitR1 ##Cox PH with CCJM

##replace the statistic of the object
fitW$statistics<-JMFitW1.clin$statistics ##td-cox with CCJM
fitR$statistics<-JMFitR1.clin$statistics ##td-cox with CCJM

##Create a new object to enable plotting of VCJM with td-cox
fitW.tve<-JMFitW1.tve.ph ##td-cox with VCJM
fitR.tve<-JMFitR1.tve.ph ##td-cox with VCJM

fitW.tve$statistics<-JMFitW1.tve$statistics ##td-cox with VCJM
fitR.tve$statistics<-JMFitR1.tve$statistics

oldpar <- par(no.readonly=TRUE) # save graphical parameters
vdist  <- hdist <- 1.2
layout(matrix(1:2, 1, 2, byrow=TRUE),widths=c(10,10),heights=c(10,10))
vdist<-4.5; hdist<-0.5
par(mar= c(vdist, 5, 4, hdist))
plot_tveffectsW(fitW.tve ,JMFitW.tve, fitW, include = "both", scale = "logHR", ccol=8, tvcol=1, conf.int.const=FALSE)
plot_tveffectsR(fitR.tve, fitR, include = "both", scale = "logHR", ccol=8, tvcol=1, conf.int.const=FALSE)
#plot_tveslope(JMFit1.tve.cvs , JMFit1, JMFit1.tve,  include = "both", scale = "slope", ccol=8, tvcol=1, conf.int.const=FALSE)
layout(matrix(1, 1, 1))
par(oldpar) #
box(which = "outer")
##########################################


##Vitamin D Adjusted Model
###########################

fitW<-JMFitW2 ##Cox PH with CCJM
fitR<-JMFitR2 ##Cox PH with CCJM

##replace the statistic of the object
fitW$statistics<-JMFitW2.clin$statistics ##td-cox with CCJM
fitR$statistics<-JMFitR2.clin$statistics ##td-cox with CCJM

##Create a new object to enable plotting of VCJM with td-cox
fitW.tve<-JMFitW2.tve.ph ##td-cox with VCJM
fitR.tve<-JMFitR2.tve.ph ##td-cox with VCJM

fitW.tve$statistics<-JMFitW2.tve$statistics ##td-cox with VCJM
fitR.tve$statistics<-JMFitR2.tve$statistics

oldpar <- par(no.readonly=TRUE) # save graphical parameters
vdist  <- hdist <- 1.2
layout(matrix(1:2, 1, 2, byrow=TRUE),widths=c(10,10),heights=c(10,10))
vdist<-4.5; hdist<-0.5
par(mar= c(vdist, 5, 4, hdist))
plot_tveffectsW(fitW.tve ,JMFitW.tve, fitW, include = "both", scale = "logHR", ccol=8, tvcol=1, conf.int.const=FALSE)
plot_tveffectsR(fitR.tve, fitR, include = "both", scale = "logHR", ccol=8, tvcol=1, conf.int.const=FALSE)
#plot_tveslope(JMFit1.tve.cvs , JMFit1, JMFit1.tve,  include = "both", scale = "slope", ccol=8, tvcol=1, conf.int.const=FALSE)
layout(matrix(1, 1, 1))
par(oldpar) #
box(which = "outer")
####################################



## Vitamin D and DMT  Adjusted Model
####################################

fitW<-JMFitW3 ##Cox PH with CCJM
fitR<-JMFitR3 ##Cox PH with CCJM

##replace the statistic of the object
fitW$statistics<-JMFitW3.clin$statistics ##td-cox with CCJM
fitR$statistics<-JMFitR3.clin$statistics ##td-cox with CCJM

##Create a new object to enable plotting of VCJM with td-cox
fitW.tve<-JMFitW3.tve.ph ##td-cox with VCJM
fitR.tve<-JMFitR3.tve.ph ##td-cox with VCJM

fitW.tve$statistics<-JMFitW3.tve$statistics ##td-cox with VCJM
fitR.tve$statistics<-JMFitR3.tve$statistics

oldpar <- par(no.readonly=TRUE) # save graphical parameters
vdist  <- hdist <- 1.2
layout(matrix(1:2, 1, 2, byrow=TRUE),widths=c(10,10),heights=c(10,10))
vdist<-4.5; hdist<-0.5
par(mar= c(vdist, 5, 4, hdist))
plot_tveffectsW(fitW.tve ,JMFitW.tve, fitW, include = "both", scale = "logHR", ccol=8, tvcol=1, conf.int.const=FALSE)
plot_tveffectsR(fitR.tve, fitR, include = "both", scale = "logHR", ccol=8, tvcol=1, conf.int.const=FALSE)
#plot_tveslope(JMFit1.tve.cvs , JMFit1, JMFit1.tve,  include = "both", scale = "slope", ccol=8, tvcol=1, conf.int.const=FALSE)
layout(matrix(1, 1, 1))
par(oldpar) #
box(which = "outer")

####################################################################
###3 months Lag Effect

oldpar <- par(no.readonly=TRUE) # save graphical parameters
vdist  <- hdist <- 1.2
layout(matrix(1:2, 1, 2, byrow=TRUE),widths=c(10,10),heights=c(10,10))
vdist<-4.5; hdist<-0.5
par(mar= c(vdist, 5, 4, hdist))
plot_tveffects(JMFitLag1.tve1,    JMFitLag1.dc,                 lag=TRUE, include = "both", scale = "logHR", ccol=8, tvcol=1, conf.int.const=FALSE)
plot_tveslope(JMFitLag1.tve.cvs1, JMFitLag1.dc, JMFitLag1.tve1, lag=TRUE, include = "both", scale = "slope", ccol=8, tvcol=1, conf.int.const=FALSE)
layout(matrix(1, 1, 1))
par(oldpar) #
box(which = "outer")


################################################
###Plot Time-Varying Marginal Treatment Effects#
################################################

max.t <- round(max(lmeFitR1$data$Obstime),0)
t <- seq(0, max.t, by = 0.01)
y <- exp(JMFitR1.clin$statistics$postMeans$betas1[2] +
           JMFitR2.clin$statistics$postMeans$betas1[6]*ns(t, 3)[,1])

plot(t, y, pch = 18, lwd=1.5,  lty = 2)
lines(t,y, pch = 19, col = "blue",lty=3, lwd=3)


TE_SS <- exp(results.marg$ss_coefs[10] + results.marg$ss_coefs[11] *
               results.marg$ss_coefs[14] * t)

## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,root.dir="C:/Users/Robert J. Tempelman/OneDrive - Michigan State University/Tempelman/Meta_analysis/Simulation")
library(tidyverse,ggplot2)
rm(list=ls())
# set default directory


## ----Simulate------------------------------------------------------------------------------------------------------------------------------------------------------------
##### THE FOLLOWING JUST SIMULATES SOME DATA BASED ON THE FOLLOWING SPECIFICATIONS   #####
##########################################################################################################
######################## BEGINNING OF SIMULATION #########################################################
##########################################################################################################
set.seed(50)  # set seed for reproducibility.
# we should all simulate exactly the same data if we use the same seed!!!

# Let's data from a number (nStudy) of studies, each being a completely randomized design involving two treatments.

setwd("C:/Users/Robert J. Tempelman/OneDrive - Michigan State University/Tempelman/Meta_analysis/Simulation")

nStudy = 25  # Number of studies
ntrt = 2     # Number of treatments 

#We'll set the variance components 

# In a 
sigmae = 5  # residual (within cow) variance
sigmacow = 5  # between-cow variance
sigmastudy = 5  # between-study variance
  # i.e., some studies are characterized by above/below average responses!
sigmastudytrt = 3  # study by trt interaction variance
  # i.e., non-zero sigmastudytrt: treatment differences may randomly differ because of differences in preparations.

(nreps =sample(seq(40,60,2),nStudy,replace=TRUE))  # simulated different number of reps per Study varying anywhere from 40 60 in increments of 2.
# ok typically studies are not really this big but wanting to demonstrate how to estimate heterogeneous residual variances per treatment in a joint analysis
# (a Bayesian analysis would be far more suitable for this.)
overallmean = 40
trteffects = seq(-1,+1,length.out=ntrt)  # treatment effects relative to overall mean
(trtlabels = toupper(letters[1:ntrt]))

Studyeffects = rnorm(sd=sqrt(sigmastudy),nStudy)  # simulated study effects

studytrtint = vector(mode="list",length= nStudy)  # simulated study x treatment effects
for (Study in 1:nStudy){
  studytrtint[[Study]]=rnorm(sd = sqrt(sigmastudytrt),ntrt)
}

# save the true study by treatment effects for later.
study_trt = unlist(studytrtint)
Study_Trtlabels = expand.grid(trtlabels,1:nStudy)
Study_Trtlabels = paste0(Study_Trtlabels$Var1,':',Study_Trtlabels$Var2)
study_trt = data.frame(Study_Trtlabels,study_trt)

Study_data <- vector(mode="list",length=nStudy)  # we'll first save the data in this list, separate study for each element

for (Study in 1:nStudy) {
  Trt = rep(trtlabels,each=nreps[Study]) # get trt label assignments for each record
  Cow = seq(1:(ntrt*nreps[Study]))  # and cow labels for each record
  Trt_effects = rep(trteffects,each=nreps[Study])  #along with their effects
  Trt_Study_effects = rep(studytrtint[[Study]],each=nreps[Study]) # and the trt*study effects
  residuals = rnorm(sd = sqrt(sigmae+sigmacow),ntrt*nreps[Study]) # and finally the residuals
  FCM = round(overallmean + Trt_effects + Studyeffects[Study] + Trt_Study_effects + residuals,2)  # generate the record and round to 2nd decimal
  Study_data[[Study]] = data.frame(Study,Trt,Cow,FCM)  #store it.
}  

overallCRDdata = bind_rows(Study_data) %>% # the entire data set
    mutate(across(!starts_with('FCM'), as.factor))
str(overallCRDdata)
write_csv(overallCRDdata,"overallCRDdata.csv")

##########################################################################################################
######################## END OF SIMULATION ###############################################################
##########################################################################################################


## ----rawdatanalysis------------------------------------------------------------------------------------------------------------------------------------------------------
# Analysis of raw data

#IPD (Individual Performance Data)
# using raw data for data analyses
# assuming equal residual variances for each study

library(emmeans)  # library useful for computing treatment means and mean differences
emm_options(pbkrtest.limit = nrow(overallCRDdata))  # allows KR degrees of freedom for specified size of dataset (default = 3000)
  # naive analysis that ignores study effects

# Common effects analysis (Ignoring Study Heterogeneity!)
fixed_analysis = lm(FCM~Trt,data=overallCRDdata)  # Fit Trt as Fixed 
summary(fixed_analysis)
fixed_Trt= emmeans(fixed_analysis,"Trt")
(fixed_effects_means = summary(fixed_Trt))
(fixed_effects_diff = pairs(fixed_Trt))

# let's model study heterogeneity
library(lme4)  # mixed model procedure in R.
# NOTE THIS MODEL PERFECTLY MATCHES THE MODEL GENERATION PROCESS
# accounts for study heterogeneity
overall_analysis = lmer(FCM~Trt+(1|Study/Trt),data=overallCRDdata)  # Fit Trt as Fixed, Study and Study*Trt as random
summary(overall_analysis)
lsmeansraw_Trt= emmeans(overall_analysis,"Trt")
(overall_meansraw = summary(lsmeansraw_Trt))
(effectsize_raw = pairs(lsmeansraw_Trt))


## ----Studyeffects--------------------------------------------------------------------------------------------------------------------------------------------------------
# This is really not important but it might be need to see how the estimated (Predicted) study effects line up with the truth
# Estimated ("Predicted") Study Effects
BLUP_Studyeffects = coef(overall_analysis)$Study[,1]-mean(coef(overall_analysis)$Study[,1])  # need to subtract the mean to express relative to zero
trt = data.frame(Studyeffects,BLUP_Studyeffects)
# Plot Estimated ("Predicted") Study Effects vs True Study Effects
ggplot(data=trt,aes(x=Studyeffects ,y=BLUP_Studyeffects)) + geom_point() + geom_abline(intercept=0,slope=1) +
  ggtitle("Predicted vs True Study Effects") +
  xlab("True effects") + 
  ylab("BLUP")



## ----studybytreatment----------------------------------------------------------------------------------------------------------------------------------------------------
# Let's do the same thing for the study by treatment effects
#Estimated study by treatment effects 
BLUP_study_trt = coef(overall_analysis)$`Trt:Study`[,1]-mean(coef(overall_analysis)$`Trt:Study`[,1]) # need to subtract the mean to express relative to zero
BLUP_study_trt = data.frame(rownames(coef(overall_analysis)$`Trt:Study`),BLUP_study_trt)
colnames(BLUP_study_trt) = c("Study_Trtlabels","BLUP_study_trt")
study_trt = study_trt %>%
   left_join(BLUP_study_trt,by=c("Study_Trtlabels"))
# Plot Estimated ("Predicted") Study by Treatment Effects vs True Study by Treatment Effects
ggplot(study_trt,aes(x=study_trt ,y=BLUP_study_trt)) + geom_point() + geom_abline(intercept=0,slope=1) +
  ggtitle("Predicted vs True Study by Treatment Effects") +
  xlab("True effects") + 
  ylab("BLUP") +
  geom_text(x=-3, y=2, label="Notice shrinkage to zero by BLUP")



## ----heteroresidual------------------------------------------------------------------------------------------------------------------------------------------------------
library(nlme)
# modeling heterogeneous residual variances across Studies and including random effects of Study and Treatment by Study
##  THIS MIGHT BE CONSIDERED THE GOLD STANDARD ANALYSIS IF INDIVIDUAL COW DATA WAS READILY AVAILABLE AND RESIDUAL VARIABILITY WAS TRULY HETEROGENEOUS ACROSS STUDIES
##  HOWEVER THIS MIGHT NOT TYPICALLY CONVERGE IF INDIVIDUAL STUDIES ARE "SMALL"
hetero.lme <- lme(FCM ~ Trt , data=overallCRDdata, random = ~ 1|Study/Trt, weights=varIdent(form = ~ 1 | Study),method="REML")  # this will not always converge
summary(hetero.lme)  
lsmeansraw_Trt_hetero= emmeans(hetero.lme,"Trt")
(overall_meansraw_hetero = summary(lsmeansraw_Trt_hetero))


# doesn't seem to be a big difference between the two (hetero versus homogeneous residual variances)
(effectsize_raw_hetero  = pairs(lsmeansraw_Trt_hetero))
effectsize_raw




## ----metafunctions-------------------------------------------------------------------------------------------------------------------------------------------------------
 #need to create some functions first.
# need to create a separate function for model
separate_model = function(df) {
  lm(FCM~Trt,data=df)
}

# function to save the mean trt difference for each Study
contrast = function(model) {
  lsmeans_Trt= emmeans(model,"Trt")
  effectsize = pairs(lsmeans_Trt)
  effectsize = as.data.frame(summary(effectsize))
  return(effectsize)
}

# function to save the trt means for each Study
lsmeans = function(model) {
  lsmeans_Trt= emmeans(model,"Trt")
  means = summary(lsmeans_Trt)
  return(means)
}



## ----metastats-----------------------------------------------------------------------------------------------------------------------------------------------------------
Study_specific_contrasts = overallCRDdata %>%
  group_by(Study) %>%
  nest() %>%
  mutate(model = map(data,separate_model)) %>%
  mutate(effectsize = map(model,contrast)) %>%
  mutate(nrec = as.numeric(map(data,nrow))) %>% 
  unnest(effectsize)    %>%
  dplyr::select(-c(data,model)) %>%
  mutate(weight = (1/SE)^2) %>%  # What i use in the mixed model
  mutate(sampvar = SE^2)         # What the r package metafor needs

head(Study_specific_contrasts)


## ----boxplot-------------------------------------------------------------------------------------------------------------------------------------------------------------
# this might not be too helpful to look at if there is a lot of differences in standard errors between studies.
boxplot(Study_specific_contrasts$estimate,main="Box plot of study-specific treatment differences")


## ----CE1-----------------------------------------------------------------------------------------------------------------------------------------------------------------
weighted.mean(Study_specific_contrasts$estimate,Study_specific_contrasts$weight)
(stderr = (sum(Study_specific_contrasts$weight))^(-1/2))


## ----CE2-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# glmmTMB is used because it allows one to fix the residual variance to 1 (just a "trick" to allow proper weightings by 1/SE^2 which already reflect residual var!)
library(glmmTMB)
# overall mean model
FM_1meta <- glmmTMB(estimate ~ 1 ,
                weights = weight,
                family = gaussian,
                REML=TRUE,
                data = Study_specific_contrasts,
                start = list(betad = log(1)),  # fix residual variance = 1
                map = list(betad = factor(NA))
)
summary(FM_1meta)


## ----CE3-----------------------------------------------------------------------------------------------------------------------------------------------------------------
library(metafor)
(res_EE = rma.uni(yi=estimate,vi=sampvar,data=Study_specific_contrasts,method="EE"))  # common effects
forest(res_EE)


## ----RM_1meta------------------------------------------------------------------------------------------------------------------------------------------------------------

library(glmmTMB)
# random effects model
RM_1meta <- glmmTMB(estimate ~ 1 + (1|Study),  # Study really involves both Study and Study*Treatment 
              weights = weight,
              family = gaussian,
              REML=TRUE,
              data = Study_specific_contrasts,
              start = list(betad = log(1)),  # fix residual variance = 1
              map = list(betad = factor(NA))
)
summary(RM_1meta)



## ----varcompest----------------------------------------------------------------------------------------------------------------------------------------------------------
vc = VarCorr(RM_1meta)  # variance components
print(vc,comp=("Variance"))



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
(res_REML = rma.uni(yi=estimate,vi=sampvar,data=Study_specific_contrasts,method="REML"))  # default
forest(res_REML)


## ----LRT1----------------------------------------------------------------------------------------------------------------------------------------------------------------

anova(res_EE,res_REML)  # likelihood ratio test to test for differences in heterogeneity.

anova(FM_1meta,RM_1meta)  # this works too.



## ----CI_PI---------------------------------------------------------------------------------------------------------------------------------------------------------------
# confirming the metafor CI on the effect size using the mixed model analysis.
Studyvar = c(((VarCorr(RM_1meta))["cond"]$cond$Study))
overallest = coef(summary(RM_1meta))$cond
cat("Confidence Interval",sep="\n")
(LCL = overallest[1,"Estimate"]- 1.96*overallest[1,"Std. Error"])
(UCL = overallest[1,"Estimate"]+ 1.96*overallest[1,"Std. Error"])

# Prediction intervals
# useful for determining uncertainty on a future experiment.
cat("Prediction Interval",sep="\n")
(LPL = overallest[1,"Estimate"]- 1.96*sqrt(overallest[1,"Std. Error"]^2 + Studyvar))
(UPL = overallest[1,"Estimate"]+ 1.96*sqrt(overallest[1,"Std. Error"]^2 + Studyvar))



## ----studymeans----------------------------------------------------------------------------------------------------------------------------------------------------------
# Suppose you have study specific means instead
# often referred to as "arm-based" (rather than contrast-based) analysis.

# We'll also save the Study_specific means in a different file.
# Their differences should correspond to the study-specific contrasts in the above file.
Study_specific_means = overallCRDdata %>%
  group_by(Study) %>%
  nest() %>%
  mutate(model = map(data,separate_model)) %>%
  mutate(lsmeans = map(model,lsmeans)) %>%
  mutate(nrec = as.numeric(map(data,nrow))) %>% 
  unnest(lsmeans)    %>%
  dplyr::select(-c(data,model,lower.CL,upper.CL)) %>%
  mutate(weights = (1/SE^2),sampvar = SE^2)

head(Study_specific_means)

boxplot(emmean~Trt,data=Study_specific_means,main="Study-specific treatment means")


## ----Model_R0------------------------------------------------------------------------------------------------------------------------------------------------------------
# this is what many animal scientists seem to be doing.
RM_0arm <- glmmTMB(emmean ~ Trt + (1|Study),
                     weights = (1/SE^2),
                     family = gaussian,
                     REML=TRUE,
                     data = Study_specific_means,
                     start = list(betad = log(1)),  # fix residual variance = 1
                     map = list(betad = factor(NA))
)
summary(RM_0arm)
lsmeans_R0_Trt= emmeans(RM_0arm,"Trt")
summary(lsmeans_R0_Trt)
(effectsize_R0_Trt = pairs(lsmeans_R0_Trt))
# this standard error is typically underreported...compare it to what we observed from analyzing the raw data!!
cat("Compare this to inferences from analyzing the raw data",sep="\n")
effectsize_raw


## ----Model_R1------------------------------------------------------------------------------------------------------------------------------------------------------------
# you need to model random inconsistency as well!   
# (this is model R1 in Madden et al. (2016))
RM_1arm <- glmmTMB(emmean ~ Trt + (1|Study/Trt),
                weights = (1/SE^2),
                family = gaussian,
                REML=TRUE,
                data = Study_specific_means,
                start = list(betad = log(1)),  # fix residual variance = 1
                map = list(betad = factor(NA))
)
# (this is also model R1 in Madden et al. (2016))
summary(RM_1arm) 
lsmeans_R1_Trt= emmeans(RM_1arm,"Trt")
summary(lsmeans_R1_Trt)
(effectsize_R1_Trt = pairs(lsmeans_R1_Trt))
cat("Compare this to inferences from analyzing the raw data",sep="\n")
effectsize_raw


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# The following is an alternative based on separate variances for each treatment
# Model R3 in Madden et al. 2016...specifies study variances to be different for each treatment (UN= unstructured)
RM_3arm <- glmmTMB(emmean ~ Trt + (Trt-1|Study),
                     weights = (1/SE^2),
                     family = gaussian,
                     REML=TRUE,
                     data = Study_specific_means,
                     start = list(betad = log(1)),  # fix residual variance = 1
                     map = list(betad = factor(NA))
)
summary(RM_3arm)
lsmeans_R3_Trt= emmeans(RM_3arm,"Trt")
summary(lsmeans_R3_Trt)
(effectsize_R3_Trt = pairs(lsmeans_R3_Trt))


## ----metafor_mv----------------------------------------------------------------------------------------------------------------------------------------------------------
#Need to use the multivariate option of metafor!!
library(metafor)
# incomplete accounting for heterogeneity -> does not model random inconsistency!
(res.mv0a <- rma.mv(emmean, sampvar, mods = ~ Trt , random = ~ 1|Study,  data=Study_specific_means)) 
# for inferences on contrasts
(res.mv0b <- rma.mv(emmean, sampvar, mods = ~ Trt -1 , random = ~ 1|Study,  data=Study_specific_means)) 
# for inferences on means

# this is equivalent to modeling Study and Study*Trt as independent random effects (as in RM_means1) 
#Model R1 in Madden et al. (2016) which is equivalent to compound symmetry
(res.mv1a <- rma.mv(emmean, sampvar, mods = ~ Trt , random = ~ Trt | Study, struct="CS", data=Study_specific_means))  
# for inferences on contrasts
(res.mv1b <- rma.mv(emmean, sampvar, mods = ~ Trt -1 , random = ~ Trt | Study, struct="CS", data=Study_specific_means)) 
# for inferences on means

#Model R2 in Madden et al. (2016)  # Heterogeneous compound symmetry
(res.mv2a <- rma.mv(emmean, sampvar, mods = ~ Trt  , random = ~ Trt | Study, struct="HCS", data=Study_specific_means)) 
# for inferences on contrasts
(res.mv2b <- rma.mv(emmean, sampvar, mods = ~ Trt -1 , random = ~ Trt | Study, struct="HCS", data=Study_specific_means)) 
# for inferences on means

# Model R3 in Madden et al. (2016)  # Unstructured
(res.mv3a <- rma.mv(emmean, sampvar, mods = ~ Trt , random = ~ Trt | Study, struct="UN", data=Study_specific_means)) 
# for inferences on contrasts
(res.mv4b <- rma.mv(emmean, sampvar, mods = ~ Trt -1 , random = ~ Trt | Study, struct="UN", data=Study_specific_means)) 
# for inferences on means



## ----LRT2----------------------------------------------------------------------------------------------------------------------------------------------------------------
anova(RM_0arm,RM_1arm,RM_3arm) # using mixed model
anova(res.mv1a,res.mv3a)


## ----SEMvSED-------------------------------------------------------------------------------------------------------------------------------------------------------------
Study_specific_SEM = Study_specific_means %>% 
   group_by(Study) %>%
   summarize(SEM=mean(SE))

Study_specific_SED = Study_specific_contrasts %>% 
   group_by(Study) %>%
   summarize(SED=mean(SE))

SEM_SED_compare = Study_specific_SEM %>%
    left_join(Study_specific_SED,by="Study") %>%
    mutate(SEM_SED_ratio = SEM/SED)

head(SEM_SED_compare)


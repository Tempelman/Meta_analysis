#' ---
#' title: "Meta analysis for Latin Square Designs"
#' author: "Robert J. Tempelman"
#' date: "2023-06-23"
#' output: html_document
#' ---
#' 
#' # This document is used to highlight considerations for meta-analysis of Latin square/crossover designs
#' 
#' 
## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE,root.dir="C:/Users/Robert J. Tempelman/OneDrive - Michigan State University/Tempelman/Meta_analysis/Simulation")
library(tidyverse,ggplot2)
# set default directory

#' 
#' ## Simulation
#' 
#' Simulate a n_trt treatment Latin square design from a number of different studies
#' 
## ----Simulate------------------------------------------------------------------------------------------------------------------------------------------------------------
##### THE FOLLOWING JUST SIMULATES SOME DATA BASED ON THE FOLLOWING SPECIFICATIONS   #####
##########################################################################################################
######################## BEGINNING OF SIMULATION #########################################################
##########################################################################################################
set.seed(100)  # set seed for reproducibility.
nStudy = 30  # Number of studies
ntrt = 3    # Number of treatments 

sigmae = 2  # residual (within cow) variance
sigmacow = 8  # between-cow variance
sigmastudy = 5  # between-study variance
sigmastudytrt = 3  # study by trt interaction variance

sigmaperiod = 2 # variance between periods although typically considered as fixed.

overallmean = 30
trteffects = c(-1,0,+1)  # effects relative to overall mean
trtlabels = toupper(letters[1:ntrt])
if (length(trteffects) != ntrt) stop ("need to change the number of treatment effects to match ntrt")
periodeffects = rnorm(sd=sqrt(sigmaperiod),nStudy*ntrt)
StudyPeriod = split(periodeffects, ceiling(seq_along(periodeffects)/ntrt))



Studyeffects = rnorm(sd=sqrt(sigmastudy),nStudy)  # simulated study effects
studytrtint = vector(mode="list",length= nStudy)  # simulated study x treatment effects
for (Study in 1:nStudy){
  studytrtint[[Study]]=rnorm(sd = sqrt(sigmastudytrt),ntrt)
}

ncows =sample(seq(ntrt*4,ntrt*10,ntrt),nStudy,replace=TRUE)  # simulated different number of replicated squares per Study
# for a balanced multiple Latin square design, ncows has to be a multiple of number of treatments

setuplatin =vector(mode="list",length=ntrt)
base = 1:ntrt
for(i in 1:ntrt){
  setuplatin[[i]]=(base+ntrt+1+i)%%ntrt+1
}  # setup for a single Latin square

Study_data <- vector(mode="list",length=nStudy)  # we'll first save the data in this list, separate study for each element
Study = 1

for (Study in 1:nStudy) {
  Trt_number = rep(unlist(setuplatin),ncows[Study]/ntrt)
  Effect_Trt = trteffects[Trt_number]
  Period = rep(1:3,ncows[Study])
  Effect_Period = StudyPeriod[[Study]][Period]
  Trt = toupper(letters[Trt_number]) # get trt label assignments for each record
  Cow = rep(1:ncows[Study],each=ntrt)  # and cow labels for each record
  Cow_effects = rnorm(sd=sqrt(sigmacow),ncows[Study])
  Effect_Cow = Cow_effects[Cow]
  Trt_Study_effects =  studytrtint[[Study]][Trt_number] # and the trt*study effects
  residuals = rnorm(sd = sqrt(sigmae),ntrt*ncows[Study]) # and finally the residuals
  FCM = round(overallmean + Effect_Trt + Effect_Period  + Effect_Cow + Studyeffects[Study] + Trt_Study_effects + residuals,2)  # generate the record and round to 2nd decimal
  Study_data[[Study]] = data.frame(Study,Trt,Cow,Period,FCM)  #store it.
}  

overallcrossoverdata = bind_rows(Study_data)  # the entire data set
cols = c("Study","Trt","Cow","Period")
overallcrossoverdata[cols] <- lapply(overallcrossoverdata [cols], factor)
str(overallcrossoverdata)
write_csv(overallcrossoverdata,"overallcrossoverdata.csv")

# save the true study by treatment effects for later.
study_trt = unlist(studytrtint)
Study_Trtlabels = expand.grid(trtlabels,1:nStudy)
Study_Trtlabels = paste0(Study_Trtlabels$Var1,':',Study_Trtlabels$Var2)
study_trt = data.frame(Study_Trtlabels,study_trt)

##########################################################################################################
######################## END OF SIMULATION ###############################################################
##########################################################################################################


#' 
#' 
#' # Analysis based on individual cow records across experiments
#' 
## ----IPD-----------------------------------------------------------------------------------------------------------------------------------------------------------------
library(lme4)  # mixed model procedure in R.
# NOTE THIS MODEL PERFECTLY MATCHES THE MODEL GENERATION PROCESS
overallcrossover_analysis = lmer(FCM~Trt+ (1|Period:Study) + (1|Cow:Study) + (1|Study/Trt),data=overallcrossoverdata)  # Fit Trt as Fixed, Study and Study*Trt as random
# Cows and Periods are nested within Study.  If Period represents a physiological stage (i.e. all cows start at roughly the same DIM), then treat as fixed.
# Might also fit Treatment by Period then as well.
summary(overallcrossover_analysis)
# fairly precise estimates of variance components....if you were able to analyze all of the data jointly.

library(emmeans)
emm_options(pbkrtest.limit = nrow(overallcrossoverdata)) 
lsmeans_Trt= emmeans(overallcrossover_analysis,"Trt")  
# so if this is taking too long, just set the emm_options(pbkrtest.limit = 500) or some reasonably low number 
(overallcrossover_means = summary(lsmeans_Trt))
(effectsize_IPD = pairs(lsmeans_Trt))  # overall effect sizes



#' 
#' 
#' The following is just for quality control checks
#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## THE FOLLOWING IS NOT REALLY IMPORTANT...JUST A QUALITY CONTROL CHECK
##Plot Estimated(Predicted) Study Effects versus True Study Effects
# Estimated ("Predicted") Study Effects
BLUP_Studyeffects = coef(overallcrossover_analysis)$Study[,1]-mean(coef(overallcrossover_analysis)$Study[,1])  # need to subtract the mean to express relative to zero
trt = data.frame(Studyeffects,BLUP_Studyeffects)
# Plot Estimated ("Predicted") Study Effects vs True Study Effects
ggplot(data=trt,aes(x=Studyeffects ,y=BLUP_Studyeffects)) + geom_point() + geom_abline(intercept=0,slope=1) +
  ggtitle("Predicted vs True Study Effects") +
  xlab("True effects") + 
  ylab("BLUP")

## THE FOLLOWING IS NOT REALLY IMPORTANT...JUST A QUALITY CONTROL CHECK
# Let's do the same thing for the study by treatment effects
#Estimated study by treatment effects 
BLUP_study_trt = coef(overallcrossover_analysis)$`Trt:Study`[,1]-mean(coef(overallcrossover_analysis)$`Trt:Study`[,1]) # need to subtract the mean to express relative to zero
BLUP_study_trt = data.frame(rownames(coef(overallcrossover_analysis)$`Trt:Study`),BLUP_study_trt)
colnames(BLUP_study_trt) = c("Study_Trtlabels","BLUP_study_trt")
study_trt = study_trt %>%
  left_join(BLUP_study_trt,by=c("Study_Trtlabels"))
# Plot Estimated ("Predicted") Study by Treatment Effects vs True Study by Treatment Effects
ggplot(study_trt,aes(x=study_trt ,y=BLUP_study_trt)) + geom_point() + geom_abline(intercept=0,slope=1) +
  ggtitle("Predicted vs True Study by Treatment Effects") +
  xlab("True effects") + 
  ylab("BLUP")

#' 
#' 
#' # Meta-analysis
#' 
#' Write a few functions that will help compute the study-specific statistics.
#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Let's conduct a separate mixed model analysis for each Study
# just like what we would anticipate for a meta-analysis
# Need to write some functions in order to do so

# need to create a separate function for model
# You could throw Trt*Period in there..maybe even consider Period as random
separate_model = function(df) {
  lmer(FCM~Trt+Period+(1|Cow),data=df)
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

#function to retrieve study-specific variance components
varcomps = function(model){
  varcomp = as.numeric(formatVC(VarCorr(model),comp="Variance")[,3])
  VClabel = (formatVC(VarCorr(model),comp="Variance")[,1])
  varcomp = data.frame(VClabel,varcomp) %>%
    pivot_wider(names_from=VClabel,values_from = varcomp)
  return(varcomp)
}




#' 
#' This is what it looks like for one study
## ----Study1--------------------------------------------------------------------------------------------------------------------------------------------------------------
Study1_analysis = lmer(FCM~Trt+Period+(1|Cow),data=filter(overallcrossoverdata,Study==1))
joint_tests(Study1_analysis)
summary(Study1_analysis)
varcomps(Study1_analysis)
lsmeans(Study1_analysis)
contrast(Study1_analysis)


#' 
#' 
## ----Study_spec_contrasts------------------------------------------------------------------------------------------------------------------------------------------------
# Let's save the contrasts (i.e. pairwise comparisons) from each study
Study_specific_contrasts = overallcrossoverdata %>%
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

#' 
#' Let's focus on the A vs B contrast for now.
#' 
#' We'll consider the common effects 
#' 
## ----CommonAvsB----------------------------------------------------------------------------------------------------------------------------------------------------------

(AB_contrast = filter(Study_specific_contrasts,contrast=="A - B"))
# The following analysis is rather common (and typically wrong!!!)
# COMMON EFFECTS MODEL 

#CE1
weighted.mean(AB_contrast$estimate,AB_contrast$weight)
(stderr = (sum(AB_contrast$weight))^(-1/2))
#compare this to what we determined from the overall analysis on the raw data
effectsize_IPD

# could go through the same exercise for the other comparisons as well.
#(AC_contrast = filter(Study_specific_contrasts,contrast=="A - C"))
#(BC_contrast = filter(Study_specific_contrasts,contrast=="B - C"))

#CE2
# # COMMON EFFECTS MODEL can also be conducted using a linear model
# glmmTMB is used because it allows one to fix the residual variance to 1 (just a "trick" to allow proper weightings by 1/SE^2 which already reflect residual var!)
library(glmmTMB)
# overall mean model
FMAvsB_1 <- glmmTMB(estimate ~ 1 ,
                weights = weight,
                family = gaussian,
                REML=TRUE,
                data = AB_contrast,
                start = list(betad = log(1)),  # fix residual variance = 1
                map = list(betad = factor(NA))
)
summary(FMAvsB_1)
# again gives the same wrong answer as above

#CE3
# Finally, onecan use metafor to get the same result:
library(metafor)
(res_EE_AvsB = rma.uni(yi=estimate,vi=sampvar,data=AB_contrast,method="EE"))  # fixed study effects
forest(res_EE_AvsB)

#' 
## ----MixedAvB------------------------------------------------------------------------------------------------------------------------------------------------------------
# an analysis that properly reflects heterogeneity between studies.
library(glmmTMB)
# random effects model
RM_1AvsB_1 <- glmmTMB(estimate ~ 1 + (1|Study),  # Study really involves both Study and Study*Treatment 
                weights = weight,
                family = gaussian,
                REML=TRUE,
                data = AB_contrast,
                start = list(betad = log(1)),  # fix residual variance = 1
                map = list(betad = factor(NA))
)
summary(RM_1AvsB_1)  
# compare this to the analysis using the raw data...point estimate is ok; standard error is off.

#using metafor
(res_REMLAvsB_1 = rma.uni(yi=estimate,vi=sampvar,data=AB_contrast,method="REML"))  # fixed study effects
# reminder as to what the corresponding estimate/se was from analyzing the raw data
effectsize_IPD



#' 
#' # Multivariate analysis using all contrast information
#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RM_mvcontrast <- glmmTMB(estimate ~ contrast -1 + (1|Study/contrast),  
                 weights = weight,
                 family = gaussian,
                 REML=TRUE,
                 data = Study_specific_contrasts,
                 start = list(betad = log(1)),  # fix residual variance = 1
                 map = list(betad = factor(NA))
)
summary(RM_mvcontrast)

# compare with analyses on individual data
effectsize_IPD


#' 
#' # Analysis based on Study specific means
#' (instead of contrasts)
#' 
## ----Studyspecmeans------------------------------------------------------------------------------------------------------------------------------------------------------
#lets save the study specific means

Study_specific_means = overallcrossoverdata %>%
  group_by(Study) %>%
  nest() %>%
  mutate(model = map(data,separate_model)) %>%
  mutate(lsmeans = map(model,lsmeans)) %>%
  mutate(nrec = as.numeric(map(data,nrow))) %>% 
  unnest(lsmeans)    %>%
  dplyr::select(-c(data,model)) %>%
  mutate(weight = (1/SE)^2) %>%     # these are not quite right!!
  mutate(sampvar = SE^2)          # these are not quite right!!
# but let's try them anyways

head(Study_specific_means)

#' 
## ----subopt--------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Let's use the suboptimal weights anyways

RM_mv_subopt <- glmmTMB(emmean ~ Trt   + (1|Study/Trt),  # 
                 weights = weight,
                 family = gaussian,
                 REML=TRUE,
                 data = Study_specific_means,
                 start = list(betad = log(1)),  # fix residual variance = 1
                 map = list(betad = factor(NA))
)
summary(RM_mv_subopt)
lsmeans_RM_mv_subopt= emmeans(RM_mv_subopt,"Trt")
(summary(lsmeans_RM_mv_subopt))
pairs(lsmeans_RM_mv_subopt)
# not too far off from what we had with raw data analysis (?)
effectsize_IPD  


#' 
#' # Let's look at study-specific variance components
#' 
## ----Study_specific_var--------------------------------------------------------------------------------------------------------------------------------------------------
Study_specific_var = overallcrossoverdata %>%
  group_by(Study) %>%
  nest() %>%
  mutate(model = map(data,separate_model)) %>%
  mutate(varcomp = map(model,varcomps))   %>%
  unnest(varcomp) %>%
  dplyr::select(-c(data,model)) 

boxplot(Study_specific_var[,-1] )

#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Let's look at the relationship between SEM, SED, and the variance components


SE_compare =  Study_specific_means %>%
    group_by(Study) %>%
    summarize(SEM = mean(SE))  %>%
    full_join(distinct(Study_specific_contrasts,SE),by="Study") %>%
    group_by(Study,SEM) %>%
    summarize(SED=mean(SE)) %>%
    full_join(Study_specific_var,by="Study") %>%
    mutate(SEM_altered = SED/sqrt(2)) %>%
    mutate(SEM_ratio = SEM/SEM_altered) %>%
    mutate(sd_ratio = sqrt((Cow+Residual)/Residual)) %>%
    mutate(weight_altered = (1/SEM_altered)^2) %>%   
    mutate(sampvar_altered = SEM_altered^2)    %>%
    mutate(weight_old = (1/SEM)^2) %>%  
    mutate(sampvar_old= SEM^2)  
head(SE_compare)

#' 
## ----weight_comparison---------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(SE_compare, aes(x=weight_old, y=weight_altered)) + 
  geom_point() 


#' 
#' 
## ----Means_analysis------------------------------------------------------------------------------------------------------------------------------------------------------
Study_specific_means2= Study_specific_means%>%
  mutate(SEM_altered = SE/sqrt(2)) %>%
  mutate(weight_altered = (1/SEM_altered)^2) %>%  # What i use in the mixed model
  mutate(sampvar_altered = SEM_altered^2)  

RM_mv2arm <- glmmTMB(emmean ~ Trt   + (1|Study/Trt),   
                 weights = weight_altered,
                 family = gaussian,
                 REML=TRUE,
                 data = Study_specific_means2,
                 start = list(betad = log(1)),  # fix residual variance = 1
                 map = list(betad = factor(NA))
)
summary(RM_mv2arm)
lsmeans_Trt_mv2arm= emmeans(RM_mv2arm,"Trt")
( summary(lsmeans_Trt_mv2arm))
pairs(lsmeans_Trt_mv2arm)

#' # # A network meta-analysis
#' # use same specifications as above:
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
trtchoose = vector(mode="list",length= nStudy)  # 

for (Study in 1:nStudy) {
  #trtchoose[[Study]] = toupper(letters[sample(1:ntrt,(ntrt-1))])
  trtchoose[[Study]] = sample(1:ntrt,(ntrt-1))
  }
Trt=unlist(trtchoose)
Study = rep(1:nStudy,each=(ntrt-1))
IB_design = data.frame(Study,Trt) %>%
  mutate(Trt = toupper(letters[Trt])) %>%
  mutate_all(as.factor)

xtabs(~Study+Trt,data=IB_design)

IB_Data = IB_design %>%
  left_join(overallcrossoverdata,by=c("Study","Trt"))


#' 
#' ## analysis of IB data
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# NOTE THIS MODEL PERFECTLY MATCHES THE MODEL GENERATION PROCESS
IB_analysis = lmer(FCM~Trt+ (1|Period:Study) + (1|Cow:Study) + (1|Study/Trt),data=IB_Data)  # Fit Trt as Fixed, Study and Study*Trt as random
# Cows and Periods are nested within Study.  If Period represents a physiological stage (i.e. all cows start at roughly the same DIM), then treat as fixed.
# Might also fit Treatment by Period then as well.
summary(IB_analysis)
lsmeans_TrtIB= emmeans(IB_analysis,"Trt")  
# so if this is taking too long, just set the emm_options(pbkrtest.limit = 500) or some reasonably low number 
(IB_analysis_means = summary(lsmeans_TrtIB))
(effectsize_IB = pairs(lsmeans_TrtIB))  # overall effect sizes


#' # Study specific contrasts
#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Study_specific_contrasts_IB = IB_Data %>%
  group_by(Study) %>%
  nest() %>%
  mutate(model = map(data,separate_model)) %>%
  mutate(effectsize = map(model,contrast)) %>%
  mutate(nrec = as.numeric(map(data,nrow))) %>% 
  unnest(effectsize)    %>%
  dplyr::select(-c(data,model)) %>%
  mutate(weight = (1/SE)^2) %>%  # What i use in the mixed model
  mutate(sampvar = SE^2) %>%         # What the r package metafor needs 
  mutate(contrast = factor(contrast))

#' 
#' # Naively only considering studies with A-B contrasts
#' 
#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Study_specific_contrasts_IB_AB = filter(Study_specific_contrasts_IB,contrast == 'A - B')
 
 # random effects model
RM_1_IB <- glmmTMB(estimate ~ 1 + (1|Study),   
                weights = (1/SE)^2,
                family = gaussian,
                REML=TRUE,
                data = Study_specific_contrasts_IB_AB,
                start = list(betad = log(1)),  # fix residual variance = 1
                map = list(betad = factor(NA))
)
summary(RM_1_IB) 


#' # can also use study-specific means to do this
#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Study_specific_means_IB = IB_Data %>%
  group_by(Study) %>%
  nest() %>%
  mutate(model = map(data,separate_model)) %>%
  mutate(lsmeans = map(model,lsmeans)) %>%
  mutate(nrec = map(data,nrow)) %>% 
  unnest(lsmeans)    %>%
  dplyr::select(-c(data,model))

SEDuse = Study_specific_contrasts_IB %>%
  group_by(Study) %>%
  summarize(SED = mean(SE)) %>%
  mutate(SEM_alter = SED/sqrt(2)) %>%
  mutate(weight_alter = (1/SEM_alter^2)) %>%
  mutate(sampvar_alter = SEM_alter^2)

Study_specific_means_IB = Study_specific_means_IB %>%
  left_join(SEDuse,by="Study")

RM_arm_IB <- glmmTMB(emmean ~ Trt   + (1|Study/Trt),   
                 weights = weight_alter,
                 family = gaussian,
                 REML=TRUE,
                 data = Study_specific_means_IB,
                 start = list(betad = log(1)),  # fix residual variance = 1
                 map = list(betad = factor(NA))
)
summary(RM_arm_IB)
lsmeans_Trt_RM_arm_IB= emmeans(RM_arm_IB,"Trt")
( summary(lsmeans_Trt_RM_arm_IB))
pairs(lsmeans_Trt_RM_arm_IB)
effectsize_IB



#' 
#' 
#' ```

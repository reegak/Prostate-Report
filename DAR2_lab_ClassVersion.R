library(corrplot)
library(MASS)
library(rJava)
library(glmulti)
library(effsize)
#install.packages("corrplot")
#install.packages("rJava")
#install.packages("effsize")
#install.packages("glmulti")
##############
# Lab Tasks: 
#   Be prepared to discuss at end of class;
#   you will not hand anything in.
#   Use PK_EDALRModeleda_ClassVersion.R, PK_diagnostics_ClassVersion, and 
#     Prostate_diagnostics_ClassVersion.R as a reference for code.
##############


## Load in prostate data set and take initial looks
Prostate = read.table("prostate.txt", header=T)
dim(Prostate)
#[1] 380 8
sum(is.na(Prostate)) # missing data: 3 in race
Prostate<-na.omit(Prostate); dim(Prostate) # 377 8
Prostate$dcaps = as.factor(Prostate$dcaps)
Prostate$dpros = as.factor(Prostate$dpros)
attach(Prostate)


## Task 1: EDA
# A reminder: in this project, we will consider interactions, but not consider transformations.
# a) binary response, capsule: consider a simple table of counts
table(capsule)
# b) binary response against categorical covariates: consider tables 
#    (can perhaps consider side-by-side bar charts, though I think contingency tables are sufficient)
table(capsule, race)
table(capsule, dpros)

# c) binary response against continuous covariates: consider summary tables and box-plots
boxplot(psa, capsule)
boxplot(age, capsule)

# d) standardized mean difference: evaluation of relationships between covariates and the response
cohen.d(age, capsule)
cohen.d(psa, capsule)
cohen.d(gleason, capsule)

## Task 2: Model building
# Logistic regression models are fit using the glm function using the logit link:
#   e.g., glm(capsule~psa+gleason, family=binomial(link=logit), data=Prostate)
Prostate = Prostate[, -1]

fit = glm(capsule~.*., family = binomial(link=logit), data = Prostate)
fit1 = glm(capsule ~  dpros + psa + gleason  , family = binomial(link = logit), data = Prostate)

# a) Stepwise model selection: include interactions, consider stepAIC for first pass
stepAIC(fit)
stepAIC(fit1)

# b) Parsimonious model: perform backward selection via p-values,
#      identify a simpler model by being strict with interaction terms. 
#      Is it that much worse than best stepwise model?
summary(fit1)


## Task 3: Model evaluation
# Prostate_diagnostics_ClassVersion.R presents sample code for a given model via explanatory variable patterns (EVPs).
# The components include:
# a) Residual plots: use the examine.logistic.reg function provided
# b) Outlier detection: evaluate EVPs for potential outlying data points
# c) HL test of overall fit: use the HLtest.R function provided; see PK_diagnostics_ClassVerion.R

# The functions for residual analysis we proposed using are as follows:
one.fourth.root=function(x){
  x^0.25
}
# make sure examine.logistic.reg.R is in your working directory or you have the right path specified
source("examine.logistic.reg.R")

# Consider model of PSA, Gleason score, and Results of digital rectal exam
dat.glm <- glm(capsule ~ psa+gleason+dpros, family = binomial, data = Prostate)
dat.mf <- model.frame(dat.glm)
## Covariate pattern: too many EVPs!
w <- aggregate(formula = capsule ~ psa+gleason+dpros, data = Prostate, FUN = sum)
n <- aggregate(formula = capsule ~ psa+gleason+dpros, data = Prostate, FUN = length)
w.n <- data.frame(w, trials = n$capsule, prop = round(w$capsule/n$capsule,2))
dim(w.n)
#[1] 342 6

# Create EVPs by binning continuous covariates
g = 5 # number of categories
psa_interval = cut(psa, quantile(psa, 0:g/g), include.lowest = TRUE)  # Creates factor with levels 1,2,...,g
levels(psa_interval)

w <- aggregate(formula = capsule ~ psa_interval+gleason+dpros, data = Prostate, FUN = sum)
n <- aggregate(formula = capsule ~ psa_interval+gleason+dpros, data = Prostate, FUN = length)
w.n <- data.frame(w, trials = n$capsule, prop = round(w$capsule/n$capsule,2))
mod.prelim1 <- glm(formula = capsule/trials ~ psa_interval+gleason+dpros,
                   family = binomial(link = logit), data = w.n, weights = trials)
save1=examine.logistic.reg(mod.prelim1, identify.points=T, scale.n=one.fourth.root, scale.cookd=sqrt)
w.n.diag1=data.frame(w.n, pi.hat=round(save1$pi.hat, 2), std.res=round(save1$stand.resid, 2), 
                     cookd=round(save1$cookd, 2), h=round(save1$h, 2))
p=length(mod.prelim1$coef) # number of parameters in model (# coefficients)
ck.out=abs(w.n.diag1$std.res)>2 | w.n.diag1$cookd>4/nrow(w.n) | w.n.diag1$h > 3*p/nrow(w.n)
extract.EVPs=w.n.diag1[ck.out, ]
extract.EVPs


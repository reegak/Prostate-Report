# Read-in data
setwd("~/SchoolWork/Stat696/Prostate")
Prostate = read.table("prostate.txt", header=T)
dim(Prostate); names(Prostate)
#[1] 380 8
sum(is.na(Prostate)) # 3 in race
Prostate <- na.omit(Prostate); dim(Prostate) # 377 8
Prostate$dcaps = as.factor(Prostate$dcaps)
Prostate$dpros = as.factor(Prostate$dpros)
Prostate$race = as.factor(Prostate$race)
Prostate = Prostate[, -1]
attach(Prostate)
head(Prostate)
summary(Prostate)
library(MASS)
library(effsize)
library(VIF)
library(xtable)
library(corrplot)
# A reminder: in this project, we will consider interactions, but not consider transformations.
# a) binary response, capsule: consider a simple table of counts
table(capsule)

# b) binary response against categorical covariates: consider tables 
#    (can perhaps consider side-by-side bar charts, though I think contingency tables are sufficient)
table(capsule, race)
table(capsule, dpros)
table(capsule, age)
# c) binary response against continuous covariates: consider summary tables and box-plots
boxplot(psa ~ capsule, main = "PSA vs Capsule", xlab = "Non Penetration vs Penetration", ylab = "PSA (mg/mL)")
boxplot(age ~ capsule, main = "Age vs Capsule", xlab = "Non Peneration vs Peneration", ylab = "Age")
boxplot(gleason~capsule, main = "Gleason Score vs Capsule", ylab = "Level of abnormality of cells", xlab = "Non Peneration vs Peneration")

# d) standardized mean difference: evaluation of relationships between covariates and the response
cohen.d(age, capsule)
cohen.d(psa, capsule)
cohen.d(gleason, capsule)
cohen.d(as.numeric(dpros), capsule)

## Task 2: Model building
# Logistic regression models are fit using the glm function using the logit link:
#   e.g., glm(capsule~psa+gleason, family=binomial(link=logit), data=Prostate)

fit = glm(capsule~ .*., family = binomial(link=logit), data = Prostate)
fit1 = glm(capsule ~  dpros + psa + gleason, family = binomial(link = logit), data = Prostate)

provars = data.frame(as.numeric(capsule), as.numeric(dpros), as.numeric(psa), as.numeric(gleason))
cpv = cor(provars)
corrplot(cpv)

# a) Stepwise model selection: include interactions, consider stepAIC for first pass
stepAIC(fit)
stepAIC(fit1)

# b) Parsimonious model: perform backward selection via p-values,
#      identify a simpler model by being strict with interaction terms. 
#      Is it that much worse than best stepwise model?
summary(fit)
summary(fit1)

## Task 3: Model evaluation
# Prostate_diagnostics_ClassVersion.R presents sample code for a given model via explanatory variable patterns (EVPs).
# The components include:
# a) Residual plots: use the examine.logistic.reg function provided
# b) Outlier detection: evaluate EVPs for potential outlying data points
# c) HL test of overall fit: use the HLtest.R function provided; see PK_diagnostics_ClassVerion.R

# Residual plots
# Load-in required functions
one.fourth.root=function(x){
  x^0.25
}
source("examine.logistic.reg.R")

# Consider model of PSA, Gleason score, and Detection of capsular involvement
dat.glm <- glm(capsule ~ psa+gleason+dcaps, family = binomial, data = Prostate)
dat.mf <- model.frame(dat.glm)
## Covariate pattern: too many EVPs!
w <- aggregate(formula = capsule ~ psa+gleason+dcaps, data = Prostate, FUN = sum)
n <- aggregate(formula = capsule ~ psa+gleason+dcaps, data = Prostate, FUN = length)
w.n <- data.frame(w, trials = n$capsule, prop = round(w$capsule/n$capsule,2))
dim(w.n)
#[1] 301 6

# Create EVPs by binning continuous covariates
g = 5 # number of categories
psa_interval = cut(psa, quantile(psa, 0:g/g), include.lowest = TRUE)  # Creates factor with levels 1,2,...,g
levels(psa_interval)

# Diagnostic plots
w <- aggregate(formula = capsule ~ psa_interval+gleason+dcaps, data = Prostate, FUN = sum)
n <- aggregate(formula = capsule ~ psa_interval+gleason+dcaps, data = Prostate, FUN = length)
w.n <- data.frame(w, trials = n$capsule, prop = round(w$capsule/n$capsule,2))
mod.prelim1 <- glm(formula = capsule/trials ~ psa_interval+gleason+dcaps,
                   family = binomial(link = logit), data = w.n, weights = trials)
save1 = examine.logistic.reg(mod.prelim1, identify.points=T, scale.n=one.fourth.root, scale.cookd=sqrt)

# Evaluation of EVPs for potential outlying sets of points
w.n.diag1 = data.frame(w.n, pi.hat=round(save1$pi.hat, 2), std.res=round(save1$stand.resid, 2), 
                     cookd=round(save1$cookd, 2), h=round(save1$h, 2))
p = length(mod.prelim1$coef) # number of parameters in model (# coefficients)
ck.out = abs(w.n.diag1$std.res)>2 | w.n.diag1$cookd>4/nrow(w.n) | w.n.diag1$h > 3*p/nrow(w.n)
extract.EVPs = w.n.diag1[ck.out, ]
extract.EVPs

# Note: EVPs purely for diagnostics, akin to histogram binning to assess
# distribution shape.  The analysis does not use this EVP binning.

betahat = formatC(signif(fit1$coeff, digits = 3), digits = 2, format = "f", flag = "#")
OR = formatC(signif(exp(fit1$coeff), digit = 3), digits = 2, format = "f", flag = "#")
SE = formatC(signif(summary(fit1)$coeff[,2], digits = 3), digits = 2, format = "f", flag = "#")
cibounds = formatC(signif(exp(confint(fit1)), digits = 3), digits = 2, format = "f", flag = "#")
pval = formatC(signif(summary(fit1)$ coeff[,4], digits = 4), format = "f", flag = "#")

x = cbind(betahat, OR, SE, pval, matrix(paste("(", cibounds[,1], "," , cibounds[,2], ")")))
colnames(x) = cbind("Coefficient", "Odds Ratio", "SE", "p-value", "95% CI on OR")
rownames(x) = cbind("intercept", "dpros2", "dpros3", "dpros4", "psa", "gleason")
inftable = xtable(x)
align(inftable) = "|l|ccccc|"
print(inftable)

###EDA###
table(race)
table(capsule)
table(capsule,race)
table(capsule,dcaps)
table(capsule,dpros)

prop.table(table(capsule))
prop.table(table(capsule,race))
prop.table(table(capsule,dcaps))
prop.table(table(capsule,dpros))


boxplot(age~capsule, data = Prostate)
boxplot(psa~capsule, data = Prostate)
boxplot(gleason~capsule, data = Prostate)
plot(dpros, capsule, data = Prostate)

cohen.d(age, capsule)
cohen.d(psa, capsule)
cohen.d(capsule, race)
cohen.d(gleason, capsule)

t.test(age[capsule==0], age[capsule==1])
t.test(psa[capsule==0], psa[capsule==1])
t.test(gleason[capsule==0], gleason[capsule==1])
table(age)
table(race)
## prediction
new.data = data.frame(dpros = as.factor(3), psa = 14.1, gleason = 7)
new.data2 = data.frame(dpros = as.factor(2), psa = 30, gleason = 8)
new.data3 = data.frame(dpros = as.factor(4), psa = 25, gleason = 9)
new.data4 = data.frame(dpros=as.factor(2), psa=15.25, gleason=6)
predict(fit1, new.data, interval = "prediction", type = "response")
predict(fit1, new.data2, interval = "prediction", type = "response")
predict(fit1, new.data3, interval = "prediction", type = "response")
predict.glm(fit1, new.data4, interval="prediction", type = "response")
chisq.test(x = dpros, y = capsule)

#diagnostics
one.fourth.root=function(x){
  x^0.25
}
source("examine.logistic.reg.R")
# Consider model of PSA, Gleason score, and Results of digital rectal exam
dat.glm <- glm(capsule ~ psa+gleason+dpros, family = binomial, data = Prostate)
dat.mf <- model.frame(dat.glm)
## Covariate pattern: too many EVPs!
w <- aggregate(formula = capsule ~ psa+gleason+dpros, data = Prostate, FUN = sum)
n <- aggregate(formula = capsule ~ psa+gleason+dpros, data = Prostate, FUN = length)
w.n <- data.frame(w, trials = n$capsule, prop = round(w$capsule/n$capsule,2))


# Create EVPs by binning continuous covariates
g = 5 # number of categories
psa_interval = cut(psa, quantile(psa, 0:g/g), include.lowest = TRUE)  # Creates factor with levels 1,2,...,g

# Diagnostic plots
v <- aggregate(formula = capsule ~ psa_interval+gleason+dpros, data = Prostate, FUN = sum)
m <- aggregate(formula = capsule ~ psa_interval+gleason+dpros, data = Prostate, FUN = length)
v.m <- data.frame(v, trials = m$capsule, prop = round(v$capsule/m$capsule,2))
mod.prelim <- glm(formula = capsule/trials ~ psa_interval+gleason+dpros,
                   family = binomial(link = logit), data = v.m, weights = trials)
save = examine.logistic.reg(mod.prelim, identify.points=T, scale.n=one.fourth.root, scale.cookd=sqrt)

v.m.diag = data.frame(v.m, pi.hat=round(save$pi.hat, 2), std.res=round(save$stand.resid, 2), 
                       cookd=round(save$cookd, 2), h=round(save$h, 2))
p = length(mod.prelim$coef) # number of parameters in model (# coefficients)
ck.out = abs(v.m.diag$std.res)>2 | v.m.diag$cookd>4/nrow(v.m) | v.m.diag$h > 3*p/nrow(v.m)
extract.EVP = v.m.diag[ck.out, ]
extract.EVP
source("C:/Users/Kelso Quan/Documents/SchoolWork/Stat696/Prostate/HLTest.R")
HLTest(mod.prelim1, 4)



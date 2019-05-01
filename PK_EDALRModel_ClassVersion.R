# Logistic regression: EDA and modeling
# NFL Placekicking data

# Load libraries
library(corrplot)
library(MASS) # for stepAIC
library(effsize) # for Cohen's d
library(glmulti) # for model selection

# Load in data and take initial looks
# 1995 NFL season
PKdata = read.table("placekick.mb.txt", header=T)
dim(PKdata)
#[1] 1438   13
numkicks = nrow(PKdata)
attach(PKdata)

# Functions to take a look at the data
# type: outdoor (1) vs. dome (0)
head(PKdata)
tail(PKdata)
sum(is.na(PKdata)) # 0: no missing data
names(PKdata)
summary(PKdata) # notice the categorical variables


################################
# EDA: Univariate analysis
################################
           
## Contingency tables for categorical variables

# response variable: good
table(good)
#good
#   0    1 
# 165 1273 
mean(good) # 88% success rate
# table of percentages
signif(prop.table(table(good)), digits = 3)

# contingency tables of good vs. categorical covariates
table(PAT, good)
#   good
#PAT   0   1
#  0 151 511
#  1  14 762

# Box-plots of continuous covariates over good
# Note that distance is measured in discrete yard units
plot(factor(good), distance, col="cyan", varwidth=T,
     ylab="Distance",xlab="Kick Success")


## Correlation plot
PKcor = cor(PKdata)
corrplot(PKcor)
# In NFL 1995, all extra points are initially attempted from the 20 yard line!
table(PAT, distance)
# Otherwise, literature does say all that really matters is distance!

## Two-sample t-tests for continuous covariates
t.test(distance[good==0], distance[good==1])

#Welch Two Sample t-test
#
#data:  distance[good == 0] and distance[good == 1]
#t = 16.904, df = 198.94, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  12.78909 16.16689
#sample estimates:
#  mean of x mean of y 
#40.35152  25.87353 

## Univariate logistic regression
fit.dist = glm(good~distance, family=binomial(link=logit), data=PKdata)
summary(fit.dist)
#Call:
#glm(formula = good ~ distance, family = binomial(link = logit), 
#    data = PKdata)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-2.7367   0.2449   0.2449   0.3830   1.6087  
#
#Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  5.781823   0.322548   17.93   <2e-16 ***
#distance    -0.114496   0.008272  -13.84   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#    Null deviance: 1024.77  on 1437  degrees of freedom
#Residual deviance:  787.06  on 1436  degrees of freedom
#AIC: 791.06
#
#Number of Fisher Scoring iterations: 6

## Chi-square test for categorical variables
chisq.test(x=home, y=good)
#Pearson's Chi-squared test with Yates' continuity correction
#data:  home and good
#X-squared = 0.32505, df = 1, p-value = 0.5686


## Cohen's d; standardized mean difference
# looking for a value above 0.2 (small), 0.5 (medium), and 0.8 (large)
cohen.d(distance, factor(good))  # d = -1.4
cohen.d(PAT, factor(good))  # d = 1.7
cohen.d(change, factor(good)) # d = 0.4


####################
# Logistic regression inferences
# Let's take a step back to review odds ratios
####################

## Odds ratio examples
fit.elap = glm(good~elap30, family=binomial(link=logit), data=PKdata)
summary(fit.elap)
exp(fit.elap$coeff) # odds ratio is e^{\beta}
# OR = 1.02 
# so for every minute further from end of half, 2% increase in odds of a successful kick

# What about when coefficient < 1
exp(fit.dist$coeff)  # 0.89
# OR = 1/0.89 = 1.12 so for every yard closer, 12% increase in odds of a successful kick
# Alternatively, negate the coefficient in the exponential
exp(-10*fit.dist$coeff) 
# OR = 3.14 so odds of success increases by 3.14 times for every 10 yards closer 
confint(fit.dist) # (-0.13, -0.10)
exp(confint(fit.dist))  # (0.88, 0.91)
1/exp(confint(fit.dist)) # (1.10, 1.14)

# What about for a categorical predictor?
fit.type = glm(good~type, family=binomial(link=logit), data=PKdata)
summary(fit.type)
exp(fit.type$coeff) 
# OR = 1.26, so 26% increase in odds of a successful kick in a domed stadium


####################
# Model building
####################

fit.all=glm(good~., family=binomial(link=logit), data=PKdata)
summary(fit.all)
#Call:
#glm(formula = good ~ ., family = binomial(link = logit), data = placekick.mb)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-2.9778   0.1679   0.2029   0.4598   1.5234  
#
#Coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  5.545e+00  9.243e-01   5.999 1.99e-09 ***
#week        -3.921e-02  2.371e-02  -1.654  0.09822 .  
#distance    -8.827e-02  1.133e-02  -7.788 6.81e-15 ***
#altitude     6.821e-05  1.080e-04   0.631  0.52779    
#home         2.216e-01  1.875e-01   1.182  0.23726    
#type         2.107e-01  3.176e-01   0.663  0.50710    
#precip      -3.371e-01  4.652e-01  -0.725  0.46866    
#wind        -7.844e-01  3.694e-01  -2.124  0.03370 *  
#change      -3.172e-01  1.966e-01  -1.614  0.10656    
#elap30       3.500e-03  1.053e-02   0.333  0.73948    
#PAT          1.090e+00  3.699e-01   2.946  0.00322 ** 
#field       -2.573e-01  2.693e-01  -0.956  0.33924    
#temp72      -9.041e-03  8.269e-03  -1.093  0.27427    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#    Null deviance: 1024.77  on 1437  degrees of freedom
#Residual deviance:  764.54  on 1425  degrees of freedom
#AIC: 790.54
#
#Number of Fisher Scoring iterations: 6


## Logistic regression modeling

# Consider all possible main effects
# Consider the following interactions:
#  distance with altitude, precip, wind, change, elap30, PAT, field, temp72
#    since each of these could have a larger effect on longer kicks than shorter kicks.
#  home with wind
#    homefield advantage: know the wind patterns better than visiting kicker
#  precip with type, field, temp72
#    precipitation can not affect kicks in a dome, has a bigger effect on grass, 
#    and is worse in colder conditions
fit.all = glm(good~week+distance+altitude+home+type+precip+wind+change+elap30+PAT+field+temp72+
                distance:altitude+distance:precip+distance:wind+distance:change+distance:elap30+distance:PAT+
                distance:field+distance:temp72+home:wind+precip:type+precip:field+field:temp72, 
              family=binomial(link=logit), data=PKdata)
stepAIC(fit.all)
#Call:  glm(formula = good ~ distance + wind + change + PAT + distance:wind + 
#    distance:PAT, family = binomial(link = logit), data = PKdata)
#
#Coefficients:
#  (Intercept)       distance           wind         change            PAT  distance:wind  
#      4.49640       -0.08069        2.92477       -0.33205        6.71190       -0.09183  
# distance:PAT  
#     -0.27171  
#
#Degrees of Freedom: 1437 Total (i.e. Null);  1431 Residual
#Null Deviance:	    1025 
#Residual Deviance: 759.7 	AIC: 773.7

# Model selected by stepAIC
final.fit = glm(good~distance+wind+change+PAT+wind:distance+PAT:distance, 
                family=binomial(link=logit), data=PKdata)
summary(final.fit)


###############################################
# Alternative function in R: glmulti
# method "h" is exhaustive search
search.aicc=glmulti(good~., data=PKdata, fitfunction="glm", level=1, method="h", crit="aic", 
                    family=binomial(link="logit"))
print(search.aicc)
aa=weightable(search.aicc)  # contains top 100 models with aicc
cbind(model=aa[1:5, 1], round(aa[1:5,2:3], 3))  # the best 5 models
# Full data set
#                                                      model    aic weights
#1                 good ~ 1 + distance + wind + change + PAT 780.764   0.031
#2          good ~ 1 + week + distance + wind + change + PAT 781.080   0.027
#3                 good ~ 1 + week + distance + change + PAT 781.123   0.026
#4                        good ~ 1 + distance + change + PAT 781.259   0.025
#5 good ~ 1 + week + distance + wind + change + PAT + temp72 781.368   0.023
#
# The values under weights are the model weights (also called "Akaike weight"). 
# From an information-theoretic perspective, the Akaike weight for a particular model 
# can be regarded as the probability that the model is the best model 
# (in a Kullback-Leibler sense of minimizing the loss of information when 
# approximating full reality by a fitted model). So, while the "best" model 
# has the highest weight/probability, its weight in this example is not 
# substantially larger than that of the second model (and also the third, fourth, 
# and so on). So, we shouldn't be all too certain here that we have really found 
# the best model. Several models are almost equally plausible (in other examples, 
# one or two models may carry most of the weight, but not here). 

fit.me = glm(good~distance+wind+change+PAT, 
             family=binomial(link=logit), data=PKdata)
summary(fit.me)

# Function to compute AIC, AICc, and BIC
AICc=function(object){
  n=length(object$y)
  r=length(object$coef)
  AICc=AIC(object)+2*r*(r+1)/(n-r-1) # finite sample size corrected AIC
  list(AIC=AIC(object), AICc=AICc, BIC=BIC(object))
}

AICc(fit.me)
# Full data
#$AIC
#[1] 780.7644
#$AICc
#[1] 780.8063
#$BIC
#[1] 807.1195
AICc(final.fit)
#$AIC
#[1] 773.7168
#$AICc
#[1] 773.7951
#$BIC
#[1] 810.6139


##############################################################
# Considering transformations of the explanatory variables
# Explanatory variable patterns (EVPs)
##############################################################

# Straight scatterplot will obviously be useless!
plot(x=distance, y=good, xlab="Kick Distance", ylab="Success of Kick")

# Transform data into binomial (convert to explanatory variable pattern (EVP) form)
# Compute for each covariate pattern the number of successful kicks (w) and the number of kicks (n)
w <- aggregate(formula = good ~ distance, data = PKdata, FUN = sum)
n <- aggregate(formula = good ~ distance, data = PKdata, FUN = length)
w.n <- data.frame(w, trials = n$good, prop = round(w$good/n$good,2))
head(w.n)
nrow(w.n)  # Number of EVPs (covariate patterns)
sum(w.n$trials)  # Number of observations
plot(x=w.n$distance, y=w.n$prop, xlab="Distance", ylab="Success Probability")
lines(lowess(x=w.n$distance, y=w.n$prop), col="blue")
abline(lm(w.n$prop~w.n$distance), col="green")

# EVPs are tougher for "truly" continuous variables as there are very few repeat values 
# At the least, need to bin the data manually
aggregate(good~elap30, data=PKdata, FUN=length) # not enough replicates in elap30
elap30cat = cut(elap30, c(-1,0,5,10,15,20,25,30)) # manually choose cut-points
summary(elap30cat)
fit.elap30cat = glm(good~elap30cat, family=binomial(link=logit))
plot(x=c(0, 2.5,7.5,12.5,17.5,22.5,27.5), y=fit.elap30cat$coefficients, 
     xlab="Intervals for # Minutes before end of half", ylab="Regression coefficient estimates")
lines(lowess(x=c(0, 2.5,7.5,12.5,17.5,22.5,27.5), y=fit.elap30cat$coefficients), col="blue")
abline(lm(fit.elap30cat$coefficients~c(0, 2.5,7.5,12.5,17.5,22.5,27.5)), col="green")
# Check out cubic fit just on elap30
e2 = elap30^2; e3 = elap30^3
summary(glm(good~elap30+e2+e3, family=binomial(link=logit)))

elap_lin = glm(good~elap30, family=binomial(link=logit))
elap_Q = glm(good~elap30+e2, family=binomial(link=logit)) 
# neither term significant in quadratic model, so probably would not even try the cubic!
elap_cubic = glm(good~elap30+e2+e3, family=binomial(link=logit))
AICc(elap_lin)
AICc(elap_Q)
AICc(elap_cubic)

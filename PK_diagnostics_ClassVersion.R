# Logistic regression: Model diagnostics, more on EVPs
# NFL Placekicking data

# Load libraries
library(binom)

# Load in data and take initial looks
PKdata = read.table("placekick.mb.txt", header=T)
#load("PKttdata.Rdata")
dim(PKdata)
#[1] 1438   13
n = dim(PKdata)[1] # sample size

# Reminder of variables
# response variable is 'good'
attach(PKdata)


# Show straight residual plot will not work in logistic regression
# Recall final model
fit.prelim1=glm(good~
                  distance + wind + change + PAT + 
                  distance:wind + distance:PAT, 
                  family=binomial(link=logit), data=PKdata)  
std.resid=rstandard(model=fit.prelim1, type="pearson")
plot(x = distance, y = std.resid, ylim = c(min(-3, std.resid), max(3, std.resid)), 
     ylab = "Standardized Pearson residuals", xlab = "Distance")
abline(h = 0, lty = "dashed", col = "gray")

# Transform data into binomial (convert to explanatory variable pattern (EVP) form)
# Compute for each covariate pattern the number of successful kicks (w) and the number of kicks (n)
w <- aggregate(formula = good ~ distance + wind + change + PAT, data = PKdata, FUN = sum)
n <- aggregate(formula = good ~ distance + wind + change + PAT, data = PKdata, FUN = length)
w.n <- data.frame(w, trials = n$good, prop = round(w$good/n$good,2))
head(w.n) # View the EVPs
nrow(w.n) # Number of EVPs (covariate patterns)
sum(w.n$trials)  # Number of observations
# preliminary model: weighted logistic regression of EVP proportions on explanatory vars in final model
mod.prelim1 <- glm(formula = good/trials ~ 
                     distance + wind + change + PAT + distance:wind + distance:PAT,
                     family = binomial(link = logit), data = w.n, weights = trials)
# weighted logistic regression results same as straight logistic regression results! 
round(summary(mod.prelim1)$coefficients, digits = 4)
round(summary(fit.prelim1)$coefficients, digits = 4)

# Now can do residual plots!
stand.resid <- rstandard(model = mod.prelim1, type = "pearson")
plot(x = w.n$distance, y = stand.resid, ylim = c(min(-3, stand.resid), max(3, stand.resid)), 
     ylab = "Standardized Pearson residuals", xlab = "Distance")
abline(h = c(3, 2, 0, -2, -3), lty = "dotted", col = "blue")
# Add a smooth to the residual plot, have to be careful of the weights
ord.dist <- order(w.n$distance)
smooth.stand <- loess(formula = stand.resid ~ distance, data = w.n, weights = trials)
lines(x = w.n$distance[ord.dist], y = predict(smooth.stand)[ord.dist], lty = "solid", col = "red")
#lines(lowess(x = w.n$distance, y = stand.resid), col="blue")


# The text "Analysis of Categorical Data with R" by Bilder and Loughin provides a 
# nice diagnostics R graphics function.  It requires a (1/4)th root function so that is included here.
one.fourth.root=function(x){
x^0.25
}
source(file="Examine.logistic.reg.R")

# Take a look at diagnostics for final PK data model.
save1=examine.logistic.reg(mod.prelim1, identify.points=T, scale.n=one.fourth.root, scale.cookd=sqrt)
names(save1)
# Let's examine the outliers/influential points idenitfied by the plots
# Store prediction, residual, Cook's D, and leverage values
w.n.diag1=data.frame(w.n, pi.hat=round(save1$pi.hat, 2), std.res=round(save1$stand.resid, 2), 
                     cookd=round(save1$cookd, 2), h=round(save1$h, 2))
p=length(mod.prelim1$coef) # number of parameters in model (# coefficients)
# Set thresholds for residuals, Cooks' distance, and leverage 
# Recall now working with EVPs, not complete data set.  
# So denominators here are the #EVPs M=124, not n=1438
ck.out=abs(w.n.diag1$std.res)>2 | w.n.diag1$cookd>4/nrow(w.n) | w.n.diag1$h > 3*p/nrow(w.n)
extract.EVPs=w.n.diag1[ck.out, ]
extract.EVPs[order(extract.EVPs$distance),]  # order by distance
#    distance wind change PAT good trials prop pi.hat std.res cookd    h
#60        18    0      1   0    1      2 0.50   0.94   -2.57  0.01 0.01
#117       20    0      0   1  605    614 0.99   0.98    0.32  0.06 0.81
#121       20    1      0   1   42     42 1.00   0.99    0.52  0.01 0.19
#123       20    0      1   1   94     97 0.97   0.98   -0.75  0.02 0.23
#101       25    1      1   0    1      2 0.50   0.94   -2.73  0.06 0.05
#119       29    0      0   1    0      1 0.00   0.73   -1.76  0.07 0.13
#120       30    0      0   1    3      4 0.75   0.65    0.83  0.31 0.76
#103       31    1      1   0    0      1 0.00   0.85   -2.43  0.03 0.03
#15        32    0      0   0   12     18 0.67   0.87   -2.67  0.06 0.05
#48        32    1      0   0    0      1 0.00   0.87   -2.62  0.02 0.02
#87        45    0      1   0    1      5 0.20   0.63   -2.03  0.02 0.03
#55        50    1      0   0    1      1 1.00   0.23    1.90  0.04 0.07

# EVP #117 has the largest h (leverage). This is due to large # of trials for this EVP. 
#   Since the corresponding Cooks' distance and residual is not overly large, no concern. 
# EVP #120 has the largest Cook's distance. 
# Notice the nonstandard distance (20 yards) for EVPs 119 and 120 for a PAT kick. 
# There are a total of 5 EVPs containing a total of 13 observations that are non-20 yard PATs. 
# To complete the analysis, we would exclude these kicks and re-run the models and diagnostics.


# HLtest.R is an R function for the classic and best known Hosmer and Lemeshow goodness-of-fit test
source("HLtest.R")
HL = HLTest(mod.prelim1, 10)  # 10 groups by default
# HL test output: Y0 are successes, Y1 are failures
cbind(HL$observed, round(HL$expect, digits=1)) # Observed and Expected table as illustration
HL # HL test results
#	Hosmer and Lemeshow gof test with 10 bins
#
#data:  mod.prelim1
#X2 = 2.3415, df = 8, p-value = 0.9687

library(ISLR)
library(ROCR)

# This assignment is based on ISLR Chapter 5, exercise 5
# Predicting whether a person will default based on annual income, monthly credit card balance, and student status
dim(Default)
n = dim(Default)[1]  # sample size
names(Default)

# Q1: How to predict a default case?

# Split data into training and testing sets and illustrate sensitivity and specificity
c = 0.5 # probability cutoff for predicting a default
set.seed(1)  # set the random number generator seed
p = 0.5 # percentage split for training and testing sets
train = sample(n, p*n)  # random sample percentage out of n; this creates the index list
length(train); head(train) # take a look at the 'train' variable index
head(Default)
# Fit logistic regression model on the training subset
fit_train = glm(default~balance+income, family=binomial(link=logit), data=Default, subset = train)
summary(fit_train)

# Make predictions on the test subset
test_probs = predict.glm(fit_train, Default, type="response")[-train]
test_class = (test_probs > c)
head(test_probs)
tail(test_class)
# Q2: How to evaluate default predictions?
# Q3: What measures to use to assess predictive performance?

table(test_class, Default$default[-train], dnn=c("Predicted", "Truth"))  # cross-classification accuracy
sum(as.numeric(test_class) == (as.numeric(Default$default[-train])-1))/(n-n*p) # accuracy
# Rows are predicted default status, columns are true default status
#test_class   No  Yes
#FALSE 4805  115
#TRUE    28   52
#sensitivity = 52/(115+52)  # percentage of true defaulters correctly identified: 31%  (true positive rate)
#specificity = 4805/(4805+28) # percentage of non-defaulters correctly identified: 99% (true negative rate)


# Q4: How can we potentially predict better?

# Function that will compute sensitivity and specificity at any given cutoff
se.sp <- function (cutoff, pred){
  sens <- performance(pred,"sens")
  spec <- performance(pred,"spec")
  num.cutoff <- which.min(abs(sens@x.values[[1]] - cutoff))
  return(list(Cutoff=sens@x.values[[1]][num.cutoff],
              Sensitivity=sens@y.values[[1]][num.cutoff], 
              Specificity=spec@y.values[[1]][num.cutoff]))
}
pred = prediction(test_probs, Default[-train,"default"])  # Prediction class from ROCR package
se.sp(0.5, pred) # Sensitivity and specificity at 0.5 cutoff

# What if increase cut-off?
se.sp(0.8, pred)
# What if decrease cut-off?
se.sp(0.3, pred)


# ROC curve illustration, presenting cutpoints on a colorized curve
pred = prediction(test_probs, Default[-train,"default"])
ROCRperf = performance(pred, "tpr", "fpr")
plot(ROCRperf, colorize=TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7),
    xlab="False positive rate (1-specificity)", ylab="True positive rate (sensitivity)") 
abline(a=0, b= 1, col="gray")
arrows(0.3,0.2,0.2,0.2,length=0.1)
text(0.3,0.2,"coin flip; no better than random", cex=1, pos=4)
# ROC curve for reporting
plot(ROCRperf)
abline(a=0, b= 1, col="gray")
# calculating AUC
auc1 = performance(pred,"auc")
# convert S4 class to vector
auc1 = unlist(slot(auc1, "y.values"))
text(0.7,0.2,paste("AUC is", signif(auc1, digits=2)), col="blue")


# Optimal cutpoint by optimizing over sensitivity and specificity simultaneously
# criterion: fpr^2 + fnr^2
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = x^2 + (1-y)^2  # minimize over fpr and fnr
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
# Compute optimal cutoff for the training set
# Present sensitivity and specificity for that optimal cutoff
test_probs = predict.glm(fit_train, Default, type="response")[train]
pred = prediction(test_probs, Default[train,"default"])
roc.perf = performance(pred, measure="tpr", x.measure="fpr")
print(opt.cut(roc.perf, pred))

# Evaluate optimal cut-point on testing set
optcut = opt.cut(roc.perf,pred)[3]
test_probs = predict.glm(fit_train, Default, type="response")[-train]
pred = prediction(test_probs, Default[-train,"default"])
roc.perf = performance(pred, measure="tpr", x.measure="fpr")
se.sp(optcut, pred)

# OptimalCutpoints package
library(OptimalCutpoints)
help(optimal.cutpoints)

# Q5: How can we create a risk score?   
Question 1 

library(ggplot2) #ggplot2 helps us to plot 2 regression lines in the same plot
set.seed(301) #this makes sure that the same results can be reproduced on any computer 
x <- c(1:99) #generating the initial set of n random numbers. Here n = 99
y <- 0.15*x+5 + rnorm(99,0,1) #generating the y variables based on the x random numbers
data <- data.frame(x, y) #creating the data set of x and y
reg1 <- lm(y~x, data=data) #creating the first regression line for our data
summary(reg1) #summary of reg1
#outlier <-  data.frame(x=c(0), y=c(-1500)) #generating an outlier to the data
updatedData <- rbind(data, c(120,-200)) #adding the outlier to the existing dataset to generate a new one
reg2 <- lm(y~x, data=updatedData) #creating the second regression line using the new dataset
summary(reg2) #summary of reg2
plot(updatedData, ylim=c(-250, 80), xlab="X Values", ylab = "Y Values", main = "Regression lines with different data points") #plotting the dataset
abline(reg1, col='red', lwd=1) #plotting the first reg line
abline(reg2, col='blue', lwd=1) #plotting the second reg line

Question 2

#loading prereqs
library("arm")
data(lalonde)
set.seed(20)
attach(data_control_only)

data_control_only <- lalonde[which(lalonde$treat==0),] #getting the data that only has control and not treatment
#creating the linear model for the control group data on the given variables
linear_model <- lm(re78 ~ age + educ + re74 + re75 + educ*re74 + educ*re75 + age*re74 + age*re75 + re74*re75, data_control_only)

#simulations
reps = 10000
simulation <- sim(linear_model, reps)
simulation


#medians
med_edu <- median(educ)
med_re74 <- median(re74)
med_re75 <- median(re75)

#creating a loop
mat1 <- matrix(, nrow = 39, ncol = 3)
store_prediction <- vector()
for (age in c(17:55)){
    Xs <- c(1, age, med_edu, med_re74, med_re75, med_edu*med_re74, med_edu*med_re75, age*med_re74, age*med_re75, med_re74*med_re75) 
    for (i in 1:reps){
      prediction <- sum((simulation@coef[i,])*Xs) + rnorm(1, 0, simulation@sigma[i])
      store_prediction <- c(store_prediction, prediction )
    
    }
    mat1[age-16,] <- c(age,  quantile(store_prediction, probs = c(0.025,0.975))) 
}
#store_prediction
#mat1

plot(x = c(1:100), y = c(1:100), type = "n", xlim = c(17,55), ylim = c(-10000,20000), 
     main = "Predicting re78 by Age using median variables", xlab = "Age", 
     ylab = "Predicted re78")

for (age in 17:55) {
  segments(
    x0 = age,
    y0 = mat1[age-16, 2],
    x1 = age,
    y1 = mat1[age-16, 3],
    lwd = 2)
}

#quantiles
quan_edu <- quantile(educ, probs = 0.9)
quan_re74 <- quantile(re74, probs = 0.9)
quan_re75 <- quantile(re75, probs = 0.9)

#creating a loop
mat2 <- matrix(, nrow = 39, ncol = 3)
store_prediction2 <- vector()
for (age in c(17:55)){
  Xz <- c(1, age, quan_edu, quan_re74, quan_re75, quan_edu*quan_re74, quan_edu*quan_re75, age*quan_re74, age*quan_re75, quan_re74*quan_re75) 
  for (i in 1:reps){
    prediction2 <- sum((simulation@coef[i,])*Xz) + rnorm(1, 0, simulation@sigma[i])
    store_prediction2 <- c(store_prediction2, prediction2 )
    
  }
  mat2[age-16,] <- c(age,  quantile(store_prediction2, probs = c(0.025,0.975))) 
}
store_prediction2
mat2

plot(x = c(1:100), y = c(1:100), type = "n", xlim = c(17,55), ylim = c(-10000,20000), 
     main = "Predicting re78 by Age using quantile variables", xlab = "Age", 
     ylab = "Predicted re78")

for (age in 17:55) {
  segments(
    x0 = age,
    y0 = mat2[age-16, 2],
    x1 = age,
    y1 = mat2[age-16, 3],
    lwd = 2)
}


table5 <- data.frame(mat1[,1], med_edu, med_re74, med_re75, mat1[,2], mat1[,3], quan_edu, quan_re74, quan_re75, mat2[,2], mat2[,3])
colnames(table5) <- c("Age", "Median Education", "Median Re74", "Median Re75", "Prediction Bound A", "Prediction Bound B", "Quantile Education", "Quantile Re74", "Quantile Re75", "Prediction A using Quantile", "Prediction B using Quantile")
 
Question 3

#seeting up vars
library(foreign)
datas <- read.dta("nsw.dta")
set.seed(123)

#model generation
model1 <- lm(re78 ~ treat, datas)
model1
summary(model1)
ci_using_formula <- confint(model1, level = 0.95)[2,]
ci_using_formula

#ci using bootstrap
store <- vector()
for (i in 1:10000){
  indexing <- sample(1:nrow(datas), nrow(datas), replace = T)
  bootstrap <- datas[indexing, ]
  bootstrapping_ci <- lm(re78 ~ treat, bootstrap)
  store <- c(store, bootstrapping_ci$coef[2])
}
ci_using_bootstrap <- quantile(store, probs = c(0.025,0.975))
ci_using_bootstrap

#histogram
hiss <- hist(store, xlab = "Coefficients of treatment", main = "Bootstrap Sampling Results", col = "green")
hiss

#table with relevant results
summary_table <- data.frame(ci_using_formula, ci_using_bootstrap)
summary_table

Question 4

rss <- vector() # residual sum of squares
tss <- vector() # total sum of squares

rsquared = function(pred_y, actual_y){ #creating a function for R Squared
  rss <- sum((pred_y - actual_y) ^ 2)
  tss <- sum((actual_y - mean(actual_y)) ^ 2)
  rsq <- 1 - rss/tss
  rsq
}

pred_y = predict(model1)
actual_y = datas$re78
function_rsq_value <- rsquared(pred_y, actual_y)
function_rsq_value
summary(model1)$r.squared

#comparing the R^2 values from the summary of the model1 from question 3, 
#vs the rsquared value from question 4 to check the accuracy
all.equal(summary(model1)$r.squared, function_rsq_value, tolerance = 1.5e-14)

Question 5

#setting up vars
library(foreign)
datas <- read.dta("nsw.dta")
set.seed(123)

#building the models
glm_fitting<- glm(treat~age+education+black+hispanic+married+nodegree+re75, data=datas, family=binomial)
glm_probabilities <- predict(glm_fitting, type="response")
View(glm_probabilities)

control_set <- vector()
treatment_set <- vector()
for (i in 1:nrow(datas)){
  if (datas$treat[i] == 0) {
    control_set <- c(control_set, glm_probabilities[i])
  }
  else {
    treatment_set <- c(treatment_set, glm_probabilities[i])
  }
}

histogram_a <- hist(treatment_set, col = "red", xlab = "Probabilities", ylab = "Treatment Group", main = "Treatment group's estimated probabilities")
histogram_b <- hist(control_set, col = "blue", xlab = "Probabilities", ylab = "Control Group", main = "Control group's estimated probabilities")


View(treatment_set)


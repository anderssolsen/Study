x <- c(1, 2) #vector
A <- matrix(c(1, 2, 3, 4), nrow = 2) #matrix
solve(A, x) #solves Ax=b (x is b?)

#random numbers and linear model
x <- runif(20, 0, 10)
y <- 2 * x - 3 + rnorm(20, sd = 1.5)
coef(lm(y ~ x))


par(mfrow = c(1, 2))
plot(x, y)
abline(lm(y ~ x), col = "red", lwd = 3)
plot(lm(y ~ x), which = 1, lwd = 3)

#Run all commands saved in a text file with:
source("myfile.R")


x <- c(1:5, -5:-1, 10, 20)
x[5:8]
x[x < 0]
x[-c(5, 9)] #all except

A <- matrix((-4):5, nrow = 2, ncol = 5)
A[A < 0]
A[A < 0] <- 0
A[2, ]
A[, c(2, 4)]

# Function with default values, note default values!
weighted.ave <- function(x, w = rep(1, length(x))) {
s1 <- sum(x * w)
s2 <- sum(w)
return(s1/s2)
}
weighted.ave(c(0, 1, 2, 3, 4), c(7, 3, 10, 17, 21))

getwd()
setwd()
save.image(file = "myStuff.RData") #save all command window and variables
load("myStuff.RData")

# load data organized with header. This is called a dataframe
# A dataframe is almost like a matrix
# Columns extracted by using $, like myData$x (for column 2 called x)
myData <- read.table("datafile.tab", header = TRUE)

x <- rep(1:5, each = 3)
f <- factor(x) #decribe categories, only 5 here
is.factor(x)
is.factor(f)

#General linear models lm()
#Generalized linear models glm()
#Mixed effects linear models lme() (in package nlme)
#Mixed effects linear models lmer() (in package lme4)
#Generalized Mixed effects linear models glmer() (in package lme4).

fit <- lm(y ~ x + f + g:h + k:z)
# this corresponds to y_i=ax_i+b(f_i)+c(g_i,h_i)+d(k_i)z_i+e_i
# which are factors? Interactions between two factors is different than factor and covariate
# Interactions are specified with f:g
# writing f*g is the same as writing f+g+f:g
# Adding -1 to the formula gets rid of the common intercept mu

x <- rep(1:5, each = 3)
y <- sin(2 * x) + rnorm(15, sd = 0.1)
f <- factor(x)
plot(x, y)
abline(lm(y ~ x), col = "red", lwd = 3)
coef(lm(y ~ x))
coef(lm(y ~ f - 1))

#Built-in distributions
#d<name>(x) Density function
#p<name>(x) Cumulated density function (probability  x )
#q<name>(p) Quantile (the point x where the probability  x is p)
#r<name>(n) Simulate random numbers from the distribution

#A basic R installation has: beta, binom, cauchy, chisq, exp, f,
#gamma, geom, hyper, logis, multinom, nbinom, norm, pois,
#signrank, t, tukey, unif, weibull, wilcox
# To get muktivariate normal (installs dmvnorm, pmvnorm, qmvnorm, rmvnorm)
install.packages("mvtnorm")
library(mvtnorm)

# Basic control flow:
A <- matrix((-4):5, nrow = 2, ncol = 5)
S <- 0 #Sum all positive elements poorly
for (i in 1:nrow(A)) {
  for (j in 1:ncol(A)) {
    if (A[i, j] > 0) {
      S <- S + A[i, j]
    }
  }
}
sum(A[A>0]) #much easier way
#In general: built in functions such as 
#(sum(), min(), max(), which(),
#which.min(), rowMeans(), colMeans(), rowSums(), colSums()

#Built in apply functions:
apply() #Use a function over one index (or more) of a matrix (array)
tapply() #Use a function within a number of groups
lapply() #Use a function for each element in a list
apply(A,2,sd) #could also be mean

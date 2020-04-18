library('splines')


# Reading data from file
worms <-read.table("c:nndatannworms.txt",header=T,row.names=1)
# uually we want to use attach to make variables accessible
# and use names to see a list of variable names. 
summary{worms}

worms[Area>3 & Slope<4,] #selects only rows where this is true
worms[order(worms[,1]),1:6] #sort by area (column one)
worms[rev(order(worms[,4])),c(4,6)] #descending order, only 2 columns as output

#Specification of models.
# y~x or y~1+x is like y_i = mu + a*x_i + e_i
# y ~ -1+x implies no intercept
# y ~ a specifies y_ij = a_j + e_ij if a is a factor
# y ~ a1+a2 additive two-sided model
# y ~ a1+a2 + a1_a2 two sided model with interaction (or a1*a2)
# a1*a2*a3 = (1+a1):(1+a2):(1+a2) = 1+a1+a2+a3+a1:a2+a1:a3+a2:a3+a1:a2:a3
# furthe3, (a1+a2+a3)^3 is same as a1*a2*a3
# Whereas (a1+a2+a3)^2 is same as a1*a2*a3 - a1:a2:a3

# we can write log(y)~sqrt(x)
# BUT we should not use ^, /, * on continuous variables, use I(x1*x2) instead

# some models
# summary(lm(..))
# anova(lm(..))
# anova(fit.H0,fit.HA) #specific hypotheses

# Tree-way anove without three-way interaction
y~N*P*K - N:P:K
# Analysis of covariance, a common slope for y and x but two intercepts, one for each gender
y~x+gender
#Split-plot anova, 3-way factorial setup but three different error variances
y~ a*b*c+Error(a/b/c)
#Including multiple polynomial regression
y~poly(x,2)+z
#Multiple regression (only with two-way interactions)
y~(x+z+w)^2
#Non-parametric model (y function of smoothed s and loss z)
y~ s(x)+lo(z)

#### NOT DONE

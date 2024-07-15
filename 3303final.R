knitr::opts_chunk$set(
	echo = FALSE,
	message = FALSE,
	warning = FALSE
)

library(rjags)
library(coda)
library(tidyverse)

a = read.csv('flu.txt', header = TRUE, sep = ',')

modelString <- "
model {
for(i in 1:N){
	y[i] ~ dbin(p[i],1);
	logit(p[i]) <- alpha[group[i]] + beta * x[i];
}

for(j in 1:g){
  alpha[j] ~ dnorm(alpha0, sigma0);
}

beta ~ dnorm(beta0, sigma10);
beta0 ~ dnorm(beta00, se1);
alpha0 ~ dnorm(alpha00, se0);
sigma0 ~ dgamma(0.01, 0.01);
sigma10 ~ dgamma(0.01, 0.01);
}
"
a$Country = factor(a$Country)
fit = glm(Infected ~ EZK , data = a, family = 'binomial')
summary(fit)

set.seed(1234)

datalist = list(y = a$Infected,
                group = a$Country,
                se1 = 0.097^2,
                se0 = 0.1385^2,
                g = length(unique(a$Country)),
                N = nrow(a),
                x = a$EZK,
                alpha00 = -0.86882,
                beta00 = 1.73273)

adaptSteps = 1000
 burnInSteps = 2000
 nChains = 4
 numSavedSteps = 8000
 thinSteps = 2
 nIter = ceiling((numSavedSteps * thinSteps)/nChains)
 
 
m1 = jags.model(textConnection(modelString), data = datalist,n.chains = nChains,n.adapt = adaptSteps)
fit.samples = coda.samples(m1,c("alpha","beta"),n.iter = nIter,n.chains = nChains, n.adapt = adaptSteps)

summary(fit.samples)


par(mar=c(2,2,2,2))
plot(window(fit.samples,start= burnInSteps+adaptSteps))

set.seed(1234)
n = 10000
x = rnorm(n, 1.6605, 2.0077) 
p1 = 1 - mean(1/(1 + exp(-x)))

n = 10000
x = rnorm(n, -0.4182, 0.2020)
p2 = mean(1/(1 + exp(-x)))
re = c(p1, p2)
names(re) = c('NonInf-EZK:1', 'Inf-EZK:0')
re

rm(list=ls())
source("00_functions.txt")
set.seed(2024)
library(multiColl)

##################################################################################

observaciones = seq(15,150,5)
media = seq(-19,19,6) # avoid zero
desv.tip = 100
total = length(observaciones)*length(media)*length(desv.tip)
total

min.CV = array(, total)
max.VIF = array(, total)
NC = array(, total)
incremento.beta = array(, total)
incremento.beta.sinCTE = array(, total)
observacion = array(, total)
medias = array(, total)
desviaciones = array(, total)
tol = 0.01
i = 1
for (n in observaciones){
  cte = array(1, n)
  for (mu in media) {
    for (dv in desv.tip){
      x1 = rnorm(n, mu, dv)
      x2 = rnorm(n, mu, dv)
      x = cbind(x1, x2)
      X = cbind(cte, x1, x2)
      #
      min.CV[i] = min(CVs(X))
      #
      max.VIF[i] = max(VIF(X))
      #
      NC[i] = CN(X)
      #
      Y = 3 + 0.4*x1 + 10*x2 + rnorm(n,0,4)
      reg = lm(Y~x)
      beta = as.matrix(reg$coefficients)
      #
      x1.p = perturb(x1, sample(-5:5,1), sample(1:5,1), tol)
      x2.p = perturb(x2, sample(-5:5,1), sample(1:5,1), tol)
      x.p = cbind(x1.p, x2.p)
      #
      reg.p = lm(Y~x.p)
      beta.p = as.matrix(reg.p$coefficients)
      #
      incremento.beta[i] = (norm(beta-beta.p,"2")/norm(beta,"2"))*100
      incremento.beta.sinCTE[i] = (norm(beta[-1]-beta.p[-1],"2")/norm(beta[-1],"2"))*100
      #
      observacion[i] = n
      medias[i] = mu
      desviaciones[i] = dv
      #
      i = i + 1
      print(i)
    }
  }
}



min(min.CV) # multicollinearity not-essential not troubling
max(min.CV)
min(max.VIF)  # multicollinearity essential not troubling
max(max.VIF)
min(NC)  # multicollinearity not troubling
max(NC)

min(incremento.beta) 
mean(incremento.beta) 
median(incremento.beta) 
max(incremento.beta) # significant changes in estimates 

min(incremento.beta.sinCTE)
mean(incremento.beta.sinCTE)
median(incremento.beta.sinCTE)
max(incremento.beta.sinCTE) # not significant changes in estimates 

#results = data.frame(min.CV, max.VIF, NC, incremento.beta, incremento.beta.sinCTE, observacion, medias, desviaciones)
#results
#write(t(results), "04_results.txt", ncolumns = 8, sep=";")
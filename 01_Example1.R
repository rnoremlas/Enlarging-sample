rm(list = ls())
source("00_functions.txt")

##################################################################################
# data reading

datos = read.table("01_Wissel.txt", header = T, sep=";")
head(datos)
datos = datos[,-5] # delete 'CP'
datos = datos[,-1] # delete 't'
head(datos)

y = datos[,1] # dependent variable
x = as.matrix(datos[,-1])  # independents variable

  # estimation of the original model

  reg = lm(y~x) 
  summary(reg) # Model 1 (second column of Table 1)

  # multicollinearity detection

  library(rvif)
  cte = rep(1,length(y))
  X = cbind(cte, x)
  RVIF(X, l_u=TRUE, l=40, intercept=TRUE, graf=TRUE) # Figure 1
  RdetR(X)
  VIF(X)
  CNs(X)
  CVs(X)

##################################################################################
# augmented model

  # calculation of the number of times to increase the sample using expression (10)
  
  n = length(y)
  k = dim(X)[2]
  hi = array(-1,k)
  for (i in 1:k){
    hi[i] = (1/n)*(((1.96/summary(reg)[[4]][i,3])^2)*(n-k)+k)
  }
  hi
  h = ceiling(max(hi[-1])) # delete intercept
  h 

  # augmented data

  data = datos
  for (i in 1:h){
    datos = rbind(datos,data)
  }
  #write(t(datos), "Wissel-augmented.txt", ncolumns=dim(datos)[2], sep=";")

  # estimation of the augmented model

  yA = datos[,1]
  xA = as.matrix(datos[,-1])

  regA = lm(yA~xA)
  summary(regA)  # Model 2 (third column of Table 1)

  # multicollinearity detection: does not change anything with respect to the lines 26-29
  
  cteA = rep(1,length(yA))
  XA = cbind(cteA, xA)
  RVIF(XA, l_u=TRUE, l=40, intercept=TRUE, graf=TRUE)
  RdetR(XA)
  VIF(XA)
  CNs(XA)
  CVs(XA)

##################################################################################
# disturbed augmented model

  # perturbation of the model data augmented
  
  tol = 0.01 # run from here (until line 125) for 0.01 = 1%, 0.02 = 2%, 0.05 = 5%
  media = 10
  dv = 10

  xA.p = matrix(-1, dim(xA)[1], dim(xA)[2])
  set.seed(2024)
  for (i in 1:dim(xA)[2]){
    xA.p[,i] = perturb(xA[,i], media, dv, tol)
  }

  # disturbed augmented data 
  
  datos.p = cbind(yA,xA.p)
  #write(t(datos.p), "Wissel-augmented-disturbed.txt", ncolumns=dim(datos.p)[2], sep=";")

  # estimation of the disturbed augmented model
  
  regA.p = lm(yA~xA.p)
  summary(regA.p) # Model 3 (fourth column of Table 1)
  
  # multicollinearity detection
  
  XA.p = cbind(cteA, xA.p)
  RVIF(XA.p, l_u=TRUE, l=40, intercept=TRUE, graf=TRUE) 
  RdetR(XA.p)
  VIF(XA.p)
  CNs(XA.p)
  CVs(XA.p)

  # percentage change in estimates

  betaA = as.double(regA$coefficients)
  betaA.p = as.double(regA.p$coefficients)
  (norm(betaA-betaA.p,"2")/norm(betaA,"2"))*100

##################################################################################
# perturbing and estimating the model repeated 1000 times

  iteraciones = 1000
  
  datosA = cbind(yA, xA) 
  disturbA = perturb_n.veces(datosA, iteraciones, media, dv, tol)
  mean(disturbA[,1]) # must coincide with 'tol'
  mean(disturbA[,2])
  min(disturbA[,2])
  max(disturbA[,2])
  var(disturbA[,2])
  sd(disturbA[,2])
  
##################################################################################
# ridge regression
  
  library(lrmest)
  
  rid(y~x, k=0) # coincide with Model 1 (second column of Table 3)
  VIFridge(x, start=0, leap=0, stop=0)
  
  # for k mitigating multicollinearity according to García, Salmerón and García (2019)
  
  detR = RdetR(X)[[2]]
  n = nrow(x)
  p = ncol(x)+1
  k_lineal =  0.01837*(1 - detR) - 0.00001262*n + 0.005678*p
  rid(y~x, k=k_lineal)  # third column of Table 3
  VIFridge(x, start=k_lineal, leap=0, stop=k_lineal)
  
  k_exp = 0.006639*exp(1 - detR) - 0.00001241*n + 0.005745*p # similar to the previous one
  rid(y~x, k=k_exp)
  VIFridge(x, start=k_exp, leap=0, stop=k_exp)
  
  k_sq =  0.7922*(1 - detR)^2 - 0.6901*(1 - detR) - 0.0000007567*n - 0.01081*p
  rid(y~x, k=k_sq)
  VIFridge(x, start=k_sq, leap=0, stop=k_sq)
  
  # for k calculated according to Hoerl, Kennard and Baldwin (1975)
  
  beta = reg$coef
  sigma2 = summary(reg)[[6]]^2
  k_HKB = p*(sigma2/crossprod(beta))
  rid(y~x, k=k_HKB) # fourth column of Table 3
  VIFridge(x, start=k_HKB, leap=0, stop=k_HKB)
  
  # plot of MSE
  
  k<-seq(0,1,0.1)
  mseval<-data.frame(rid(y~x, k))
  mseval
  k<-seq(0.1,1,0.1)
  plot(rid(y~x, k), main=c("Plot of MSE of ridge regression estimator"), type="b", lwd=3, col="Blue", xlab="k") # MSE plot example 1: Figure 2
  
  # for the k that minimizes the MSE
  
  k_min <- mseval[order(mseval[,2]),][1,1] # similar to the k_HKB
  rid(y~x, k=k_min)
  VIFridge(x, start=k_min, leap=0, stop=k_min)
  
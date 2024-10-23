rm(list = ls())
source("00_functions.txt")

##################################################################################
# data reading

datos = read.table("02_KleinGoldberger.txt", header = F, sep=";")
head(datos)
C = datos[,1] 
I = datos[,2]
InA = datos[,3]
IA = datos[,4]

# estimation of the original model

reg = lm(C~I+InA+IA)
summary(reg) # Model 1 (second column of Table 4)

# multicollinearity detection

  library(rvif)
  cte = rep(1,length(C))
  x = as.matrix(datos[,-1])
  X = cbind(cte, x)
  RVIF(X, l_u=TRUE, l=40, intercept=TRUE, graf=TRUE) # Figure 3
  RdetR(X)
  VIF(X)
  CNs(X)
  CVs(X)

##################################################################################
# augmented model

  # calculation of the number of times to increase the sample using expression (10)
  
  n = length(C)
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
  #write(t(datos), "KleinGoldberger-augmented.txt", ncolumns=dim(datos)[2], sep=";")
  
  # estimation of the augmented model
  
  yA = datos[,1]
  xA = as.matrix(datos[,-1])
  
  regA = lm(yA~xA)
  summary(regA)  # Model 2 (third column of Table 2)
  
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
  
  tol = 0.01 # run from here for 0.01 = 1%, 0.02 = 2%, 0.05 = 5%
  media = 10
  dv = 10
  
  xA.p = matrix(-1, dim(xA)[1], dim(xA)[2])
  set.seed(2024)
  for (i in 1:dim(xA)[2]){
    xA.p[,i] = perturb(xA[,i], media, dv, tol)
  }
  
  # disturbed augmented data 
  
  datos.p = cbind(yA,xA.p)
  #write(t(datos.p), "KleinGoldberger-augmented-disturbed.txt", ncolumns=dim(datos.p)[2], sep=";")
  
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
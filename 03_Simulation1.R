rm(list = ls())
source("00_functions.txt")

set.seed(2024)

##################################################################################
# Model 1: base data

obs = 50

x2 = round(rnorm(obs, 1100, 100), digits=2) # incomes
# x \in [a,b] \rightarrow y = A + ((B-A)/(b-a))*(x-a) \in [A, B]
a = min(x2)
b = max(x2)
A = 23
B = 32 
x3 = round(A + ((B-A)/(b-a))*(x2-a), digits=1) + sample(c(-0.1, 0, 0.1), obs, replace=T) # age
u = rnorm(obs, 0, 4)
y = round(3 + 0.4*x2 + 10*x3 + u, digits=1) # consumption

# OLS estimation

  reg = lm(y~x2+x3)
  summary(reg)

  cor.test(y, x2) # simple linear correlation significantly different from zero
  cor.test(y, x3) # simple linear correlation significantly different from zero

# multicollinearity detection

  library(rvif)
  cte = rep(1,length(y))
  X = cbind(cte, x2, x3)
  RVIF(X, l_u=TRUE, l=40, intercept=TRUE, graf=TRUE) # Figure 5
  RdetR(X)
  VIF(X)
  CNs(X)
  CVs(X)
  
# disturbed model  
  
  tol = 0.01 # 0.01 = 1%, 0.02 = 2%, 0.05 = 5%
  media = 10
  dv = 10
  
  Xp = matrix(-1, nrow(X), ncol(X)-1)
  for (i in 1:(ncol(X)-1)){
    Xp[,i] = perturb(X[,i+1], media, dv, tol)
  }
  
  # OLS estimation
  
  x2p = Xp[,1]
  x3p = Xp[,2]
  reg_p = lm(y~x2p+x3p)
  summary(reg_p)
  
  # percentage change in estimates
  
  beta = as.double(reg$coefficients)
  beta_p = as.double(reg_p$coefficients)
  variation1 = (norm(beta-beta_p,"2")/norm(beta,"2"))*100
  variation1
  variation1bis = (norm(beta[-1]-beta_p[-1],"2")/norm(beta[-1],"2"))*100
  variation1bis
  
  # model comparison
  
  library(stargazer)
  stargazer(reg, reg_p, type="text", title="Model comparison")
  
##################################################################################
# Model 2: "more of the same" is added to Model 1

obs1 = 50
  
  x2a = round(rnorm(obs1, 1100, 100), digits=2) # incomes
  x3a = round(A + ((B-A)/(b-a))*(x2a-a), digits=1) + sample(c(-0.1, 0, 0.1), obs1, replace=T) # age
  ua = rnorm(obs1, 0, 4)
  ya = round(3 + 0.4*x2a + 10*x3a + ua, digits=1) # consumption
  
  # OLS estimation
  
  x2a1 = c(x2, x2a)
  x3a1 = c(x3, x3a)
  ya1 = c(y, ya)
  
  reg_a1 = lm(ya1~x2a1+x3a1)
  summary(reg_a1)
  
  cor.test(ya1, x2a1) # simple linear correlation significantly different from zero
  cor.test(ya1, x3a1) # simple linear correlation significantly different from zero
  
  # multicollinearity detection
  
  cte = rep(1,length(ya1))
  Xa1 = cbind(cte, x2a1, x3a1)
  RVIF(Xa1, l_u=TRUE, l=40, intercept=TRUE, graf=TRUE) 
  RdetR(Xa1)
  VIF(Xa1)
  CNs(Xa1)
  CVs(Xa1)
  
  # disturbed model  
  
  tol = 0.01 # 0.01 = 1%, 0.02 = 2%, 0.05 = 5%
  media = 10
  dv = 10
  
  Xp = matrix(-1, nrow(Xa1), ncol(Xa1)-1)
  for (i in 1:(ncol(Xa1)-1)){
    Xp[,i] = perturb(Xa1[,i+1], media, dv, tol)
  }
  
  # OLS estimation
  
  x2a1p = Xp[,1]
  x3a1p = Xp[,2]
  reg_a1p = lm(ya1~x2a1p+x3a1p)
  summary(reg_a1p)
  
  # percentage change in estimates
  
  beta_a1 = as.double(reg_a1$coefficients)
  beta_a1p = as.double(reg_a1p$coefficients)
  variation2 = (norm(beta_a1-beta_a1p,"2")/norm(beta_a1,"2"))*100
  variation2
  variation2bis = (norm(beta_a1[-1]-beta_a1p[-1],"2")/norm(beta_a1[-1],"2"))*100
  variation2bis
  
  # model comparison
  
  stargazer(reg, reg_p, reg_a1, reg_a1p, type="text", title="Model comparison")  
  
##################################################################################
# Model 3: other information than the initial one is added to Model 1
  
  x2a = round(rnorm(obs1, 1800, 300), digits=2) # incomes
  x3a = round(rnorm(obs1, 45, 15), digits=1) # age
  ya = round(3 + 0.4*x2a + 10*x3a + ua, digits=1) # consumption
  
  # OLS estimation
  
  x2a2 = c(x2, x2a)
  x3a2 = c(x3, x3a)
  ya2 = c(y, ya)
  
  reg_a2 = lm(ya2~x2a2+x3a2)
  summary(reg_a2)
  
  cor.test(ya2, x2a2) # simple linear correlation significantly different from zero
  cor.test(ya2, x3a2) # simple linear correlation significantly different from zero
  
  # multicollinearity detection
  
  cte = rep(1,length(ya2))
  Xa2 = cbind(cte, x2a2, x3a2)
  RVIF(Xa2, l_u=TRUE, l=40, intercept=TRUE, graf=TRUE) # Figure 5
  RdetR(Xa2)
  VIF(Xa2)
  CNs(Xa2)
  CVs(Xa2)
  
  # disturbed model  
  
  tol = 0.01 # 0.01 = 1%, 0.02 = 2%, 0.05 = 5%
  media = 10
  dv = 10
  
  Xp = matrix(-1, nrow(Xa2), ncol(Xa2)-1)
  for (i in 1:(ncol(Xa2)-1)){
    Xp[,i] = perturb(Xa2[,i+1], media, dv, tol)
  }
  
  # OLS estimation
  
  x2a2p = Xp[,1]
  x3a2p = Xp[,2]
  reg_a2p = lm(ya2~x2a2p+x3a2p)
  summary(reg_a2p)
  
  # percentage change in estimates
  
  beta_a2 = as.double(reg_a2$coefficients)
  beta_a2p = as.double(reg_a2p$coefficients)
  variation3 = (norm(beta_a2-beta_a2p,"2")/norm(beta_a2,"2"))*100
  variation3
  variation3bis = (norm(beta_a2[-1]-beta_a2p[-1],"2")/norm(beta_a2[-1],"2"))*100
  variation3bis
  
  # model comparison
  
  stargazer(reg, reg_p, reg_a1, reg_a1p, reg_a2, reg_a2p, type="text", title="Model comparison")  # R2Model5 = summary(reg_a2)[[8]]
  
  c(variation1, variation2, variation3)
  c(variation1bis, variation2bis, variation3bis)
  c(max(VIF(X)), max(VIF(Xa1)), max(VIF(Xa2)))
  c(CN(X), CN(Xa1), CN(Xa2))
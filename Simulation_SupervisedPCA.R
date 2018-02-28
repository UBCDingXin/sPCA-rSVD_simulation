rm(list=ls())
library(MASS)
library(Matrix)
library(Rcpp)
library(elasticnet)
library(superpc)
library(glmnet)

setwd("C:/Users/dingx/ownCloud/Research/STAT548/Fourth paper_Gaby/Simulation")
source("Sparse PCA via Regularized Low Rank Matrix Approximation ----- rsvd.spca.Rfun.R")
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("sPCA-rSVD.cpp")

## Setups
nSim = 100 #number of simulations
n= 1000 # sample size
p = 100 # number of variables
ntest = n
prop_nonzero = 0.1; #proportion of non-zero entries in the parameter vector
varnum = prop_nonzero * p; #number of non-zero entries in the loading vector
CV = 5
varnum_grid = round(seq(from = 0.05 * p, to = 0.2 * p, length.out = 10))

MSE = array(0, c(nSim, 9)); #6 methods' MSE plus 3 CV versions

set.seed(1)
#sparse eigenvectors
vw1 = array(0,c(p,1))
nonzero_indx = sort(sample(x = 1:p, size = round(prop_nonzero * p), replace = FALSE));
vw1[nonzero_indx] = rnorm(n = length(nonzero_indx), mean = 10, sd = 5)
v1 = vw1 / (norm(vw1, type="2")) #first eigenvector
#singular values
C = c(200,seq(from = 1, to = 0.1, length.out = p-1))
#generate covariance matrix Sigma according to Gram-Schmidt orthogonalization
# GSO = gram( v.ini = vw1 , c.seq = C )
GSO = gramCpp( v_ini = vw1 , c_seq = C )
Sigma = GSO$s.mat
# v1 = GSO$v.mat[,1]

for (nn in 1:nSim){
  print(nn)
  set.seed(nn)
  ## ---------------------------------------------------------------------------------
  # data generation
  #generate X from normal with covariance Sigma
  X = mvrnorm(n = n+ntest, mu = rep(0, p), Sigma = Sigma )
  pc1 = X %*% v1; #first principal component
  #----------------------------------------------------
  #true model
  Y = 5 + pc1 + rnorm(n+ntest) #generate Y based on PC1
  #training set and testing set
  Xtrain = X[1:n,]; Xtest = X[(n+1):(n+ntest),];
  Ytrain = Y[1:n]; Ytest = Y[(n+1):(n+ntest)];
  
  ###############################################
  # 1. oracle OLS
  fit_lm = lm(Ytrain~pc1[1:n])
  Ypred = cbind(rep(1,ntest),pc1[(n+1):(n+ntest)]) %*% fit_lm$coefficients
  MSE[nn,1] = mean((Ypred-Ytest)^2)
  
  ###############################################
  # 2. Lasso regression
  fit_lass = cv.glmnet(Xtrain, Ytrain, family="gaussian", alpha = 1)
  Ypred2 = cbind(rep(1,ntest),Xtest) %*% coef(fit_lass,s=fit_lass$lambda.min)
  MSE[nn,2] = mean((Ypred2-Ytest)^2)
  
  ###############################################
  # 3. Supervised PCA
  #superpc
  featurenames1 <- paste("feature",as.character(1:n),sep="")
  featurenames2 <- paste("feature",as.character(1:ntest),sep="")
  data = list(x=t(Xtrain), y=Ytrain, featurenames = featurenames1)
  data.test = list(x=t(Xtest), y=Ytest, featurenames = featurenames2)
  train.obj = superpc.train(data, type="regression")
  cv.obj = superpc.cv(train.obj, data, n.components = 1)
  threshold = cv.obj$thresholds[colMeans(cv.obj$scor) == max(colMeans(cv.obj$scor))]
  fit.cts = superpc.predict(train.obj, data, data.test, 
                            threshold=threshold, n.components=1, prediction.type="continuous")
  MSE[nn,3] = mean((fit.cts$v.pred - Ytest)^2)
  
  ###############################################
  # 4. sPCA-rSVD-soft
  pca_soft_train = sPCArSVD(X = Xtrain, u = svd(Xtrain)$u[,1],
                            v = svd(Xtrain)$v[,1]*svd(Xtrain)$d[1],
                            varnum = varnum, type = 1);
  # pca_soft_train = sPCArSVDcv(X = Xtrain, varnum_grid = varnum_grid, type = 1, CV = CV)
  pc1_train = Xtrain %*% pca_soft_train$v / norm(pca_soft_train$v, "2")
  fit_tmp = lm(Ytrain~pc1_train) #regress Ytrain on pc1_train
  pca_soft_test = sPCArSVD(X = Xtest, u = svd(Xtest)$u[,1],
                           v = svd(Xtest)$v[,1]*svd(Xtest)$d[1],
                           varnum = varnum, type = 1)
  # pca_soft_test = sPCArSVDcv(X = Xtest, varnum_grid = varnum_grid, type = 1, CV = CV)
  pc1_test = Xtest %*% pca_soft_test$v / norm(pca_soft_test$v, "2")
  Ypred4 = cbind(rep(1,ntest),pc1_test) %*% fit_tmp$coefficients
  MSE[nn,4] = mean((Ypred4 - Ytest)^2)
  
  ###############################################
  # 5. sPCA-rSVD-hard
  pca_hard_train = sPCArSVD(X = Xtrain, u = svd(Xtrain)$u[,1],
                            v = svd(Xtrain)$v[,1]*svd(Xtrain)$d[1],
                            varnum = varnum, type = 2);
  # pca_hard_train = sPCArSVDcv(X = Xtrain, varnum_grid = varnum_grid, type = 2, CV = CV)
  pc1_train = Xtrain %*% pca_hard_train$v / norm(pca_hard_train$v, "2")
  fit_tmp = lm(Ytrain~pc1_train) #regress Ytrain on pc1_train
  pca_hard_test = sPCArSVD(X = Xtest, u = svd(Xtest)$u[,1],
                           v = svd(Xtest)$v[,1]*svd(Xtest)$d[1],
                           varnum = varnum, type = 2)
  # pca_hard_test = sPCArSVDcv(X = Xtest, varnum_grid = varnum_grid, type = 2, CV = CV)
  pc1_test = Xtest %*% pca_hard_test$v / norm(pca_hard_test$v, "2")
  Ypred5 = cbind(rep(1,ntest),pc1_test) %*% fit_tmp$coefficients
  MSE[nn,5] = mean((Ypred5 - Ytest)^2)

  ###############################################
  # 6. sPCA-rSVD-SCAD
  pca_SCAD_train = sPCArSVD(X = Xtrain, u = svd(Xtrain)$u[,1],
                            v = svd(Xtrain)$v[,1]*svd(Xtrain)$d[1],
                            varnum = varnum, type = 3);
  # pca_SCAD_train = sPCArSVDcv(X = Xtrain, varnum_grid = varnum_grid, type = 3, CV = CV)
  pc1_train = Xtrain %*% pca_SCAD_train$v / norm(pca_SCAD_train$v, "2")
  fit_tmp = lm(Ytrain~pc1_train) #regress Ytrain on pc1_train
  pca_SCAD_test = sPCArSVD(X = Xtest, u = svd(Xtest)$u[,1],
                           v = svd(Xtest)$v[,1]*svd(Xtest)$d[1],
                           varnum = varnum, type = 3)
  # pca_SCAD_test = sPCArSVDcv(X = Xtest, varnum_grid = varnum_grid, type = 3, CV = CV)
  pc1_test = Xtest %*% pca_SCAD_test$v / norm(pca_SCAD_test$v, "2")
  Ypred6 = cbind(rep(1,ntest),pc1_test) %*% fit_tmp$coefficients
  MSE[nn,6] = mean((Ypred6 - Ytest)^2)
  
  ###############################################
  # 7. sPCA-rSVD-soft CV
  pca_soft_train = sPCArSVDcv(X = Xtrain, varnum_grid = varnum_grid, type = 1, CV = CV)
  pc1_train = Xtrain %*% pca_soft_train$v / norm(pca_soft_train$v, "2")
  fit_tmp = lm(Ytrain~pc1_train) #regress Ytrain on pc1_train
  pca_soft_test = sPCArSVDcv(X = Xtest, varnum_grid = varnum_grid, type = 1, CV = CV)
  pc1_test = Xtest %*% pca_soft_test$v / norm(pca_soft_test$v, "2")
  Ypred4 = cbind(rep(1,ntest),pc1_test) %*% fit_tmp$coefficients
  MSE[nn,7] = mean((Ypred4 - Ytest)^2)
  
  ###############################################
  # 8. sPCA-rSVD-hard CV
  pca_hard_train = sPCArSVDcv(X = Xtrain, varnum_grid = varnum_grid, type = 2, CV = CV)
  pc1_train = Xtrain %*% pca_hard_train$v / norm(pca_hard_train$v, "2")
  fit_tmp = lm(Ytrain~pc1_train) #regress Ytrain on pc1_train
  pca_hard_test = sPCArSVDcv(X = Xtest, varnum_grid = varnum_grid, type = 2, CV = CV)
  pc1_test = Xtest %*% pca_hard_test$v / norm(pca_hard_test$v, "2")
  Ypred5 = cbind(rep(1,ntest),pc1_test) %*% fit_tmp$coefficients
  MSE[nn,8] = mean((Ypred5 - Ytest)^2)
  
  ###############################################
  # 9. sPCA-rSVD-SCAD CV
  pca_SCAD_train = sPCArSVDcv(X = Xtrain, varnum_grid = varnum_grid, type = 3, CV = CV)
  pc1_train = Xtrain %*% pca_SCAD_train$v / norm(pca_SCAD_train$v, "2")
  fit_tmp = lm(Ytrain~pc1_train) #regress Ytrain on pc1_train
  pca_SCAD_test = sPCArSVDcv(X = Xtest, varnum_grid = varnum_grid, type = 3, CV = CV)
  pc1_test = Xtest %*% pca_SCAD_test$v / norm(pca_SCAD_test$v, "2")
  Ypred6 = cbind(rep(1,ntest),pc1_test) %*% fit_tmp$coefficients
  MSE[nn,9] = mean((Ypred6 - Ytest)^2)
  
}#nSim simulations

# MSE[which(is.na(MSE),arr.ind=TRUE)] = colMeans(MSE,na.rm=TRUE)[which(is.na(MSE),arr.ind=TRUE)[2]]

output = cbind(colMeans(MSE,na.rm=TRUE), apply(MSE, MARGIN=2,FUN=min), 
               apply(MSE, MARGIN=2,FUN=max), apply(MSE,MARGIN=2,FUN=sd))
Sys.setenv(JAVA_HOME='C:/Program Files/Java/jdk1.8.0_131/jre')
library(XLConnect)
library(xlsx)
wb <- XLConnect::loadWorkbook("tmp.xlsx", create = TRUE)
sheetname="Sheet1"
XLConnect::writeWorksheet(wb,round(output,digits=2),sheetname,startRow = 1, startCol = 1, header = FALSE)
XLConnect::saveWorkbook(wb)

# save.image("superPC.Rdata")
























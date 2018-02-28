rm(list=ls())
library(MASS)
library(Matrix)
library(Rcpp)
library(elasticnet)
library(superpc)

# setwd()
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("sPCA-rSVD.cpp")

set.seed(548)
###################################################################
#Example 1
nSim = 100 #number of simulations
n= 300 # sample size
p = 10 # number of variables
q = 2 # number of sparse eigenvectors
CV = 5 #number of folders in cross-validation
varnum = 6 #number of non-zero entries in the loading vector
varnum_grid = 2:10
prop_outliers = 0.1
outliers_mean_prop = 5
outliers_sd = 1

#sparse eigenvectors
vw1 = c(rep(1,4), rep(0,4), rep(0.9,2))
vw2 = c(rep(0,4), rep(1,4), -0.3, 0.3)
v1 = vw1 / (norm(vw1, type="2"))
v2 = vw2 / (norm(vw2, type="2"))
eigenvectors = cbind(v1 , v2)

#threshold for simple thresholding method; 
thr_simple = 0.1

#singular values
C = diag(c(200, 100, 50, 50, 6, 5, 4, 3, 2, 1))

#generate covariance matrix Sigma according to Gram-Schmidt orthogonalization
# GSO = gram( v.ini = eigenvectors , c.seq = diag(C) )
GSO = gramCpp( v_ini = eigenvectors , c_seq = diag(C) )
Sigma = GSO$s.mat
# eigenvectors = GSO$v.mat[,1:2] #extract the first two eigenvectors after Gram-Schmidt orthogonalization
# eigenvectors[,2] = - eigenvectors[,2]

#initilization
#to store first two loading vectors (corresponding to the first PCs)
Load2_pca = Load2_soft = Load2_hard = Load2_SCAD = 
  Load2_Simple = Load2_SPCAk2 = Load2_SPCAk1 = array(0, c(p,2))
Load2_soft_cv = Load2_hard_cv = Load2_SCAD_cv = array(0, c(p,2))
#to store medain angle, correctness and incorrectness (true negative and false negative)
err2_pca = err2_soft = err2_hard = err2_SCAD = 
  err2_Simple = err2_SPCAk2 = err2_SPCAk1 = array(0, c(nSim,8)) 
err2_soft_cv = err2_hard_cv = err2_SCAD_cv = array(0, c(nSim,8)) 
#to store selected varnum in cross-validation
sel_varnum_soft = sel_varnum_hard = sel_varnum_SCAD = array(0, c(nSim,2))
#store a sparse outlier matrix
S = matrix(rep(0,n*p),nrow = n)


for (nn in 1:nSim){
  print(nn)
  if (nn>54){
    set.seed(nn+1000)
  }else{set.seed(nn+2000)}
  
  #keep generating data matrix X until it's full rank
  rankX = 0
  while(rankX != n && rankX != p){ #check whether X is full rank or not
    X_raw = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    rankX = rankMatrix(X_raw)[1]
  }
  
  #plus sparse outliers matrix
  S_nonzero_indx = sort(sample(x = 1:(n*p), size = round(prop_outliers * n*p), replace = FALSE)); 
  Outliers = rnorm(n = length(S_nonzero_indx),
                   mean = outliers_mean_prop*max(abs(c(max(X_raw),min(X_raw)))),
                   sd = outliers_sd);
  S[S_nonzero_indx] = Outliers;
  X_raw = X_raw + S; 

  # X = X_raw
  #center each column of X
  X_mean = colMeans(X_raw)
  X = X_raw - matrix(rep(X_mean, nrow(X_raw)), nrow = nrow(X_raw), byrow = TRUE)
  
  #1111111111111111111111111111111111111111111111111111111111111
  #standard PCA
  pca_std = prcomp(X)
  Load2_pca = pca_std$rotation[,1:2] #extract loading vectors
  #angle for v1
  err2_pca[nn, 1] = acos( (crossprod(eigenvectors[,1], Load2_pca[,1])) / 
                            (norm(eigenvectors[,1],"2")*norm(Load2_pca[,1],"2")) ) * 180 / pi
  #true negative rate for v1 (%)
  err2_pca[nn, 2] = CalCorrect(eigenvectors[,1], Load2_pca[,1])[1]
  #false negative rate for v1 (%)
  err2_pca[nn, 3] = CalCorrect(eigenvectors[,1], Load2_pca[,1])[2]
  #angle for v2
  err2_pca[nn, 4] = acos( (crossprod(eigenvectors[,2], Load2_pca[,2])) / 
                            (norm(eigenvectors[,2],"2")*norm(Load2_pca[,2],"2")) ) * 180 / pi
  #true negative rate for v2 (%)
  err2_pca[nn, 5] = CalCorrect(eigenvectors[,2], Load2_pca[,2])[1]
  #false negative rate for v2 (%)
  err2_pca[nn, 6] = CalCorrect(eigenvectors[,2], Load2_pca[,2])[2]
  # #l2-norm between v1 and hat(v1)
  # err2_pca[nn, 7] = norm(eigenvectors[,1]-Load2_pca[,1],"2")
  # #l2-norm between v2 and hat(v2)
  # err2_pca[nn, 8] = norm(eigenvectors[,2]-Load2_pca[,2],"2")
  
  #222222222222222222222222222222222222222222222222222222222222
  #SPCA_SOFT
  pca_soft1 = sPCArSVD(X = X, u = svd(X)$u[,1], 
                       v = svd(X)$v[,1]*svd(X)$d[1],
                       varnum = varnum, type = 1)
  R1 = X - tcrossprod(pca_soft1$u,pca_soft1$v) #residual matrix
  pca_soft2 = sPCArSVD(X = R1, u = svd(R1)$u[,1], 
                       v = svd(R1)$v[,1]*svd(R1)$d[1],
                       varnum = varnum, type = 1)
  Load2_soft[,1] = pca_soft1$v/norm(pca_soft1$v, "2") #extract loading vector 1
  Load2_soft[,2] = pca_soft2$v/norm(pca_soft2$v, "2") #extract loading vector 2
  tmp_v1 = Load2_soft[,1]
  tmp_v2 = Load2_soft[,2]
  #angle for v1
  err2_soft[nn, 1] = acos( (crossprod(eigenvectors[,1], tmp_v1)) / 
                             (norm(eigenvectors[,1],"2")*norm(tmp_v1,"2")) ) * 180 / pi
  #true negative rate for v1 (%)
  err2_soft[nn, 2] = CalCorrect(eigenvectors[,1], tmp_v1)[1]
  #false negative rate for v1 (%)
  err2_soft[nn, 3] = CalCorrect(eigenvectors[,1], tmp_v1)[2]
  #angle for v2
  err2_soft[nn, 4] = acos( (crossprod(eigenvectors[,2], tmp_v2)) / 
                             (norm(eigenvectors[,2],"2")*norm(tmp_v2,"2")) ) * 180 / pi
  #true negative rate for v2 (%)
  err2_soft[nn, 5] = CalCorrect(eigenvectors[,2], tmp_v2)[1]
  #false negative rate for v2 (%)
  err2_soft[nn, 6] = CalCorrect(eigenvectors[,2], tmp_v2)[2]
  # #l2-norm between v1 and hat(v1)
  # err2_soft[nn, 7] = norm(eigenvectors[,1]-tmp_v1,"2")
  # #l2-norm between v2 and hat(v2)
  # err2_soft[nn, 8] = norm(eigenvectors[,2]-tmp_v2,"2")
  
  
  #333333333333333333333333333333333333333333333333333333333333333333
  #SPCA_HARD
  pca_hard1 = sPCArSVD(X = X, u = svd(X)$u[,1], 
                       v = svd(X)$v[,1]*svd(X)$d[1],
                       varnum = varnum, type = 2)
  R1 = X - tcrossprod(pca_hard1$u,pca_hard1$v) #residual matrix
  pca_hard2 = sPCArSVD(X = R1, u = svd(R1)$u[,1], 
                       v = svd(R1)$v[,1]*svd(R1)$d[1],
                       varnum = varnum, type = 2)
  Load2_hard[,1] = pca_hard1$v/norm(pca_hard1$v, "2") #extract loading vector 1
  Load2_hard[,2] = pca_hard2$v/norm(pca_hard2$v, "2") #extract loading vector 2
  tmp_v1 = Load2_hard[,1]
  tmp_v2 = Load2_hard[,2]
  #angle for v1
  err2_hard[nn, 1] = acos( (crossprod(eigenvectors[,1], tmp_v1)) / 
                             (norm(eigenvectors[,1],"2")*norm(tmp_v1,"2")) ) * 180 / pi
  #true negative rate for v1 (%)
  err2_hard[nn, 2] = CalCorrect(eigenvectors[,1], tmp_v1)[1]
  #false negative rate for v1 (%)
  err2_hard[nn, 3] = CalCorrect(eigenvectors[,1], tmp_v1)[2]
  #angle for v2
  err2_hard[nn, 4] = acos( (crossprod(eigenvectors[,2], tmp_v2)) / 
                             (norm(eigenvectors[,2],"2")*norm(tmp_v2,"2")) ) * 180 / pi
  #true negative rate for v2 (%)
  err2_hard[nn, 5] = CalCorrect(eigenvectors[,2], tmp_v2)[1]
  #false negative rate for v2 (%)
  err2_hard[nn, 6] = CalCorrect(eigenvectors[,2], tmp_v2)[2]
  # #l2-norm between v1 and hat(v1)
  # err2_hard[nn, 7] = norm(eigenvectors[,1]-tmp_v1,"2")
  # #l2-norm between v2 and hat(v2)
  # err2_hard[nn, 8] = norm(eigenvectors[,2]-tmp_v2,"2")
  
  #4444444444444444444444444444444444444444444444444444444444444444444
  #SPCA_SCAD
  pca_SCAD1 = sPCArSVD(X = X, u = svd(X)$u[,1], 
                       v = svd(X)$v[,1]*svd(X)$d[1],
                       varnum = varnum, type = 3)
  R1 = X - tcrossprod(pca_SCAD1$u,pca_SCAD1$v) #residual matrix
  pca_SCAD2 = sPCArSVD(X = R1, u = svd(R1)$u[,1], 
                       v = svd(R1)$v[,1]*svd(R1)$d[1],
                       varnum = varnum, type = 3)
  Load2_SCAD[,1] = pca_SCAD1$v/norm(pca_SCAD1$v, "2") #extract loading vector 1
  Load2_SCAD[,2] = pca_SCAD2$v/norm(pca_SCAD2$v, "2") #extract loading vector 2
  tmp_v1 = Load2_SCAD[,1]
  tmp_v2 = Load2_SCAD[,2]
  #angle for v1
  err2_SCAD[nn, 1] = acos( (crossprod(eigenvectors[,1], tmp_v1)) / 
                             (norm(eigenvectors[,1],"2")*norm(tmp_v1,"2")) ) * 180 / pi
  #true negative rate for v1 (%)
  err2_SCAD[nn, 2] = CalCorrect(eigenvectors[,1], tmp_v1)[1]
  #false negative rate for v1 (%)
  err2_SCAD[nn, 3] = CalCorrect(eigenvectors[,1], tmp_v1)[2]
  #angle for v2
  err2_SCAD[nn, 4] = acos( (crossprod(eigenvectors[,2], tmp_v2)) / 
                             (norm(eigenvectors[,2],"2")*norm(tmp_v2,"2")) ) * 180 / pi
  #true negative rate for v2 (%)
  err2_SCAD[nn, 5] = CalCorrect(eigenvectors[,2], tmp_v2)[1]
  #false negative rate for v2 (%)
  err2_SCAD[nn, 6] = CalCorrect(eigenvectors[,2], tmp_v2)[2]
  # #l2-norm between v1 and hat(v1)
  # err2_SCAD[nn, 7] = norm(eigenvectors[,1]-tmp_v1,"2")
  # #l2-norm between v2 and hat(v2)
  # err2_SCAD[nn, 8] = norm(eigenvectors[,2]-tmp_v2,"2")
  
  
  #55555555555555555555555555555555555555555555555555555555555555555555
  #Simple thresholding method
  tmp_v1 = Load2_pca[,1]
  tmp_v2 = Load2_pca[,2]
  tmp_v1[abs(tmp_v1) < thr_simple] = 0
  tmp_v2[abs(tmp_v2) < thr_simple] = 0
  #angle for v1
  err2_Simple[nn, 1] = acos( (crossprod(eigenvectors[,1], tmp_v1)) / 
                               (norm(eigenvectors[,1],"2")*norm(tmp_v1,"2")) ) * 180 / pi
  #true negative rate for v1 (%)
  err2_Simple[nn, 2] = CalCorrect(eigenvectors[,1], tmp_v1)[1]
  #false negative rate for v1 (%)
  err2_Simple[nn, 3] = CalCorrect(eigenvectors[,1], tmp_v1)[2]
  #angle for v2
  err2_Simple[nn, 4] = acos( (crossprod(eigenvectors[,2], tmp_v2)) / 
                               (norm(eigenvectors[,2],"2")*norm(tmp_v2,"2")) ) * 180 / pi
  #true negative rate for v2 (%)
  err2_Simple[nn, 5] = CalCorrect(eigenvectors[,2], tmp_v2)[1]
  #false negative rate for v2 (%)
  err2_Simple[nn, 6] = CalCorrect(eigenvectors[,2], tmp_v2)[2]
  # #l2-norm between v1 and hat(v1)
  # err2_Simple[nn, 7] = norm(eigenvectors[,1]-tmp_v1,"2")
  # #l2-norm between v2 and hat(v2)
  # err2_Simple[nn, 8] = norm(eigenvectors[,2]-tmp_v2,"2")
  
  
  #66666666666666666666666666666666666666666666666666666666666666666666
  #SPCA_k2
  Load2_SPCAk2 = spca(X, K = 2, para = c(varnum, varnum), type="predictor", sparse = "varnum")$loadings #extract loading vectors
  tmp_v1 = Load2_SPCAk2[,1]
  tmp_v2 = Load2_SPCAk2[,2]
  #angle for v1
  err2_SPCAk2[nn, 1] = acos( (crossprod(eigenvectors[,1], tmp_v1)) / 
                               (norm(eigenvectors[,1],"2")*norm(tmp_v1,"2")) ) * 180 / pi
  #true negative rate for v1 (%)
  err2_SPCAk2[nn, 2] = CalCorrect(eigenvectors[,1], tmp_v1)[1]
  #false negative rate for v1 (%)
  err2_SPCAk2[nn, 3] = CalCorrect(eigenvectors[,1], tmp_v1)[2]
  #angle for v2
  err2_SPCAk2[nn, 4] = acos( (crossprod(eigenvectors[,2], tmp_v2)) / 
                               (norm(eigenvectors[,2],"2")*norm(tmp_v2,"2")) ) * 180 / pi
  #true negative rate for v2 (%)
  err2_SPCAk2[nn, 5] = CalCorrect(eigenvectors[,2], tmp_v2)[1]
  #false negative rate for v2 (%)
  err2_SPCAk2[nn, 6] = CalCorrect(eigenvectors[,2], tmp_v2)[2]
  # #l2-norm between v1 and hat(v1)
  # err2_SPCAk2[nn, 7] = norm(eigenvectors[,1]-tmp_v1,"2")
  # #l2-norm between v2 and hat(v2)
  # err2_SPCAk2[nn, 8] = norm(eigenvectors[,2]-tmp_v2,"2")
  
  
  #77777777777777777777777777777777777777777777777777777777777777777777
  #SPCA_k1
  Load2_SPCAk1 = spca(X, K = 1, para = varnum, type="predictor", sparse = "varnum")$loadings #extract loading vectors
  tmp_v1 = Load2_SPCAk1
  #angle for v1
  err2_SPCAk1[nn, 1] = acos( (crossprod(eigenvectors[,1], tmp_v1)) / 
                               (norm(eigenvectors[,1],"2")*norm(tmp_v1,"2")) ) * 180 / pi
  #true negative rate for v1 (%)
  err2_SPCAk1[nn, 2] = CalCorrect(eigenvectors[,1], tmp_v1)[1]
  #false negative rate for v1 (%)
  err2_SPCAk1[nn, 3] = CalCorrect(eigenvectors[,1], tmp_v1)[2]
  # #l2-norm between v1 and hat(v1)
  # err2_SPCAk1[nn, 7] = norm(eigenvectors[,1]-tmp_v1,"2")
  
  
  # #-----------------------------------------------------------------------
  # #cross-validation
  # #88888888888888888888888888888888888888888888888888888888888888888888
  # #SPCA_SOFT_cv
  # pca_soft_cv1 = sPCArSVDcv(X = X, varnum_grid = varnum_grid, type = 1, CV = CV)
  # R1 = X - tcrossprod(pca_soft_cv1$u,pca_soft_cv1$v) #residual matrix
  # pca_soft_cv2 = sPCArSVDcv(X = R1, varnum_grid = varnum_grid, type = 1, CV = CV)
  # Load2_soft_cv[,1] = pca_soft_cv1$v/norm(pca_soft_cv1$v, "2") #extract loading vector 1
  # Load2_soft_cv[,2] = pca_soft_cv2$v/norm(pca_soft_cv2$v, "2") #extract loading vector 2
  # tmp_v1 = Load2_soft_cv[,1]
  # tmp_v2 = Load2_soft_cv[,2]
  # #angle for v1
  # err2_soft_cv[nn, 1] = acos( (crossprod(eigenvectors[,1], tmp_v1)) /
  #                               (norm(eigenvectors[,1],"2")*norm(tmp_v1,"2")) ) * 180 / pi
  # #true negative rate for v1 (%)
  # err2_soft_cv[nn, 2] = CalCorrect(eigenvectors[,1], tmp_v1)[1]
  # #false negative rate for v1 (%)
  # err2_soft_cv[nn, 3] = CalCorrect(eigenvectors[,1], tmp_v1)[2]
  # #angle for v2
  # err2_soft_cv[nn, 4] = acos( (crossprod(eigenvectors[,2], tmp_v2)) /
  #                               (norm(eigenvectors[,2],"2")*norm(tmp_v2,"2")) ) * 180 / pi
  # #true negative rate for v2 (%)
  # err2_soft_cv[nn, 5] = CalCorrect(eigenvectors[,2], tmp_v2)[1]
  # #false negative rate for v2 (%)
  # err2_soft_cv[nn, 6] = CalCorrect(eigenvectors[,2], tmp_v2)[2]
  # #l2-norm between v1 and hat(v1)
  # err2_soft_cv[nn, 7] = norm(eigenvectors[,1]-tmp_v1,"2")
  # #l2-norm between v2 and hat(v2)
  # err2_soft_cv[nn, 8] = norm(eigenvectors[,2]-tmp_v2,"2")
  # 
  # #selected varnum
  # sel_varnum_soft[nn,1] = pca_soft_cv1$varnum_min
  # sel_varnum_soft[nn,2] = pca_soft_cv2$varnum_min
  # 
  # #99999999999999999999999999999999999999999999999999999999999999999999
  # #SPCA_HARD_cv
  # pca_hard_cv1 = sPCArSVDcv(X = X, varnum_grid = varnum_grid, type = 2, CV = CV)
  # R1 = X - tcrossprod(pca_hard_cv1$u,pca_hard_cv1$v) #residual matrix
  # pca_hard_cv2 = sPCArSVDcv(X = R1, varnum_grid = varnum_grid, type = 2, CV = CV)
  # Load2_hard_cv[,1] = pca_hard_cv1$v/norm(pca_hard_cv1$v, "2") #extract loading vector 1
  # Load2_hard_cv[,2] = pca_hard_cv2$v/norm(pca_hard_cv2$v, "2") #extract loading vector 2
  # tmp_v1 = Load2_hard_cv[,1]
  # tmp_v2 = Load2_hard_cv[,2]
  # #angle for v1
  # err2_hard_cv[nn, 1] = acos( (crossprod(eigenvectors[,1], tmp_v1)) /
  #                               (norm(eigenvectors[,1],"2")*norm(tmp_v1,"2")) ) * 180 / pi
  # #true negative rate for v1 (%)
  # err2_hard_cv[nn, 2] = CalCorrect(eigenvectors[,1], tmp_v1)[1]
  # #false negative rate for v1 (%)
  # err2_hard_cv[nn, 3] = CalCorrect(eigenvectors[,1], tmp_v1)[2]
  # #angle for v2
  # err2_hard_cv[nn, 4] = acos( (crossprod(eigenvectors[,2], tmp_v2)) /
  #                               (norm(eigenvectors[,2],"2")*norm(tmp_v2,"2")) ) * 180 / pi
  # #true negative rate for v2 (%)
  # err2_hard_cv[nn, 5] = CalCorrect(eigenvectors[,2], tmp_v2)[1]
  # #false negative rate for v2 (%)
  # err2_hard_cv[nn, 6] = CalCorrect(eigenvectors[,2], tmp_v2)[2]
  # #l2-norm between v1 and hat(v1)
  # err2_hard_cv[nn, 7] = norm(eigenvectors[,1]-tmp_v1,"2")
  # #l2-norm between v2 and hat(v2)
  # err2_hard_cv[nn, 8] = norm(eigenvectors[,2]-tmp_v2,"2")
  # 
  # #selected varnum
  # sel_varnum_hard[nn,1] = pca_hard_cv1$varnum_min
  # sel_varnum_hard[nn,2] = pca_hard_cv2$varnum_min
  # 
  # #101010101010101010101010101010101010101010101010101010101010101010101010
  # #SPCA_SCAD_cv
  # pca_SCAD_cv1 = sPCArSVDcv(X = X, varnum_grid = varnum_grid, type = 3, CV = CV)
  # R1 = X - tcrossprod(pca_SCAD_cv1$u,pca_SCAD_cv1$v) #residual matrix
  # pca_SCAD_cv2 = sPCArSVDcv(X = R1, varnum_grid = varnum_grid, type = 3, CV = CV)
  # Load2_SCAD_cv[,1] = pca_SCAD_cv1$v/norm(pca_SCAD_cv1$v, "2") #extract loading vector 1
  # Load2_SCAD_cv[,2] = pca_SCAD_cv2$v/norm(pca_SCAD_cv2$v, "2") #extract loading vector 2
  # tmp_v1 = Load2_SCAD_cv[,1]
  # tmp_v2 = Load2_SCAD_cv[,2]
  # #angle for v1
  # err2_SCAD_cv[nn, 1] = acos( (crossprod(eigenvectors[,1], tmp_v1)) /
  #                               (norm(eigenvectors[,1],"2")*norm(tmp_v1,"2")) ) * 180 / pi
  # #true negative rate for v1 (%)
  # err2_SCAD_cv[nn, 2] = CalCorrect(eigenvectors[,1], tmp_v1)[1]
  # #false negative rate for v1 (%)
  # err2_SCAD_cv[nn, 3] = CalCorrect(eigenvectors[,1], tmp_v1)[2]
  # #angle for v2
  # err2_SCAD_cv[nn, 4] = acos( (crossprod(eigenvectors[,2], tmp_v2)) /
  #                               (norm(eigenvectors[,2],"2")*norm(tmp_v2,"2")) ) * 180 / pi
  # #true negative rate for v2 (%)
  # err2_SCAD_cv[nn, 5] = CalCorrect(eigenvectors[,2], tmp_v2)[1]
  # #false negative rate for v2 (%)
  # err2_SCAD_cv[nn, 6] = CalCorrect(eigenvectors[,2], tmp_v2)[2]
  # #l2-norm between v1 and hat(v1)
  # err2_SCAD_cv[nn, 7] = norm(eigenvectors[,1]-tmp_v1,"2")
  # #l2-norm between v2 and hat(v2)
  # err2_SCAD_cv[nn, 8] = norm(eigenvectors[,2]-tmp_v2,"2")
  # 
  # #selected varnum
  # sel_varnum_SCAD[nn,1] = pca_SCAD_cv1$varnum_min
  # sel_varnum_SCAD[nn,2] = pca_SCAD_cv2$varnum_min
}#end for



#merge all results in a table
output = array(0, c(7, 8))
output[1,] = c(median(err2_pca[,1]), colMeans(err2_pca[,2:3]) * 100, mean(err2_pca[,1]),
               median(err2_pca[,4]), colMeans(err2_pca[,5:6]) * 100, mean(err2_pca[,4]))
output[2,] = c(median(err2_soft[,1]), colMeans(err2_soft[,2:3]) * 100, mean(err2_soft[,1]),
               median(err2_soft[,4]), colMeans(err2_soft[,5:6]) * 100, mean(err2_soft[,4]))
output[3,] = c(median(err2_hard[,1]), colMeans(err2_hard[,2:3]) * 100, mean(err2_hard[,1]),
               median(err2_hard[,4]), colMeans(err2_hard[,5:6]) * 100, mean(err2_hard[,4]))
output[4,] = c(median(err2_SCAD[,1]), colMeans(err2_SCAD[,2:3]) * 100, mean(err2_SCAD[,1]),
               median(err2_SCAD[,4]), colMeans(err2_SCAD[,5:6]) * 100, mean(err2_SCAD[,4]))
output[5,] = c(median(err2_Simple[,1]), colMeans(err2_Simple[,2:3]) * 100, mean(err2_Simple[,1]),
               median(err2_Simple[,4]), colMeans(err2_Simple[,5:6]) * 100,  mean(err2_Simple[,4]))
output[6,] = c(median(err2_SPCAk2[,1]), colMeans(err2_SPCAk2[,2:3]) * 100, mean(err2_SPCAk2[,1]),
               median(err2_SPCAk2[,4]), colMeans(err2_SPCAk2[,5:6]) * 100, mean(err2_SPCAk2[,4]))
output[7,c(1:4)] = c(median(err2_SPCAk1[,1]), colMeans(err2_SPCAk1[,2:3]) * 100, mean(err2_SPCAk1[,1]))

output



Sys.setenv(JAVA_HOME='C:/Program Files/Java/jdk1.8.0_131/jre')
library(XLConnect)
library(xlsx)
wb <- XLConnect::loadWorkbook("tmp.xlsx", create = TRUE)
sheetname="Sheet1"
XLConnect::writeWorksheet(wb,round(output,digits=2),sheetname,startRow = 1, startCol = 1, header = FALSE)
XLConnect::saveWorkbook(wb)



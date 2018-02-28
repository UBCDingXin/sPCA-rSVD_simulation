//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp:plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
using namespace Rcpp;
using namespace arma;

arma::vec signVec(arma::vec x){
  arma::vec output(x.n_elem,fill::zeros);
  for (unsigned int i = 0; i < x.n_elem; i++){
    if (x(i) > 0){
      output(i) = 1;
    }else if(x(i) < 0){
      output(i) = -1;
    }else {
      output(i) = 0;
    }
  }
  return output;
}//end signVec


//function to implement the soft, hard, SCAD thresholding rule
arma::vec thresh(arma::vec z, int type, int varnum, double a = 3.7){
  /*
   z: argument
   type: thresholding rule type
   delta: thresholding level
   a: default choice for SCAD penalty
   */
  int n = z.n_elem;
  arma::vec output, tmp;
  double delta;
  
  if (varnum < n){
    tmp = sort(abs(z),"ascend");
    delta = tmp(n - varnum - 1);
  }else {
    delta = min(abs(z)) / 2;
  }//
  
  // soft 
  if (type == 1){
    arma::vec tmp(n, fill::zeros);
    for (int i = 0; i < n; i++){
      if (abs(z(i)) >= delta){tmp(i) = 1;}
    }
    output = signVec(z) % tmp % ( abs(z) - delta );
  }//end soft
  
  //hard
  if (type == 2){
    arma::vec tmp(n, fill::zeros);
    for (int i = 0; i < n; i++){
      if (abs(z(i)) > delta){tmp(i) = 1;}
    }
    output = z % tmp;
  }//end hard
  
  //SCAD
  if (type == 3){
    arma::vec tmp1(n, fill::zeros), tmp2(n, fill::zeros), 
          tmp3(n, fill::zeros), tmp4(n, fill::zeros), tmp5(n, fill::zeros);
    for (int i = 0; i < n; i ++){
      if (abs(z(i)) >= delta){tmp1(i) = 1;}
      if (abs(z(i)) <= 2 * delta){tmp2(i) = 1;}
      if (2 * delta < abs(z(i))){tmp3(i) = 1;}
      if (abs(z(i)) <= a * delta){tmp4(i) = 1;}
      if (abs(z(i)) > a * delta){tmp5(i) = 1;}
    }
    output = signVec(z) % tmp1 % (abs(z)-delta) % tmp2 
      + ((a-1)*z-signVec(z)*a*delta) / (a-2) % tmp3 % tmp4 + z % tmp5;
  }//end SCAD
  
  return output;
}//end for thresh


//different penalty function
double penalty(arma::vec z, int type, int varnum, double a = 3.7){
  /*
  z: argument
  type: thresholding rule type
  delta: thresholding level
  a: default choice for SCAD penalty
  */
  int n = z.n_elem;
  arma::vec tmp;
  double delta, output = 0;
  
  if (varnum < n){
    tmp = sort(abs(z),"ascend");
    delta = tmp(n - varnum - 1);
  }else {
    delta = min(abs(z)) / 2;
  }//
  
  // soft 
  if (type == 1){
    output = 2 * delta * sum(abs(z));
  }//end soft
  
  //hard
  if (type == 2){
    arma::vec tmp(n, fill::zeros);
    for (int i = 0; i < n; i++){
      if (abs(z(i)) != 0){tmp(i) = 1;}
    }
    output = sum( pow(delta, 2) * tmp );
  }//end hard
  
  //SCAD
  if (type == 3){
    arma::vec tmp1(n, fill::zeros), tmp2(n, fill::zeros), 
    tmp3(n, fill::zeros);
    for (int i = 0; i < n; i ++){
      if (abs(z(i)) <= delta){tmp1(i) = 1;}
      if ((abs(z(i)) <= a * delta)&&(abs(z(i)) > delta)){tmp2(i) = 1;}
      if (abs(z(i)) > a * delta){tmp3(i) = 1;}
     }
    output = sum( 2*delta*abs(z)*tmp1 
                    - ( pow(z,2) - 2*a*delta*abs(z) + pow(delta,2) ) / (a-1) * tmp2 
                    + (a + 1) * pow(delta,2) * tmp3 );
  }//end SCAD
  
  return output;
}//end penalty function



//---------------------------------------------------------------------------------------
// function for sparse PCA via regularized SVD
// [[Rcpp::export]]
List sPCArSVD(arma::mat X, arma::vec u, arma::vec v,
              int varnum, int type = 1,
              int niter = 5000, double err = 1e-5, bool trace = false){
  /*
   X: data (predictor) matrix
   u, v: initial values. For example, u=svd(X)$u, v=svd(X)$v*svd(X)$d
   varnum: tuning parameter, # of nonzero loadings (i.e. complement of sparsity);
   type: penalty type; 1 for soft, 2 for hard, 3 for SCAD with a=3.7
   niter: # of maximum iterations
   err: threshold to declare converge
  */
  
  double u_dis = 1, v_dis = 1; //store distance between old u and new u
  int iter = 0; //iteration number
  vec v1, u1;
  double length_u1;//length of u1
  
  while ((u_dis > err) || (v_dis > err)){
    iter ++;
    
    //update v and u
    v1 = X.t() * u;
    v1 = thresh(v1, type, varnum, 3.7);//update v1
    u1 = X * v1;
    length_u1 = as_scalar(sqrt(sum(pow(u1, 2))));
    if(length_u1 > 0){
      u1 = u1 / length_u1;
    }
    // u1 = u1 / length_u1;
    u_dis = sqrt(sum(pow(u1 - u, 2)));
    v_dis = sqrt(sum(pow(v1 - v, 2)));
    
    if (iter > niter){
      std::cout<<"Type"<<type<<";Fail to converge!"<<std::endl;
      break;
    }// print iter
    u = u1;
    v = v1;
  }//end while
  
  if (trace){
    std::cout<<"iter:"<<iter<<"u_dis"<<u_dis<<"v_dis"<<v_dis<<std::endl;
  }
  
  return List::create(_["u"] = u1, _["v"] = v1);
}//end main function


//---------------------------------------------------------------------------------------
//Cross-validation of sPCArSVD
// [[Rcpp::export]]
List sPCArSVDcv(arma::mat X, arma::vec varnum_grid, 
                int type = 1, int CV = 5,
                int niter = 1000, double err = 1e-3, bool trace = false){
  /*
   X: data (predictor) matrix n * p
   varnum_grid: a grid of # non-zeros entries in the loading vector
   type: penalty type; 1 for soft, 2 for hard, 3 for SCAD with a=3.7
   CV: number of folders in cross-validation
   niter: # of maximum iterations
   err: threshold to declare converge
   */
  int n = X.n_rows;
  int p = X.n_cols;
  int nvarnum = varnum_grid.n_elem;
  // vec varnum_grid(p + 1); //grid of # variables (p + 1 sparsity); from p to 0
  // for (int i = 0; i < p + 1; i++){
  //   varnum_grid(i) = p - i;
  // }//sparsity from 0 to p (p+1 parameters)
  
  int ntrainCVm = floor(static_cast<double>(n) / CV); // number of samples in each fold for training 
  int ntrainCVlast = n - ntrainCVm * (CV - 1); // if not even one fold may have more or less samples.
  vec nkp(CV);
  for (int i = 0; i < CV; i++){
    if (i != CV - 1){
      nkp(i) = static_cast<double>(1) / (ntrainCVm * p);
    }else {
      nkp(i) = static_cast<double>(1) / (ntrainCVlast * p);;
    }
  }
  
  NumericVector errs_gross; //vector to store CV scores for each sparsity
  mat Xtrain, Xtrain1, Xtrain2, Xtest;
  int varnum;
  List CVfit;
  double errtmp;
  NumericVector vtmp;
  vec uk, v_k;
  mat U, V;
  vec s;
  
  //cross-validation
  for (int j = 0; j < nvarnum; j++){ 
    varnum = varnum_grid[j];
    NumericVector CVs; // vector to store CV scores
    for (int k = 0; k < CV; k++){
      // std::cout<<k<<";"<<j<<std::endl;
      //determine the training set and testing set in CV
      if ((k > 0) && (k < CV - 1)){
        Xtrain1 = X.rows(0, ntrainCVm * k - 1);
        Xtrain2 = X.rows(ntrainCVm * (k + 1), n - 1);
        Xtrain = join_cols(Xtrain1, Xtrain2);
        Xtest = X.rows(ntrainCVm * k, ntrainCVm * (k + 1) - 1);
      }else if (k == 0){
        Xtrain = X.rows(ntrainCVm, n - 1);
        Xtest = X.rows(0, ntrainCVm - 1);
      }else{
        Xtrain = X.rows(0, n - ntrainCVlast - 1);
        Xtest = X.rows(n - ntrainCVlast, n - 1);
      }
      svd( U, s, V, Xtrain );
      // std::cout<<V.col(0)<<std::endl;
      CVfit = sPCArSVD(Xtrain, U.col(0), V.col(0) * s(0), varnum, type); // fit sPCA on traning set
      vtmp = CVfit["v"]; //extract v_wave
      v_k = vtmp;
      uk = Xtest * v_k; // project Xtest
      if (as_scalar(sqrt(sum(pow(uk, 2)))) > 0){
        uk = uk / as_scalar(sqrt(sum(pow(uk, 2))));
      }
      // uk = uk / as_scalar(sqrt(sum(pow(uk, 2))));
      errtmp = as_scalar( sqrt(sum(sum(pow(Xtest - uk * v_k.t(), 2))) ) * nkp(k) ) ;
      CVs.push_back(errtmp);
    }//for loop over CV folds
    errs_gross.push_back(sum(CVs));
  }//for loop over each sparsity
  
  int indx_min = which_min(errs_gross); // find the index of minimum
  int varnum_min = varnum_grid(indx_min); // value of varnum that leads to the minimum CV errors
  // fit sPCArSVD based on varnum_min
  svd( U, s, V, X );
  List final_fit = sPCArSVD(X, U.col(0), V.col(0) * s(0), varnum_min, type);
  
  return List::create(_["CV_scores"] = errs_gross, _["varnum_min"] = varnum_min,
                      _["u"] = final_fit["u"], _["v"] = final_fit["v"]);
}//end cv version


//---------------------------------------------------------------------------------------
//calculate the correctly identified zeros and incrrectly identified zeros
// [[Rcpp::export]]
NumericVector CalCorrect(NumericVector ix, NumericVector iy, double tol = 1e-15){
  // x, y: true vector and estimated vector
  // tol: the threshold of determining numerical zero
  int n = ix.size(); //the length of x and y must be the same 
  /*
   * true negative rate = # correctly identified negatives / # negatives in true eigenvector
   * false negative rate = # incorrectly identified negatives / # positives in true eigenvector
   */
  NumericVector x, y;
  x = abs(ix);
  y = abs(iy);
  double tneg = 0, fneg = 0; //true negative rate and false negative rate
  NumericVector output(2);
  for (int i = 0; i < n; i++){
    if ( x[i] < tol && y[i] < tol ){ //true negative
    // if (abs(x[i]) < tol && abs(y[i]) < tol ){ //true negative
      tneg =  tneg + 1.0;
    }
    if ( x[i] > tol && y[i] < tol ){ //false negative
    // if (abs(x[i]) > tol && abs(y[i]) < tol ){ //false negative
      fneg = fneg + 1.0;
    }
  }//end for loop
  // std::cout<<tneg<<";"<<sum(abs(x) < tol)<<std::endl;
  tneg = static_cast<double>(tneg) / sum(abs(x) < tol);
  fneg = static_cast<double>(fneg) / sum(abs(x) > tol);
  output[0] = tneg;
  output[1] = fneg;
  return output;
  // return tneg;
}//end CalCorrect




//---------------------------------------------------------------------------------------
//function to generate a covariance matrix according to Gram-Schmidt orthogonalization
// [[Rcpp::export]]
List gramCpp(arma::mat v_ini, arma::vec c_seq){
  unsigned int n = v_ini.n_rows;
  int m = v_ini.n_cols;
  mat v_mat(n,n,fill::zeros);
  bool flag = false;
  mat Q, R, s_mat;
  
  if (m == 1){
    v_mat.col(0) = v_ini;
  }else{
    v_mat.cols(0,(m-1)) = v_ini;
  }

  while(!flag){
    v_mat.cols(m,n-1) = randu(n,n-m);
    flag = (rank(v_mat) == n);
    // qr(Q,R,v_mat);
    qr_econ(Q,R,v_mat);
  }
  
  s_mat = Q * diagmat(c_seq) * Q.t();
  
  return List::create(_["s.mat"] = s_mat, _["v.mat"] = Q);
}//ends gram







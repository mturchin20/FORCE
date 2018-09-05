// load Rcpp
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
arma::mat GetLinearKernel(arma::mat X){
    double p = X.n_rows;
    return X.t()*X/p;
}

// [[Rcpp::export]]
arma::mat ComputePCs(arma::mat X,int top = 10){
    mat U;
    vec s;
    mat V;
    svd(U,s,V,X);
    
    mat PCs = U*diagmat(s);
    return PCs.cols(0,top-1);
}

// [[Rcpp::export]]
arma::mat droprows(arma::mat X, uvec j){
    vec ind(X.n_rows); ind.ones();
    ind.elem(j-1) = zeros(j.n_elem);
    return X.rows(arma::find(ind == 1));
}

////////////////////////////////////////////////////////////////////////////

//Below is a function for InterPath looking for additive effects for pathways

////////////////////////////////////////////////////////////////////////////


// List InterPath_Add(mat X,vec y,List regions,int cores = 1){
// [[Rcpp::export]]
List InterPath(mat X,vec y,mat GSM,mat G,mat Z,List regions,int n,int nsnp,int cores = 1){ 
    int i;
//    const int n = X.n_cols;
//    const int nsnp = X.n_rows;
    const int p = regions.size();
//    const int p = 1;
    const int q = Z.n_rows;

    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector pve(p);
    mat Lambda(n,p);
//    mat K2;   

    //Pre-compute the Linear GSM
//    mat GSM = GetLinearKernel(X);
    
//    omp_set_num_threads(cores);
    for(i=0; i<p; i++){
        //Pre-compute the Linear GSM
        uvec j = regions[i];
        
        //Compute K covariance matrices
        mat K = (GSM*nsnp-G*j.n_elem)/(nsnp-j.n_elem-1);
        
        //Transform K using projection M
//        mat b = zeros(n,j.n_elem+1);
//        b.col(0) = ones<vec>(n); b.cols(1,j.n_elem) = trans(X.rows(j-1));
        mat b = zeros(n,q+1);
        b.col(0) = ones<vec>(n); b.cols(1,q) = Z.t();
	mat Kc = K-b*btb_inv*(b.t()*K)-(K*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(K*b))*btb_inv*b.t();
        vec yc = (eye<mat>(n,n)-(b*btb_inv)*b.t())*y;
        
        //Compute the quantities q and S
        vec q = zeros(2); //Create k-vector q to save
        mat S = zeros(2,2); //Create kxk-matrix S to save
        
        q(0) = as_scalar(yc.t()*Kc*yc);
        q(1) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);
        
        S(0,0) = as_scalar(accu(Kc.t()%Kc));
        S(0,1) = as_scalar(accu(Kc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        S(1,0) = S(0,1);
        S(1,1) = as_scalar(accu(trans(eye<mat>(n,n)-(b*btb_inv)*b.t())%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        //Compute delta and Sinv
        mat Sinv = inv(S);
        vec delta = Sinv*q;
        
        //Save point estimates and SE of the epistasis component
        sigma_est(i) = delta(0);
        
        //Find the eigenvalues of the projection matrix
        vec evals;
        eig_sym(evals,(Sinv(0,0)*Kc+Sinv(0,1)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*q(1)/S(1,1));
        Lambda.col(i) = evals;
        
        //Compute the PVE
        pve(i) = delta(0)/accu(delta);
    }
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("Eigenvalues") = Lambda,Rcpp::Named("PVE") = pve);
}

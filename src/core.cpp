#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
// #include <omp.h>
//#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]

void init_lda_v2(arma::umat& wp, arma::umat& dp, arma::uvec& ztot,
                 arma::umat weight, arma::mat& prior, arma::uvec d,
                 arma::uvec w, arma::uvec& z, int N, int T, int knowndiseases){
    int i = 0, t, ttt;
    arma::vec probs(T);
    arma::ivec count(T);
    double ss=0.0, currprob, U;
    int maxWeight = weight.max();
    
    for(i=0;i<N;i++){
        probs = prior.col(d[i]);
        ss = accu(probs);
        
        if(ss>1){
            probs = probs/ss;
        }else{
            for(t=knowndiseases;t<T;t++){
                probs[t] = (1-ss)/(T-knowndiseases);
            }
        }
        
        U = arma::randu();
        currprob = probs[0];
        ttt = 0;
        while(U>currprob){
            ttt++;
            currprob +=  probs[ttt];
        }
        z[i] = ttt;

        wp(ttt, w[i]) += weight(ttt,w[i]);
        dp(ttt, d[i]) += weight(ttt,w[i]);
        ztot[ttt] += weight(ttt,w[i]);
    }
    
    for(i=0;i<prior.n_cols;i++){
        prior.col(i) = prior.col(i)*maxWeight;
        for(t=knowndiseases;t<T;t++){
            prior(t,i) = maxWeight*0.5;
        }
    } 
}


// [[Rcpp::export]]
arma::umat lda_rcpp(arma::uvec d, arma::uvec w, arma::uvec z,
                    arma::umat weight, arma::mat prior, double alpha,
                    double beta, int T, int knowndiseases,
                    int burnin, int ITER){

    int N = d.n_elem;
    int W = weight.n_cols;
    int D = prior.n_cols;
    arma::uword word, doc, tt, ttt;
    int i, iter, t ;
    arma::vec probs(T);
    arma::uvec ztot(T);
    ztot.zeros();

    arma::umat res(T,D+W);
    res.zeros();

    double Z=0.0, U, currprob;
    double Wbeta = W*beta;

    arma::umat wp(T,W);
    wp.zeros();
    arma::umat dp(T,D);
    dp.zeros();
    
    
    Rcout << "Initiation... \n";
    
    init_lda_v2(wp, dp, ztot, weight, prior, d, w, z, N, T, knowndiseases);
    
    
    Rcout << "wp[2,3]  = "<< wp(2,3) << "\n";
    Rcout << "dp[2,3]  = "<< dp(2,3) << "\n";
    Rcout << "ztot[2]  = "<< ztot[2] << "\n";
    Rcout << "weight[2,3]  =  "<< weight(2,3) << "\n";
    Rcout << "prior[2,3]  = "<< prior(2,3) << "\n";
    
    Rcout << "Initiation done... \n";
    
    for(iter=0;iter<ITER+burnin;iter++){
        
        for(i=0;i<N;i++){
            word = w[i];
            doc = d[i];
            
            tt = z[i];
            
            wp(tt,word) -= weight(tt,word);
            dp(tt,doc) -= weight(tt,word);
            ztot[tt] -= weight(tt,word);
            
            Z = 0.0;
            for(t=0;t<T;t++){
                probs[t] = 1.0*(dp(t,doc) + prior(t,doc))*(wp(t,word) + beta)/(ztot[t] + Wbeta);
                Z = Z+probs[t];
            }
            U = Z*arma::randu();
            currprob = probs[0];
            ttt = 0;
            while(U>currprob){
                ttt++;
                currprob +=  probs[ttt];
            }
            z[i] = ttt;
            wp(ttt,word) += weight(ttt,word);
            dp(ttt,doc) += weight(ttt,word);
            ztot[ttt] += weight(ttt,word);
        }
        
        if(iter>=burnin){
            //res.slice(iter-burnin) = dp;
            for(t=0;t<T;t++){
               for(i=0;i<D;i++){
                   res(t,i) += dp(t,i);
               }
               for(i=0;i<W;i++){
                   res(t,D+i) += wp(t,i);
               }
            }
            Rcout << "Iteration " << iter-burnin << "\n";
        }
        else{
            Rcout << "Burnin " << iter << "\n";
        }
    }
    
    return res;    
}


// [[Rcpp::export]]
arma::mat lda_pred_rcpp(arma::mat weight, arma::mat X, 
                        arma::mat prior, arma::mat phi, double err=0.1){
    unsigned int N = prior.n_rows;
    unsigned int knowndiseases =prior.n_cols;
    unsigned int D = phi.n_rows;
    arma::mat pred_mat(N,D,arma::fill::zeros);
    arma::mat post_i =  phi;
    
    for(unsigned int i=0; i<N; i++){
        arma::vec prior_i(D,arma::fill::zeros);
        prior_i.subvec(0,knowndiseases-1) = prior.row(i).t();
        if (arma::accu(prior_i) > 0){
            prior_i /= arma::accu(prior_i);
            post_i = phi.each_col()%prior_i;
            post_i.each_row() /= arma::sum(post_i);
            arma::vec z_i = (post_i%weight)*X.row(i).t();
            arma::vec old(D, arma::fill::zeros);
            while(arma::any(arma::abs(z_i-old) >=err) ){
                old = z_i;
                prior_i.fill(1.0);
                prior_i.subvec(0,knowndiseases-1) = prior.row(i).t();
                prior_i += z_i;
                prior_i /= arma::accu(prior_i);
                post_i = phi.each_col()%prior_i;
                post_i.each_row() /= arma::sum(post_i);
                z_i = (post_i%weight)*X.row(i).t();
            }
            pred_mat.row(i)= z_i.t();
        }
        else{
            pred_mat.row(i).zeros();
        }
    }
    
    return pred_mat;
}


// arma::mat lda_pred_rcpp_MP(arma::mat weight, arma::mat X,
//                            arma::mat prior, arma::mat phi,
//                            unsigned int mcores=3, double err=0.1){
// 
//     omp_set_num_threads(mcores);
//     unsigned int N = prior.n_rows;
//     unsigned int knowndiseases =prior.n_cols;
//     unsigned int D = phi.n_rows;
//     arma::mat pred_mat(D,N, arma::fill::zeros);
//     //arma::mat post_i =  phi;
// 
// #pragma omp parallel for schedule(static)
//     for(unsigned int i=0; i<N; i++){
//         arma::vec prior_i(D,arma::fill::zeros);
//         prior_i.subvec(0,knowndiseases-1) = prior.row(i).t();
//         prior_i /= arma::accu(prior_i);
//         arma::mat post_i = phi.each_col()%prior_i;
//         post_i.each_row() /= arma::sum(post_i);
//         arma::vec z_i = (post_i%weight)*X.row(i).t();
//         arma::vec old(D, arma::fill::zeros);
//         while(arma::any(arma::abs(z_i-old) >=err) ){
//             old = z_i;
//             prior_i.fill(1.0);
//             prior_i.subvec(0,knowndiseases-1) = prior.row(i).t();
//             prior_i += z_i;
//             prior_i /= arma::accu(prior_i);
//             post_i = phi.each_col()%prior_i;
//             post_i.each_row() /= arma::sum(post_i);
//             z_i = (post_i%weight)*X.row(i).t();
//         }
// 
//         pred_mat.col(i)= z_i;
//     }
// 
//     return pred_mat.t();
// }

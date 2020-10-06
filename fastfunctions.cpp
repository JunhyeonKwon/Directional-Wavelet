#include <RcppArmadillo.h>
#include <algorithm>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
arma::uvec std_setdiff(arma::uvec& x, arma::uvec& y) {
  
  std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;
  
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));
  
  return arma::conv_to<arma::uvec>::from(out);
}

// [[Rcpp::export]]
arma::uvec std_union(arma::uvec a, arma::uvec b){
  std::vector<int> first = arma::conv_to<std::vector<int>>::from(arma::sort(a));
  std::vector<int> second = arma::conv_to<std::vector<int>>::from(arma::sort(b));
  std::vector<int> output;
  
  std::set_union(first.begin(),first.end(),second.begin(),second.end(),
                 std::inserter(output,output.begin()));
  
  arma::uvec result = conv_to<arma::uvec>::from(output);
  return(result);
}

// [[Rcpp::export]]
arma::uvec SelectSplatForward_arma(arma::uvec basis_cand_idx, arma::uvec basis_init_cidx, const arma::mat & eval_mat_obs, uword basis_max_num, const arma::mat & data, bool verbose){

  basis_cand_idx = basis_cand_idx - 1;
  
  if(basis_init_cidx.is_empty()){ // not sure this will work
    
    arma::vec resid_init_tmp(basis_cand_idx.n_elem);
    resid_init_tmp.fill(0);
    
    for(uword basis_cidx = 0; basis_cidx < basis_cand_idx.n_elem; basis_cidx++){
      
      uword basis_cand_idx_tmp = basis_cand_idx(basis_cidx);
      arma::mat design_mat_init = eval_mat_obs.col(basis_cand_idx_tmp);
      arma::colvec ones(design_mat_init.n_rows, fill::ones);
      design_mat_init.insert_cols(0, ones);
      arma::mat mat_id_r(design_mat_init.n_rows, design_mat_init.n_rows, fill::eye);
      arma::mat mat_id_c(design_mat_init.n_cols, design_mat_init.n_cols, fill::eye);
      
      arma::mat hat_mat_tmp = design_mat_init * solve(design_mat_init.t() * design_mat_init, mat_id_c) * design_mat_init.t();
      arma::mat resid_mat_tmp = mat_id_r - hat_mat_tmp;
      arma::vec resid_vec_tmp = resid_mat_tmp * data.col(2);
      resid_init_tmp(basis_cidx) = arma::accu(resid_vec_tmp % resid_vec_tmp);
      
    }
    
    basis_init_cidx = resid_init_tmp.index_min();
    
  }
  
  arma::uvec basis_selected_idx_old = basis_cand_idx(basis_init_cidx);
  arma::uvec basis_selected_idx_new = basis_selected_idx_old;
  
  arma::mat design_mat_tmp = eval_mat_obs.cols(basis_selected_idx_old);
  
  arma::colvec ones(design_mat_tmp.n_rows, fill::ones);
  design_mat_tmp.insert_cols(0, ones);
  
  // AIC calculation
  arma::colvec coef = arma::solve(design_mat_tmp, data.col(2));
  arma::colvec fit_val = design_mat_tmp * coef;
  arma::colvec resid = data.col(2) - fit_val;
  double aic_old = design_mat_tmp.n_rows * log(2 * datum::pi * arma::as_scalar(resid.t() * resid) / design_mat_tmp.n_rows) + design_mat_tmp.n_rows + 2 * (design_mat_tmp.n_cols+1);
  
  uword n_var = basis_selected_idx_old.n_elem;
  
  while(n_var < basis_max_num){
    
    if(verbose){
      cout << "n.var " << n_var << " / AIC:" << aic_old << endl;
    }
    
    arma::uvec remained_idx = std_setdiff(basis_cand_idx, basis_selected_idx_old);
    arma::vec aic_vec(remained_idx.n_elem);
    
    for(int i = 0; i < remained_idx.n_elem; i++){
      arma:: uvec basis_selected_idx_tmp = resize(basis_selected_idx_old, basis_selected_idx_old.n_elem + 1, 1);
      basis_selected_idx_tmp[basis_selected_idx_tmp.n_elem - 1] = remained_idx[i];
      
      arma::mat design_mat_tmp = eval_mat_obs.cols(basis_selected_idx_tmp);
      
      arma::colvec ones(design_mat_tmp.n_rows, fill::ones);
      design_mat_tmp.insert_cols(0, ones);

      
      // AIC calculation
      arma::colvec coef = arma::solve(design_mat_tmp, data.col(2));
      arma::colvec fit_val = design_mat_tmp * coef;
      arma::colvec resid = data.col(2) - fit_val;
      double aic_tmp = design_mat_tmp.n_rows * log(2 * datum::pi * arma::as_scalar(resid.t() * resid) / design_mat_tmp.n_rows) + design_mat_tmp.n_rows + 2 * (design_mat_tmp.n_cols+1);
      
      aic_vec[i] = aic_tmp;
    }
    
    arma::uword aic_idx_min = aic_vec.index_min();
    double aic_new = aic_vec[aic_idx_min];
    
    
    arma::uvec basis_selected_idx_new = resize(basis_selected_idx_old, basis_selected_idx_old.n_elem + 1, 1);
    basis_selected_idx_new[basis_selected_idx_new.n_elem - 1] = remained_idx[aic_idx_min];
    
   if(aic_new > aic_old){
     if(verbose){
       cout << "no more variable will be added.\n\n" << endl;
     }
     break;
   }else{
     aic_old = aic_new;
     basis_selected_idx_old = basis_selected_idx_new;
     n_var += 1;
   }
    
    
  }

  return basis_selected_idx_old + 1;
}



// [[Rcpp::export]]
double fastLmAIC(const arma::vec & y, const arma::mat & X) {
  
  int n = X.n_rows, k = X.n_cols;
  
  arma::colvec coef = arma::solve(X, y);
  arma::colvec resid = y - X*coef;
  
  double aic = n * log(2 * datum::pi * arma::as_scalar(resid.t() * resid) / n) + n + 2 * (k + 1);
  
  return aic;
}





// [[Rcpp::export]]
List GetPointDeviationFromCurve_arma(arma::vec p, arma::mat point_seg){
 
  arma::vec point0 = p;
  double dist_min = datum::inf;
  int sign_min = 1;
  double lambda_min = datum::nan;
  // arma::vec lambda_vec = 
  
  arma::vec point_seg_dist(point_seg.n_rows-1);
  
  for(arma::uword i = 0; i < point_seg_dist.n_elem; i++){
    point_seg_dist[i] = (point_seg(i+1, 0) - point_seg(i, 0)) * (point_seg(i+1, 0) - point_seg(i, 0)) +
      (point_seg(i+1, 1) - point_seg(i, 1)) * (point_seg(i+1, 1) - point_seg(i, 1));
    point_seg_dist[i] = sqrt(point_seg_dist[i]);
  }
  
  arma::vec tmpvec(1);
  tmpvec(0) = 0;
  point_seg_dist = join_cols(tmpvec, point_seg_dist);
  
  arma::vec point_seg_dist_cumsum = point_seg_dist;
  for(uword i = 1; i < point_seg_dist_cumsum.n_elem - 1; i++){
    for(uword j = i + 1; j < point_seg_dist_cumsum.n_elem; j++){
      point_seg_dist_cumsum(j) += point_seg_dist(i);
    }
  }
  
  arma::vec lambda_vec = point_seg_dist_cumsum;
  double dist_first = lambda_vec(1);
  for(int i = 0; i < lambda_vec.n_elem; i++){
    lambda_vec(i) = lambda_vec(i) - dist_first;
  }
  
  for(arma::uword seg_idx = 1; seg_idx < point_seg.n_rows; seg_idx++){
    arma::vec point1 = arma::vectorise(point_seg.row(seg_idx - 1));
    arma::vec point2 = arma::vectorise(point_seg.row(seg_idx));
    
    if(all(point1 == point0)){
      dist_min = 0;
      sign_min = 0;
      lambda_min = lambda_vec(seg_idx - 1);
      break;
    }
    if(all(point2 == point0)){
      dist_min = 0;
      sign_min = 0;
      lambda_min = lambda_vec(seg_idx);
      break;
    }
    
    double r_tmp = arma::accu((point0 - point1) % (point2 - point1)) / arma::accu((point2 - point1) % (point2 - point1));
    
    double dist_tmp;
    double lambda_tmp;

    if(all(point1 == point2) || r_tmp <= 0){
      dist_tmp = arma::accu((point0 - point1) % (point0 - point1));
      lambda_tmp = lambda_vec(seg_idx - 1);
    }else if(r_tmp >= 1){
      dist_tmp = arma::accu((point0 - point2) % (point0 - point2));
      lambda_tmp = lambda_vec(seg_idx);
    }else{
      dist_tmp = arma::accu((point0 - point1) % (point0 - point1)) -
        r_tmp * r_tmp * arma::accu((point1 - point2) % (point1 - point2));
      lambda_tmp = (1 - r_tmp) * lambda_vec(seg_idx - 1) +
        r_tmp * lambda_vec(seg_idx);
    }
    
    if(dist_tmp < 0){
      dist_tmp = 0;
    }else{
      dist_tmp = sqrt(dist_tmp);
    }
    
    int sign_tmp = sign((point0(1) - point1(1)) - ((point2(1) - point1(1)) /
                        (point2(0) - point1(0))*(point0(0)-point1(0))));
    
    if(dist_tmp < dist_min){
      dist_min = dist_tmp;
      sign_min = sign_tmp;
      lambda_min = lambda_tmp;
    }
    
  }
 
  
 return List::create(Named("lambda") = lambda_min,
                     Named("dev") = dist_min * sign_min); 
}



/*** R
SelectSplatForward = function(basis.cand.idx, basis.init.cidx = NULL, eval.mat.obs, basis.max.num, data, prop.var = 0.8, BIC = FALSE, verbose = FALSE, scaling = FALSE, centering = TRUE, parallelization = FALSE){
  
  AIC.old = AIC.new = 0
  
  if(is.null(basis.init.cidx)){
    
    resid.tmp = numeric(length(basis.cand.idx))
    
    for(basis.cidx in 1:length(basis.cand.idx)){
      resid.tmp.tmp = diag(x=1, nrow=dim(eval.mat.obs)[1]) - 1/sum(eval.mat.obs[,basis.cand.idx[basis.cidx]]^2) * eval.mat.obs[,basis.cand.idx[basis.cidx]] %*% t(eval.mat.obs[,basis.cand.idx[basis.cidx]])
      resid.tmp.tmp = resid.tmp.tmp %*% data[,3]
      resid.tmp[basis.cidx] = sum(resid.tmp.tmp^2)
    }
    
    basis.init.cidx = which.min(resid.tmp)
    
  }
  
  basis.selected.idx.old = basis.selected.idx.new = basis.cand.idx[basis.init.cidx]
  
  design.matrix.tmp = eval.mat.obs[, basis.selected.idx.old]
  
  # pcr.tmp = pcr.jh(response = data[,3], predictor = as.matrix(design.matrix.tmp, ncol=1), prop.var = prop.var, scaling = scaling, BIC = BIC, centering = centering)
  # AIC.old = pcr.tmp$info.crit
  lmfit.tmp = lm(data[,3]~design.matrix.tmp)
  AIC.old = AIC(lmfit.tmp)
  
  n.var = length(basis.selected.idx.old)
  
  while(n.var < basis.max.num){
    if(verbose){
      cat("n.var:", n.var, "/ AIC:", AIC.old,"\n")
    }
    
    remained.cidx  = which(!is.element(basis.cand.idx, basis.selected.idx.old))
    
    if(parallelization){
      AIC.vec = foreach(basis.cidx = remained.cidx, .combine = 'c') %dopar% {
        basis.selected.idx.tmp = c(basis.selected.idx.old, basis.cand.idx[basis.cidx])
        design.matrix.tmp = eval.mat.obs[, basis.selected.idx.tmp]
        
        fastLmAIC(data[,3], design.matrix.tmp)
      }
    }else{
      
      AIC.vec = vector()
      
      for( basis.cidx in remained.cidx ){
        
        basis.selected.idx.tmp = c(basis.selected.idx.old, basis.cand.idx[basis.cidx])
        
        design.matrix.tmp = eval.mat.obs[, basis.selected.idx.tmp]
        
        # pcr.tmp = pcr.jh(response = data[,3], predictor = design.matrix.tmp, prop.var = prop.var, scaling = scaling, BIC = BIC)
        # AIC.vec[length(AIC.vec)+1] = pcr.tmp$info.crit
        
        # lmfit = lm(data[,3]~design.matrix.tmp)
        AIC.vec[length(AIC.vec) + 1] = fastLmAIC(data[,3], design.matrix.tmp)
      }  
    }
    
    
    AIC.new = min(AIC.vec)
    
    remained.cidx.min = remained.cidx[which.min(AIC.vec)]
    basis.selected.idx.new = c(basis.selected.idx.old, basis.cand.idx[remained.cidx.min])
    
    if(AIC.new > AIC.old){
      if(verbose){
        cat("no more variable will be added.\n\n")
      }
      break
    }else{
      AIC.old = AIC.new
      basis.selected.idx.old = basis.selected.idx.new
      if(is.infinite(AIC.old)){
        break
      }
    }
    
    n.var = length(basis.selected.idx.old)
  }
  
  return(basis.selected.idx.old)
}
*/







// // [[Rcpp::export]]
// List GetDeviationFromCurve_arma(arma::mat x, arma::mat point_seg){
//   arma::mat point_seg_diff(point_seg.n_rows - 1, point_seg.n_cols);
//   
//   arma::uvec indicator_uniq(point_seg_diff.n_rows);
//   for(arma::uword i = 0; i < point_seg_diff.n_rows; i++){
//     point_seg_diff.row(i) = point_seg.row(i + 1) - point_seg.row(i);
//     if(any(point_seg_diff.row(i) != 0)){
//       indicator_uniq(i) = 1;
//     }else{
//       indicator_uniq(i) = 0;
//     }
//   }
//   
//   arma::uvec idx_uniq = find(indicator_uniq != 0);
//   
//   arma::mat point_seg_uniq(idx_uniq.n_elem, point_seg.n_cols);
//   point_seg_uniq = point_seg.rows(idx_uniq);
// 
//   // int count_uniq = 0;
//   // for(arma::uword i = 0; i < point_seg_diff.n_rows; i++){
//   //   point_seg_diff.row(i) = point_seg.row(i + 1) - point_seg.row(i);
//   //   if(any(point_seg_diff.row(i) != 0)){
//   //     count_uniq++;
//   //   }
//   // }
//   // 
//   // arma::uvec idx_uniq(count_uniq);
//   // if(count_uniq != point_seg_diff.n_rows){
//   //   int count_uniq2 = 0;
//   //   for(arma::uword i = 0; i < point_seg_diff.n_rows; i++){
//   //     if(any(point_seg_diff.row(i) != 0)){
//   //       idx_uniq(count_uniq2) = i;
//   //       count_uniq2++;
//   //     }
//   //   }
//   // }
//   // 
//   // arma::mat point_seg_diff_uniq(idx_uniq.n_elem, point_seg_diff.n_cols);
//   // point_seg_diff_uniq = point_seg_diff.rows(idx_uniq);
//   
//   arma::vec dev_vec(x.n_rows);
//   arma::vec lambda_vec(x.n_rows);
//   
//   for(arma::uword i = 0; i < x.n_rows; i++){
//     cout << "Here 0" << endl;
//     List list_tmp = GetPointDeviationFromCurve_arma(x.row(i), point_seg_uniq);
//     cout << "Here 1" << endl;
//     dev_vec(i) = list_tmp["dev"];
//     cout << "Here 2" << endl;
//     lambda_vec(i) = list_tmp["lambda"];
//     cout << "Here 3" << endl;
//   }
//   
//   return List::create(Named("lambda") = lambda_vec,
//                       Named("dev") = dev_vec);
// }





// // [[Rcpp::export]]
// List fastLm(const arma::vec & y, const arma::mat & X) {
// 
//   int n = X.n_rows, k = X.n_cols;
// 
//   arma::colvec coef = arma::solve(X, y);
//   arma::colvec resid = y - X*coef;
// 
//   double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
//   arma::colvec stderrest =
//     arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );
// 
//   return List::create(Named("coefficients") = coef,
//                       Named("stderr")       = stderrest);
// }
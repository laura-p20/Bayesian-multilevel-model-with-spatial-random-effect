// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath> 

using namespace Rcpp;
using namespace arma;

// ====================================================================
// A. Helper Function: Multivariate Normal Sampler
// ====================================================================
// Essential for the Beta samplers, replacing mvtnorm::rmvnorm
arma::vec rmvnorm_arma(const arma::vec& mu, const arma::mat& Sigma) {
  int p = mu.n_elem;
  // Cholesky decomposition: Sigma = L * L.t()
  arma::mat L = arma::chol(Sigma, "lower");
  // Standard normal noise vector Z
  arma::vec Z = arma::randn(p);
  // Sample: mu + L * Z
  return mu + L * Z;
}

// ====================================================================
// B. phi_sampler_rcpp
// Optimizes: Residue calculation (rowsum), neighbor summation (sapply)
// ====================================================================

// [[Rcpp::export]]
arma::vec phi_sampler_rcpp(arma::vec phi, const arma::vec& kappa2_jk, 
                           const arma::vec& d_jk, double tau_phi, 
                           const arma::vec& n_jk, const Rcpp::List& W, 
                           const arma::vec& r_ijk_no_phi, // Pre-calculated r_ijk without phi component
                           const arma::ivec& mun_idx) { // Municipal index (1 to m) for each student
  
  int m = phi.n_elem;
  
  // 1. Calculate R_jk = rowsum(r_ijk_no_phi) / kappa2_jk (Efficient rowsum)
  arma::vec R_jk_sum(m, arma::fill::zeros);
  for (int i = 0; i < r_ijk_no_phi.n_elem; ++i) {
    // Indices are 1-based in R, so use [mun_idx[i] - 1]
    R_jk_sum[mun_idx[i] - 1] += r_ijk_no_phi[i];
  }
  arma::vec R_jk = R_jk_sum / kappa2_jk;
  
  // 2. Calculate Sum of Neighbors (Suma_l phi_l)
  arma::vec sum_phi_neighbors(m, arma::fill::zeros);
  for (int j = 0; j < m; ++j) {
    Rcpp::IntegerVector neighbors = W[j]; // W is R list of 1-based indices
    double sum_ph = 0.0;
    for (int l : neighbors) {
      sum_ph += phi[l - 1]; // Convert to 0-based
    }
    sum_phi_neighbors[j] = sum_ph;
  }
  
  // 3. Sample from N(mean, variance)
  arma::vec b2 = (d_jk / tau_phi) + (n_jk / kappa2_jk);
  arma::vec a = (sum_phi_neighbors / tau_phi) + R_jk;
  
  arma::vec mean = a / b2;
  arma::vec variance = 1.0 / b2;
  
  // Generate phi
  for (int j = 0; j < m; ++j) {
    phi[j] = R::rnorm(mean[j], std::sqrt(variance[j]));
  }
  
  // 4. Apply restriction (sum(phi)=0)
  phi = phi - arma::mean(phi);
  
  return phi;
}

// ====================================================================
// C. betaE_sampler_rcpp, betaM_sampler_rcpp, betaD_sampler_rcpp
// Optimizes: Matrix multiplication (t(X)%*%(X*inv_kappa)) and inversion (solve(b))
// Note: All three Beta samplers share the same highly optimized structure.
// ====================================================================

// [[Rcpp::export]]
arma::vec betaE_sampler_rcpp(const arma::mat& X_ijk, const arma::vec& r_jk_E, 
                             const arma::vec& inv_kappa, const arma::vec& mu_E, 
                             double sigma2_E, int pe) {
  
  // 1. Calculate B_lik = X^T * D^-1 * X (D^-1 is diag(inv_kappa))
  // X_weighted = X * D^-1 (efficiently calculated via column-wise multiplication)
  arma::mat X_weighted = X_ijk.each_col() % inv_kappa; 
  arma::mat B_lik = X_ijk.t() * X_weighted;
  
  // Matriz de Precisión: B = Sigma_prior^-1 + B_lik
  arma::mat B = (1.0 / sigma2_E) * arma::eye(pe, pe) + B_lik;
  
  // 2. Calculate vector a
  arma::vec a_lik = X_weighted.t() * r_jk_E;
  arma::vec a = (1.0 / sigma2_E) * mu_E + a_lik;
  
  // 3. Solve for mu: mu = B^-1 * a (using efficient linear solver)
  arma::vec mu = arma::solve(B, a);
  
  // 4. Muestreo Multivariado Normal
  arma::mat Sigma = arma::inv_sympd(B); // Covariance Matrix
  
  return rmvnorm_arma(mu, Sigma);
}

// [[Rcpp::export]]
arma::vec betaM_sampler_rcpp(const arma::mat& Z_jk, const arma::vec& r_jk_M, 
                             const arma::vec& inv_kappa, const arma::vec& mu_M, 
                             double sigma2_M, int pm) {
  
  arma::mat Z_weighted = Z_jk.each_col() % inv_kappa; 
  arma::mat B_lik = Z_jk.t() * Z_weighted;
  
  arma::mat B = (1.0 / sigma2_M) * arma::eye(pm, pm) + B_lik;
  
  arma::vec a_lik = Z_weighted.t() * r_jk_M;
  arma::vec a = (1.0 / sigma2_M) * mu_M + a_lik;
  
  arma::vec mu = arma::solve(B, a);
  arma::mat Sigma = arma::inv_sympd(B); 
  
  return rmvnorm_arma(mu, Sigma);
}

// [[Rcpp::export]]
arma::vec betaD_sampler_rcpp(const arma::mat& W_k, const arma::vec& r_jk_D, 
                             const arma::vec& inv_kappa, const arma::vec& mu_D, 
                             double sigma2_D, int pd) {
  
  arma::mat W_weighted = W_k.each_col() % inv_kappa; 
  arma::mat B_lik = W_k.t() * W_weighted;
  
  arma::mat B = (1.0 / sigma2_D) * arma::eye(pd, pd) + B_lik;
  
  arma::vec a_lik = W_weighted.t() * r_jk_D;
  arma::vec a = (1.0 / sigma2_D) * mu_D + a_lik;
  
  arma::vec mu = arma::solve(B, a);
  arma::mat Sigma = arma::inv_sympd(B); 
  
  return rmvnorm_arma(mu, Sigma);
}

// ====================================================================
// D. kappa2jk_sampler_rcpp
// Optimizes: Explicit R 'for' loop and 'rowsum' calculation
// ====================================================================

// [[Rcpp::export]]
arma::vec kappa2jk_sampler_rcpp(double nu_kappa, const arma::vec& n_jk, 
                                const arma::vec& r_jk, const arma::vec& kappa2_k_mun_map, 
                                const arma::ivec& mun_idx, int m) {
  
  arma::vec kappa2jk(m);
  
  // 1. r_jk_groupped = rowsum(r_jk*r_jk) (Efficient aggregation)
  arma::vec r_sq = r_jk % r_jk;
  arma::vec r_jk_groupped(m, arma::fill::zeros);
  for (int i = 0; i < r_sq.n_elem; ++i) {
    r_jk_groupped[mun_idx[i] - 1] += r_sq[i]; // Asume mun_idx is 1-indexado
  }
  
  // 2. Parameters for Inverse Gamma (1/Gamma)
  arma::vec a_kappa2_jk = (nu_kappa + n_jk) / 2.0;
  // kappa2_k_mun_map is the required rep(kappa2_k, m_k) vector created in R
  arma::vec b_kappa2_jk = (nu_kappa * kappa2_k_mun_map + r_jk_groupped) / 2.0;
  
  // 3. Muestreo de la inversa Gamma (1/Gamma(a, 1/b))
  for (int idx = 0; idx < m; ++idx) {
    double draw = 1.0 / R::rgamma(a_kappa2_jk[idx], 1.0 / b_kappa2_jk[idx]);
    // Safety check (pmax(1e-9, ...))
    kappa2jk[idx] = std::max(1e-9, draw);
  }
  
  return kappa2jk;
}
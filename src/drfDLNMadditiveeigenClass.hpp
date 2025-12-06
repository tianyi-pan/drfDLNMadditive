#ifndef DRFDLNMADDITIVEEIGENCLASS_HPP
#define DRFDLNMADDITIVEEIGENCLASS_HPP


#ifdef _OPENMP
  #include <omp.h>
#endif

#include <lambda_lanczos.hpp>
using lambda_lanczos::LambdaLanczos;


#include <cmath>

#include <iostream>
using namespace std;



// ************** PART 1: Define some functions **********************
void choleskyAD(Eigen::MatrixXd& L) {
  // Eigen::MatrixXd will be overwritten; lower triangle will be its Cholesky. Only the lower triangle is computed/stored
  int s = L.cols();
  for (int k = 0; k < s; k++) {
    // (a) define pivot
    L(k,k) = sqrt(L(k,k));

    // (b) adjust lead column
    for (int j = k+1; j < s; j++) L(j,k) /= L(k,k);
    for (int j = k+1; j < s; j++)
      for (int i = j; i < s; i++) L(i,j) -= L(i,k) * L(j,k);
  }
}

Eigen::MatrixXd invertL(Eigen::MatrixXd &L) {
  // inverse of a lower triangular matrix
  int n = L.cols();
  Eigen::MatrixXd M(n, n);
  M.setZero();
  for (int i = 0; i < n; i++)
  {
    M(i,i) = 1.0 / L(i,i);
    for (int j = 0; j < i; j++)
    {
      for (int k = j; k < i; k++) M(i,j) += L(i,k) * M(k,j);
      M(i,j) = -M(i,j) / L(i,i);
    }
  }
  return M;
}

// check whether there is nan in the input vector
bool hasNaN(Eigen::VectorXd vec) {
    for (int i = 0; i < vec.size(); i++) {
      if( std::isnan(vec(i))) return true; // has nan
    }
    return false; // No nan
}

// TODO: the bspline function evaluated at the points outside the boundaries are incorret!
// Bspline(l=0) = 0,0,0,0... It should be a linear function of l, not always equal to 0.
int knotindexEigen(double x,const Eigen::VectorXd t) {
  int q = t.size();
  int k=0;
  if (x < t(0)) return -1;
  while(x>=t(k)){
    k++;
    if (k >= q) break;
  }

  return k-1;
}

double weight(double x, const Eigen::VectorXd& t,int i,int k) {
  if (t(i+k-1) != t(i-1))
    return((x - t(i-1))/(t(i+k-1)-t(i-1)));
  return 0.;
}


double Bspline(double x, int j, const Eigen::VectorXd& t,int p) {
  // Evaluate the jth B-spline
  // B_p(x) of order p (degree p-1) at x
  if (p==1)
    return(x>=t(j-1) && x<t(j+1-1));

  double w1 = weight(x,t,j,p-1);
  double w2 = weight(x,t,j+1,p-1);
  double b1 = Bspline(x,j,t,p-1);
  double b2 = Bspline(x,j+1,t,p-1);

  return w1*b1 + (1.-w2)*b2;
}

double Bspline1st(double x, int j, const Eigen::VectorXd& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  double bb1 = Bspline(x,j+1,t,p-1);
  double bb2 = Bspline(x,j,t,p-1);

  double ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  double ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

double Bspline2nd(double x, int j, const Eigen::VectorXd& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  double bb1 = Bspline1st(x,j+1,t,p-1);
  double bb2 = Bspline1st(x,j,t,p-1);

  double ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  double ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

Eigen::VectorXd Bsplinevec(double x, const Eigen::VectorXd& t,int p) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  // int k = knotindexEigen(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return b;
}

Eigen::VectorXd BsplinevecCon(double x, const Eigen::VectorXd& t, int p, Eigen::MatrixXd Z) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Eigen::VectorXd Bsplinevec1st(double x, const Eigen::VectorXd& t,int p) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  // int k = knotindexEigen(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return b;
}

// [[Rcpp::export]]
Eigen::VectorXd BsplinevecCon1st(double x, const Eigen::VectorXd& t, int p, Eigen::MatrixXd Z) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Eigen::VectorXd Bsplinevec2nd(double x, const Eigen::VectorXd& t,int p) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return b;
}

// [[Rcpp::export]]
Eigen::VectorXd BsplinevecCon2nd(double x, const Eigen::VectorXd& t,int p, Eigen::MatrixXd Z) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}





// https://github.com/tminka/lightspeed/blob/master/digamma.m
double lgamma1st (double x) {
  const double pi = 3.141592653589793238462643383279;
  const double large = 9.5;
  const double d1 = -0.5772156649015328606065121;
  const double d2 = pi*pi/6.0;
  const double small = 1e-6;
  const double s3 = 1.0/12.0;
  const double s4 = 1.0/120.0;
  const double s5 = 1.0/252.0;
  const double s6 = 1.0/240.0;
  const double s7 = 1.0/132.0;
  const double s8 = 691.0/32760.0;
  const double s9 = 1.0/12.0;
  const double s10 = 3617.0/8160.0;

  // Use de Moivre's expansion if x >= large = 9.5
  // calculate lgamma1st(x+10)
  double xplus10 = x + 10.0;
  double y = 0.0;
  double r = 1.0 / xplus10;
  y += log(xplus10) - 0.5 * r;
  r = r * r;
  y = y - r * ( s3 - r * ( s4 - r * (s5 - r * (s6 - r * s7))));

  // lgamma1st(x+10) = (1/x + 1/(x+1) + ... + 1/(x+9)) + lgamma1st(x)
  y = y - 1.0/x - 1.0/(x+1.0) - 1.0/(x+2.0) - 1.0/(x+3.0) - 1.0/(x+4.0) - 1.0/(x+5.0) - 1.0/(x+6.0) - 1.0/(x+7.0) - 1.0/(x+8.0) - 1/(x+9);

  return y;
}

int containsNaN(const Eigen::MatrixXd& mat) {
    int out = 0;
    for (int i = 0; i < mat.rows(); ++i) {
        for (int j = 0; j < mat.cols(); ++j) {
            if (std::isnan(mat(i, j))) {
              out++;
              // return true;
            }
        }
    }
    return out;
}



// https://github.com/tminka/lightspeed/blob/master/trigamma.m

double lgamma2nd (double x) {
  const double pi = 3.141592653589793238462643383279;
  const double c = pi*pi/6;
  const double c1 = -2.404113806319188570799476;
  const double b2 =  1.0/6.0;
  const double b4 = -1.0/30.0;
  const double b6 =  1.0/42.0;
  const double b8 = -1.0/30.0;
  const double b10 = 5.0/66.0;

  // TO DO: % Reduce to trigamma(x+n) where ( X + N ) >= large.

  double z = 1./(x*x);
  double y = 0.5*z + (1.0 + z*(b2 + z*(b4 + z*(b6 + z*(b8 + z*b10))))) / x;

  // std::cout << "x" << CppAD::Value(x) << std::endl;
  // std::cout << "trigamma" << CppAD::Value(y) << std::endl;
  return y;
}

// **************** PART 2: g(mu) = DL term + linear term + smooth term *************************
class Model {
  // The DLNM model

private:
  // DATA
  const Eigen::VectorXd y; // Response
  const std::vector<Eigen::MatrixXd> B_inner_list;
  const std::vector<Eigen::VectorXd> knots_f_list;
  const std::vector<Eigen::VectorXd> knots_w_list;


  const Eigen::MatrixXd Xfix; // fixed effects
  const Eigen::MatrixXd Xrand; // random effects
  const std::vector<Eigen::MatrixXd> Zf_list;

  const Eigen::VectorXd Xoffset; // offset


public:

  // DATA
  int n;
  int kE;
  int kEp;
  int kl; // knots for marginal lag
  int kx; // knots for marginal x
  int L; // max Lag
  int kbetaR;
  int kbetaF;
  int mE;
  int p; // number of smooth terms in Xrand

  std::vector<Eigen::MatrixXd> SfR_list;
  std::vector<Eigen::MatrixXd> SwR_list;

  const std::vector<Eigen::MatrixXd> Sf_list; // penalty matrix for f(E)
  const std::vector<Eigen::MatrixXd> Sw_list;


  Eigen::MatrixXd Blag; // lag basis
  Eigen::VectorXd Lseq; // sequence of lags


  const Eigen::VectorXd r; // rank of each smooth

  // PARAMETERS
  Eigen::VectorXd alpha_f;
  double log_theta;
  Eigen::VectorXd log_smoothing_f;
  Eigen::VectorXd log_smoothing_w;

  Eigen::VectorXd betaF; // parameters for fixed effects
  Eigen::VectorXd betaR; // parameters for random effects
  Eigen::VectorXd logsmoothing; // log smoothing parameters for random effects

  // Components generated
  double theta;
  Eigen::VectorXd smoothing_f;
  Eigen::VectorXd smoothing_w;
  Eigen::VectorXd smoothing;
  Eigen::MatrixXd Bf_matrix;
  Eigen::VectorXd Bf;
  Eigen::VectorXd Bf_tmp;
  Eigen::VectorXd eta;
  Eigen::VectorXd eta_remaining; // remaining terms = Xfix * betaF + Xrand * betaR
  Eigen::VectorXd mu; // log(mu) = eta + eta_remaining + Xoffset
  double NegLogL; // NegativeLogLikelihood value


  // Components for derivatives

  Eigen::MatrixXd dlogmu_df_mat;
  Eigen::MatrixXd dlogmu_dbetaR_mat;
  Eigen::MatrixXd dlogmu_dbetaF_mat;
  Eigen::MatrixXd dlogmu_dindex_par_mat;
  Eigen::VectorXd dlogdensity_dmu_vec;
  Eigen::MatrixXd dmu_df_mat;
  Eigen::MatrixXd dmu_dbetaR_mat;
  Eigen::MatrixXd dmu_dbetaF_mat;
  Eigen::MatrixXd dmu_dindex_par_mat;
  Eigen::MatrixXd dindex_par_dcon_index_par_mat;
  Eigen::VectorXd gr_index_par_vec;
  Eigen::VectorXd d2logdensity_dmudmu_vec;
  // std::vector<Eigen::MatrixXd> d2mu_dfdf_list;
  // std::vector<Eigen::MatrixXd> d2mu_dbetaRdbetaR_list;
  // std::vector<Eigen::MatrixXd> d2mu_dbetaFdbetaF_list;
  std::vector<Eigen::MatrixXd> d2index_par_dcon_index_pardcon_index_par_list;
  Eigen::MatrixXd he_index_par_mat;
  Eigen::MatrixXd he_alpha_f_index_par_mat;
  double dlogdensity_dtheta_scalar;
  double d2logdensity_dthetadtheta_scalar;
  Eigen::VectorXd d2logdensity_dmudtheta_vec;



  // full gradient
  Eigen::VectorXd gr_alpha_f_vec;
  Eigen::VectorXd gr_betaR_vec;
  Eigen::VectorXd gr_betaF_vec;
  Eigen::VectorXd gr_con_index_par_vec;
  Eigen::VectorXd gr_log_smoothing_f_vec;
  Eigen::VectorXd gr_log_smoothing_w_vec;
  double gr_log_theta_scalar;
  Eigen::VectorXd gr_logsmoothing_vec;

  Eigen::VectorXd gr_s_u_vec;
  Eigen::VectorXd gr_s_par_vec;

  // full hessian
  Eigen::MatrixXd he_alpha_f_mat;
  Eigen::MatrixXd he_betaR_mat;
  Eigen::MatrixXd he_betaF_mat;
  Eigen::MatrixXd he_con_index_par_mat;
  Eigen::MatrixXd he_alpha_f_con_index_par_mat;
  Eigen::MatrixXd he_alpha_f_betaF_mat;
  Eigen::MatrixXd he_alpha_f_betaR_mat;
  Eigen::MatrixXd he_betaR_betaF_mat;
  Eigen::MatrixXd he_con_index_par_betaF_mat;
  Eigen::MatrixXd he_con_index_par_betaR_mat;
  Eigen::MatrixXd he_log_smoothing_f_mat;
  Eigen::MatrixXd he_log_smoothing_w_mat;
  Eigen::MatrixXd he_logsmoothing_mat;
  double he_log_theta_scalar;
  Eigen::MatrixXd he_alpha_f_log_smoothing_f_mat;
  Eigen::MatrixXd he_alpha_f_log_smoothing_w_mat;
  Eigen::MatrixXd he_betaR_logsmoothing_mat;
  Eigen::VectorXd he_alpha_f_log_theta_vec;
  Eigen::VectorXd he_betaR_log_theta_vec;
  Eigen::VectorXd he_betaF_log_theta_vec;
  Eigen::VectorXd he_con_index_par_log_theta_vec;

  Eigen::MatrixXd he_s_u_mat;
  Eigen::MatrixXd he_s_par_u_mat;


  // third derivatives
  Eigen::MatrixXd he_s_u_mat_inv;
  Eigen::VectorXd d3logdensity_dmudmudlog_theta_vec;
  Eigen::VectorXd d2logdensity_dmudlog_theta_vec;
  Eigen::VectorXd d3logdensity_dmudmudmu_vec;
  Eigen::VectorXd Tri_alpha_f_vec;
  Eigen::VectorXd dlogmu_df_mat_t_vec;
  Eigen::VectorXd dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec;
  Eigen::VectorXd dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec;
  Eigen::VectorXd dlogmu_dbetaR_mat_t_vec;
  Eigen::VectorXd dlogmu_dbetaF_mat_t_vec;
  Eigen::VectorXd dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec;

  Eigen::VectorXd Tri_alpha_f2_betaR_vec;
  Eigen::VectorXd Tri_alpha_f2_betaF_vec;
  double Tri_alpha_f2_log_theta_scalar;

  Eigen::VectorXd Tri_betaR_vec;
  Eigen::VectorXd Tri_betaR2_alpha_f_vec;
  Eigen::VectorXd Tri_betaR2_betaF_vec;
  double Tri_betaR2_log_theta_scalar;

  Eigen::VectorXd Tri_betaF_vec;
  Eigen::VectorXd Tri_betaF2_alpha_f_vec;
  Eigen::VectorXd Tri_betaF2_betaR_vec;
  double Tri_betaF2_log_theta_scalar;

  double Tri_alpha_f_betaR_log_theta_scalar;
  Eigen::VectorXd Tri_alpha_f_betaR_alpha_f_vec;
  Eigen::VectorXd Tri_alpha_f_betaR_betaR_vec;
  Eigen::VectorXd Tri_alpha_f_betaR_betaF_vec;

  double Tri_alpha_f_betaF_log_theta_scalar;
  Eigen::VectorXd Tri_alpha_f_betaF_alpha_f_vec;
  Eigen::VectorXd Tri_alpha_f_betaF_betaR_vec;
  Eigen::VectorXd Tri_alpha_f_betaF_betaF_vec;

  double Tri_betaR_betaF_log_theta_scalar;
  Eigen::VectorXd Tri_betaR_betaF_alpha_f_vec;
  Eigen::VectorXd Tri_betaR_betaF_betaR_vec;
  Eigen::VectorXd Tri_betaR_betaF_betaF_vec;

  // To compute AIC
  double NegLogL_l; // NegativeLogLikelihood without penalty
  // matrix for I (hessian of log likelihood without penalty)
  Eigen::MatrixXd I_alpha_f_mat;
  Eigen::MatrixXd I_betaR_mat;
  Eigen::MatrixXd I_mat;
  Eigen::MatrixXd IS_mat; // I_mat with smoothing penalty

  Eigen::MatrixXd K_alpha_f_mat;
  Eigen::MatrixXd K_betaR_mat;
  Eigen::MatrixXd K_betaF_mat;
  Eigen::VectorXd K_log_theta_vec;

  Eigen::MatrixXd Kleft;
  Eigen::MatrixXd Khat; // Khat = Kleft * Kleft.transpose()


  // NCV
  Eigen::MatrixXd I_alpha_f_i_mat;
  Eigen::MatrixXd I_betaR_i_mat;
  Eigen::MatrixXd he_betaF_i_mat;
  Eigen::MatrixXd he_alpha_f_betaF_i_mat;
  Eigen::MatrixXd he_alpha_f_betaR_i_mat;
  Eigen::MatrixXd he_betaR_betaF_i_mat;

  Eigen::VectorXd he_alpha_f_log_theta_i_vec;
  Eigen::VectorXd he_betaR_log_theta_i_vec;
  Eigen::VectorXd he_betaF_log_theta_i_vec;
  Eigen::VectorXd he_con_index_par_log_theta_i_vec;


  // results for profile likelihood
  int converge; // 0: converge. 99: not converge

  double PLg; // modelobj.PLg = g.maxCoeff() in inner()

  Eigen::VectorXd logdetSplus_vec; // log of product of positive eigenvalues of penalty matrix.
  std::vector<Eigen::MatrixXd> Sinv_list; // inverse of penalty matrix
  Eigen::VectorXd eigvalS; // eigenvalues of penalty matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigS;

  // Constructor
  Model(const Eigen::VectorXd& y_,
        const std::vector<Eigen::MatrixXd>& B_inner_list_,
        const std::vector<Eigen::VectorXd>& knots_f_list_,
        const std::vector<Eigen::VectorXd>& knots_w_list_,
        const std::vector<Eigen::MatrixXd>& Sw_list_,
        const std::vector<Eigen::MatrixXd>& SwR_list_,
        const std::vector<Eigen::MatrixXd>& Sf_list_,
        const std::vector<Eigen::MatrixXd>& SfR_list_,
        const Eigen::MatrixXd& Xrand_,
        const Eigen::MatrixXd& Xfix_,
        const std::vector<Eigen::MatrixXd>& Zf_list_,
        const Eigen::VectorXd& Xoffset_,
        const Eigen::VectorXd& r_,
        Eigen::VectorXd& alpha_f_,
        double log_theta_,
        const Eigen::VectorXd& log_smoothing_f_,
        const Eigen::VectorXd& log_smoothing_w_,
        Eigen::VectorXd& betaR_,
        Eigen::VectorXd& betaF_,
        Eigen::VectorXd& logsmoothing_):
    y(y_), B_inner_list(B_inner_list_), knots_f_list(knots_f_list_), knots_w_list(knots_w_list_), Sw_list(Sw_list_), SwR_list(SwR_list_), Sf_list(Sf_list_), SfR_list(SfR_list_), Xrand(Xrand_), Xfix(Xfix_), Zf_list(Zf_list_), Xoffset(Xoffset_), r(r_),
    alpha_f(alpha_f_), log_theta(log_theta_), log_smoothing_f(log_smoothing_f_), log_smoothing_w(log_smoothing_w_), betaR(betaR_), betaF(betaF_), logsmoothing(logsmoothing_) {

      n = y.size(); // sample size
      mE = B_inner_list.size(); // number of exposures, which is the length of

      kE = alpha_f.size();
      kEp = kE / mE; // kEp = kx * kl

      kx = knots_f_list.at(0).size() - 4 - 1; // dim for a single marginal x. one for identifiability constraint (intercept).
      kl = knots_w_list.at(0).size() - 4; // dim for a single marginal lag
      L = B_inner_list.at(0).cols()-1; // max lag. ncol = L + 1
      kbetaR = betaR.size();
      kbetaF = betaF.size();

      p = r.size();




      Bf.resize(kE);
      Bf.setZero();
      Bf_tmp.resize(kE);
      Bf_tmp.setZero();

      Bf_matrix.resize(n, kE);


      theta = exp(log_theta);
    //   smoothing_f = exp(log_smoothing_f);
      smoothing_f.resize(mE);
      for (int i = 0; i < mE; i++) smoothing_f(i) = exp(log_smoothing_f(i));
    //   smoothing_w = exp(log_smoothing_w);
      smoothing_w.resize(mE);
      for (int i = 0; i < mE; i++) smoothing_w(i) = exp(log_smoothing_w(i));

      smoothing.resize(p);
      for (int i = 0; i < p; i++) smoothing(i) = exp(logsmoothing(i));





      Lseq.resize(L+1); // Set the length of Lseq to L+1
      for (int l = 0; l < (L+1); l++) {
        Lseq(l) = (double) l;
      }

      Blag.resize(kl, L+1); // basis for lag
      for (int j = 0; j < (L+1); j++) {
        Blag.col(j) = Bsplinevec(Lseq(j), knots_w_list.at(0), 4);
      }

      eta.resize(n);
      eta_remaining.resize(n);
      mu.resize(n);
      // START Bf_matrix and eta
      // B_inner.resize(n, L+1);
      for (int i = 0; i < n; i++) {
        Bf.setZero();
        for (int m = 0; m < mE; m++) {
            Bf_tmp.setZero();
            for (int j = 0; j < (L+1); j++) {
                Eigen::VectorXd bx = BsplinevecCon(B_inner_list.at(m)(i, j), knots_f_list.at(m), 4, Zf_list.at(m));
                for (int ii = 0; ii < kx; ii++) {
                    // TODO: optimize the computation of Bf_tmp
                    Bf_tmp.segment(m*kEp + ii*kl, kl) = bx(ii) * Blag.col(j);
                }
                Bf += Bf_tmp;
            }
        }
        Bf_matrix.row(i) = Bf;
        eta(i) = Bf.dot(alpha_f);
        eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
        mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
      }

      // END Bf_matrix and eta





      // Initialize the derivative components and NegativeLogLikelihood


      gr_s_u_vec.resize(kE+kbetaR+kbetaF);
      he_s_u_mat.resize(kE+kbetaR+kbetaF, kE+kbetaR+kbetaF);
      he_s_u_mat_inv.resize(kE+kbetaR+kbetaF, kE+kbetaR+kbetaF);
      I_mat.resize(kE+kbetaR+kbetaF + 1, kE+kbetaR+kbetaF + 1);
      IS_mat.resize(kE+kbetaR+kbetaF + 1, kE+kbetaR+kbetaF + 1);

      gr_s_par_vec.resize(1+2*mE+p);
      he_s_par_u_mat.resize(1+2*mE+p, kE+kbetaR+kbetaF);
      




      // calculate logdetSplus
      logdetSplus_vec.resize(mE);
      logdetSplus_vec.setZero();
      for (int m = 0; m < mE; m++) {
        eigS.compute(smoothing_f(m)*SfR_list.at(m) + smoothing_w(m)*SwR_list.at(m), Eigen::EigenvaluesOnly);
        logdetSplus_vec(m) = 0.0;
        eigvalS = eigS.eigenvalues();
        // std::cout << "eigvalS: " << eigvalS.transpose() << std::endl;
        for (int i = 0; i < eigvalS.size(); i++) {
            if (eigvalS(i) > 1e-8) logdetSplus_vec(m) += log(eigvalS(i));
        }
        Eigen::MatrixXd Sinv_tmp = (smoothing_f(m)*SfR_list.at(m) + smoothing_w(m)*SwR_list.at(m)).ldlt().solve(Eigen::MatrixXd::Identity(SfR_list.at(m).rows(), SfR_list.at(m).cols()));
        Sinv_list.push_back(Sinv_tmp);
      }



    //   derivative_coef();
    //   derivative_he();
    //   derivative_full();

    //   NegativeLogLikelihood();

    }


  // Functions to set parameters
  void setAlphaF(const Eigen::VectorXd alpha_f_) {
    alpha_f = alpha_f_;

    for (int i = 0; i < n; i++) {
      eta(i) = Bf_matrix.row(i).dot(alpha_f);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
  }




  void setBetaF(const Eigen::VectorXd betaF_) {
    betaF = betaF_;
    for (int i = 0; i < n; i++) {
      eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
  }
  void setBetaR(const Eigen::VectorXd betaR_) {
    betaR = betaR_;
    for (int i = 0; i < n; i++) {
      eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
  }
  void setLogTheta(const double log_theta_) {
    log_theta = log_theta_;
    theta = exp(log_theta);
  }

  void setLogSmoothingF(const Eigen::VectorXd log_smoothing_f_) {
    // clean Sinv_list
    Sinv_list = std::vector<Eigen::MatrixXd>();

    log_smoothing_f = log_smoothing_f_;
    for (int i = 0; i < mE; i++) smoothing_f(i) = exp(log_smoothing_f(i));
    // calculate logdetSplus
    for (int m = 0; m < mE; m++) {
        eigS.compute(smoothing_f(m)*SfR_list.at(m) + smoothing_w(m)*SwR_list.at(m), Eigen::EigenvaluesOnly);
        logdetSplus_vec(m) = 0.0;
        eigvalS = eigS.eigenvalues();
        // std::cout << "eigvalS: " << eigvalS.transpose() << std::endl;
        for (int i = 0; i < eigvalS.size(); i++) {
            if (eigvalS(i) > 1e-8) logdetSplus_vec(m) += log(eigvalS(i));
        }
        Eigen::MatrixXd Sinv_tmp = (smoothing_f(m)*SfR_list.at(m) + smoothing_w(m)*SwR_list.at(m)).ldlt().solve(Eigen::MatrixXd::Identity(SfR_list.at(m).rows(), SfR_list.at(m).cols()));
        Sinv_list.push_back(Sinv_tmp);
    }
  }

  void setLogSmoothingW(const Eigen::VectorXd log_smoothing_w_) {
    // clean Sinv_list
    Sinv_list = std::vector<Eigen::MatrixXd>();

    log_smoothing_w = log_smoothing_w_;
    for (int i = 0; i < mE; i++) smoothing_w(i) = exp(log_smoothing_w(i));
    // calculate logdetSplus
    for (int m = 0; m < mE; m++) {
        eigS.compute(smoothing_f(m)*SfR_list.at(m) + smoothing_w(m)*SwR_list.at(m), Eigen::EigenvaluesOnly);
        logdetSplus_vec(m) = 0.0;
        eigvalS = eigS.eigenvalues();
        // std::cout << "eigvalS: " << eigvalS.transpose() << std::endl;
        for (int i = 0; i < eigvalS.size(); i++) {
            if (eigvalS(i) > 1e-8) logdetSplus_vec(m) += log(eigvalS(i));
        }
        Eigen::MatrixXd Sinv_tmp = (smoothing_f(m)*SfR_list.at(m) + smoothing_w(m)*SwR_list.at(m)).ldlt().solve(Eigen::MatrixXd::Identity(SfR_list.at(m).rows(), SfR_list.at(m).cols()));
        Sinv_list.push_back(Sinv_tmp);
    }
  }

  void setLogsmoothing(const Eigen::VectorXd logsmoothing_) { // log smoothing parameters for remaining terms
    logsmoothing = logsmoothing_;
    for (int i = 0; i < p; i++) smoothing(i) = exp(logsmoothing(i));
  }


  double get_p_i (Eigen::VectorXd alpha_f_i, Eigen::VectorXd betaR_i, 
                  Eigen::VectorXd betaF_i, double log_theta_i,
                  int i) {
      
      // Bf_matrix.row(i) = Bf;
      double eta_i = Bf_matrix.row(i).dot(alpha_f_i);
      double eta_remaining_i = Xfix.row(i).dot(betaF_i) + Xrand.row(i).dot(betaR_i);
      double mu_i = exp(eta_i + eta_remaining_i + Xoffset(i));

      double theta_i = exp(log_theta_i);
      double log_p_i = lgamma(y(i) + theta_i) - lgamma(theta_i) - lgamma(y(i) + 1) -
                                    theta_i * log(1 + mu_i/theta_i) +
                                    y(i)*( eta_i + eta_remaining_i + Xoffset(i) - log_theta_i - log(1 + mu_i/theta_i) );
      return exp(log_p_i);

  }

  double get_D_i (Eigen::VectorXd alpha_f_i, Eigen::VectorXd betaR_i, 
                  Eigen::VectorXd betaF_i, double log_theta_i,
                  int i) {
      
      // Bf_matrix.row(i) = Bf;
      double eta_i = Bf_matrix.row(i).dot(alpha_f_i);
      double eta_remaining_i = Xfix.row(i).dot(betaF_i) + Xrand.row(i).dot(betaR_i);
      double mu_i = exp(eta_i + eta_remaining_i + Xoffset(i));

      double theta_i = exp(log_theta_i);
      double log_p_i = lgamma(y(i) + theta_i) - lgamma(theta_i) - lgamma(y(i) + 1) -
                                    theta_i * log(1 + mu_i/theta_i) +
                                    y(i)*( eta_i + eta_remaining_i + Xoffset(i) - log_theta_i - log(1 + mu_i/theta_i) );
      return -1.0*log_p_i; // loss function negative log likelihood

  }

  // get private members
  Eigen::MatrixXd getXfix () {
    return Xfix;
  }
  Eigen::MatrixXd getXrand () {
    return Xrand;
  }
  Eigen::VectorXd getXoffset() {
    return Xoffset;
  }


  // Function to update derivatives.
  // RUN the function derivative_coef(), derivative_he() and derivative_full() after update parameters.
  // update derivatives related to spline coefficients alpha_f, and betaR and betaF
  void derivative_coef() {
    dlogmu_df_mat = dlogmu_df();
    dlogmu_dbetaR_mat = dlogmu_dbetaR();
    dlogmu_dbetaF_mat = dlogmu_dbetaF();
    dlogdensity_dmu_vec = dlogdensity_dmu();
    dmu_df_mat = dmu_df();

    dmu_dbetaR_mat = dmu_dbetaR();
    dmu_dbetaF_mat = dmu_dbetaF();
    d2logdensity_dmudmu_vec = d2logdensity_dmudmu();


    dlogdensity_dtheta_scalar = dlogdensity_dtheta();
    d2logdensity_dmudtheta_vec = d2logdensity_dmudtheta();

    // obtain gradient
    gr_alpha_f_vec = gr_alpha_f();
    gr_betaR_vec = gr_betaR();
    gr_betaF_vec = gr_betaF();
    // obtain hessian
    he_alpha_f_mat = he_alpha_f();
    he_betaR_mat = he_betaR();
    he_betaF_mat = he_betaF();
    he_alpha_f_betaF_mat = he_alpha_f_betaF();
    he_alpha_f_betaR_mat = he_alpha_f_betaR();
    he_betaR_betaF_mat = he_betaR_betaF();
  }

  // update full gradient and hessian of alpha_f and betaR and betaF
  void derivative_he () {
    gr_s_u_vec << gr_alpha_f_vec, gr_betaR_vec, gr_betaF_vec;


    he_s_u_mat.setZero();


    he_s_u_mat.block(0, 0, kE, kE)  = he_alpha_f_mat;


    he_s_u_mat.block(kE, kE, kbetaR, kbetaR) = he_betaR_mat;
    he_s_u_mat.block(kE+kbetaR, kE+kbetaR, kbetaF, kbetaF) = he_betaF_mat;

    he_s_u_mat.block(0,kE,kE,kbetaR) = he_alpha_f_betaR_mat;

    he_s_u_mat.block(0,kE+kbetaR,kE,kbetaF) = he_alpha_f_betaF_mat;
    he_s_u_mat.block(kE, kE+kbetaR, kbetaR, kbetaF) = he_betaR_betaF_mat;

    he_s_u_mat.block(kE,0,kbetaR,kE) = he_alpha_f_betaR_mat.transpose();
    he_s_u_mat.block(kE+kbetaR,0,kbetaF,kE) = he_alpha_f_betaF_mat.transpose();
    he_s_u_mat.block(kE+kbetaR, kE, kbetaF, kbetaR) = he_betaR_betaF_mat.transpose();

    // make it symmetric. Comment out ...
    he_s_u_mat = (he_s_u_mat + he_s_u_mat.transpose())/2.0;
  }

  // update derivatives related to overdispersion and smoothing parameters
  // Full derivative for LAML
  void derivative_full () {
    // obtain full gradient

    gr_log_smoothing_f_vec = gr_log_smoothing_f();
    gr_log_smoothing_w_vec = gr_log_smoothing_w();
    gr_log_theta_scalar = gr_log_theta();
    gr_logsmoothing_vec = gr_logsmoothing();



    // u represents spline coefficient alpha_f, and betaR and betaF
    // par represents overdispersion and smoothing parameters

    gr_s_par_vec << gr_log_theta_scalar, gr_log_smoothing_f_vec, gr_log_smoothing_w_vec, gr_logsmoothing_vec;


    // obtain full hessian
    he_alpha_f_log_smoothing_f_mat = he_alpha_f_log_smoothing_f();
    he_alpha_f_log_smoothing_w_mat = he_alpha_f_log_smoothing_w();
    he_betaR_logsmoothing_mat = he_betaR_logsmoothing();
    he_alpha_f_log_theta_vec = he_alpha_f_log_theta();
    he_betaR_log_theta_vec = he_betaR_log_theta();
    he_betaF_log_theta_vec = he_betaF_log_theta();


    he_s_par_u_mat.setZero();

    he_s_par_u_mat.row(0) << he_alpha_f_log_theta_vec.transpose(), he_betaR_log_theta_vec.transpose(), he_betaF_log_theta_vec.transpose();
    he_s_par_u_mat.block(1, 0, mE, kE) = he_alpha_f_log_smoothing_f_mat.transpose();
    he_s_par_u_mat.block(1+mE, 0, mE, kE) = he_alpha_f_log_smoothing_w_mat.transpose();
    he_s_par_u_mat.block(1+2*mE, kE, p, kbetaR) = he_betaR_logsmoothing_mat.transpose();
  }

  void derivative_third (int nthreads_derivH = 1) {
    
    he_s_u_mat_inv = he_s_u_mat.ldlt().solve(Eigen::MatrixXd::Identity(kE+kbetaR+kbetaF, kE+kbetaR+kbetaF));
    


    
    
    if(nthreads_derivH > 1) {
        #ifdef _OPENMP
            omp_set_num_threads(nthreads_derivH);
        #endif
        #ifdef _OPENMP
            // OpenMP available
            // std::cout << "OpenMP is ON with " << nthreads_derivH << " threads for third derivative computation." << std::endl;
        #else
            std::cout << "Warnings: OpenMP is OFF. Install OpenMP to enable parallel processing." << std::endl;
        #endif

        Eigen::setNbThreads(1);
        dlogmu_df_mat_t_vec.resize(n);
        // Eigen::Map<const Eigen::VectorXd> vAff(he_s_u_mat_inv.block(0, 0, kE, kE).data(), kE * kE);
        dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec.resize(n);
        // Eigen::Map<const Eigen::VectorXd> vAfR(he_s_u_mat_inv.block(0,kE,kE,kbetaR).data(), kE * kbetaR);
        dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec.resize(n);
        // Eigen::Map<const Eigen::VectorXd> vAfF(he_s_u_mat_inv.block(0,kE+kbetaR,kE,kbetaF).data(), kE * kbetaF);
        dlogmu_dbetaR_mat_t_vec.resize(n);
        // Eigen::Map<const Eigen::VectorXd> vARR(he_s_u_mat_inv.block(kE, kE, kbetaR, kbetaR).data(), kbetaR * kbetaR);
        dlogmu_dbetaF_mat_t_vec.resize(n);
        // Eigen::Map<const Eigen::VectorXd> vAFF(he_s_u_mat_inv.block(kE+kbetaR, kE+kbetaR, kbetaF, kbetaF).data(), kbetaF * kbetaF);
        dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec.resize(n);
        // Eigen::Map<const Eigen::VectorXd> vARF(he_s_u_mat_inv.block(kE, kE+kbetaR, kbetaR, kbetaF).data(), kbetaR * kbetaF);
        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            // std::cout << "Thread " << omp_get_thread_num() << " is processing row " << i << std::endl;
            Eigen::MatrixXd tmpff = dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i);
            dlogmu_df_mat_t_vec(i) = (he_s_u_mat_inv.block(0, 0, kE, kE).array() * tmpff.array()).sum();

            Eigen::MatrixXd tmpfR = dlogmu_df_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i);
            dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec(i) = (he_s_u_mat_inv.block(0,kE,kE,kbetaR).array() * tmpfR.array()).sum();
            
            Eigen::MatrixXd tmpfF = dlogmu_df_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i);
            dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec(i) = (he_s_u_mat_inv.block(0,kE+kbetaR,kE,kbetaF).array() * tmpfF.array()).sum();
            
            Eigen::MatrixXd tmpRR = dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i);
            dlogmu_dbetaR_mat_t_vec(i) = (he_s_u_mat_inv.block(kE, kE, kbetaR, kbetaR).array() * tmpRR.array()).sum();
            
            Eigen::MatrixXd tmpFF = dlogmu_dbetaF_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i);
            dlogmu_dbetaF_mat_t_vec(i) = (he_s_u_mat_inv.block(kE+kbetaR, kE+kbetaR, kbetaF, kbetaF).array() * tmpFF.array()).sum();
            
            Eigen::MatrixXd tmpRF = dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i);
            dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec(i) = (he_s_u_mat_inv.block(kE, kE+kbetaR, kbetaR, kbetaF).array() * tmpRF.array()).sum();
        }
    } else {
        // std::cout << "single" << std::endl;
        dlogmu_df_mat_t_vec.resize(n);
        for (int i = 0; i < n; i++) {
            Eigen::MatrixXd tmp = dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i);
            // Eigen::Map<const Eigen::VectorXd> vB(tmp.data(), tmp.size());
            // Eigen::Map<const Eigen::VectorXd> vA(he_s_u_mat_inv.block(0, 0, kE, kE).data(), kE * kE);
            // dlogmu_df_mat_t_vec(i) = vA.dot(vB);
            dlogmu_df_mat_t_vec(i) = (he_s_u_mat_inv.block(0, 0, kE, kE).array() * tmp.array()).sum();
        }

        dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec.resize(n);
        for (int i = 0; i < n; i++) {
            Eigen::MatrixXd tmp = dlogmu_df_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i);
            // Eigen::Map<const Eigen::VectorXd> vA(he_s_u_mat_inv.block(0,kE,kE,kbetaR).data(), kE * kbetaR);
            // Eigen::Map<const Eigen::VectorXd> vB(tmp.data(), tmp.size());
            // dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec(i) = vA.dot(vB);
            dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec(i) = (he_s_u_mat_inv.block(0,kE,kE,kbetaR).array() * tmp.array()).sum();
        }

        dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec.resize(n);
        for (int i = 0; i < n; i++) {
            Eigen::MatrixXd tmp = dlogmu_df_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i);
            // Eigen::Map<const Eigen::VectorXd> vA(he_s_u_mat_inv.block(0,kE+kbetaR,kE,kbetaF).data(), kE * kbetaF);
            // Eigen::Map<const Eigen::VectorXd> vB(tmp.data(), tmp.size());
            // dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec(i) = vA.dot(vB);
            dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec(i) = (he_s_u_mat_inv.block(0,kE+kbetaR,kE,kbetaF).array() * tmp.array()).sum();
        }

        dlogmu_dbetaR_mat_t_vec.resize(n);
        for (int i = 0; i < n; i++) {
            Eigen::MatrixXd tmp = dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i);
            // Eigen::Map<const Eigen::VectorXd> vA(he_s_u_mat_inv.block(kE, kE, kbetaR, kbetaR).data(), kbetaR * kbetaR);
            // Eigen::Map<const Eigen::VectorXd> vB(tmp.data(), tmp.size());
            // dlogmu_dbetaR_mat_t_vec(i) = vA.dot(vB);
            dlogmu_dbetaR_mat_t_vec(i) = (he_s_u_mat_inv.block(kE, kE, kbetaR, kbetaR).array() * tmp.array()).sum();
        }
        dlogmu_dbetaF_mat_t_vec.resize(n);
        for (int i = 0; i < n; i++) {
            Eigen::MatrixXd tmp = dlogmu_dbetaF_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i);
            // Eigen::Map<const Eigen::VectorXd> vA(he_s_u_mat_inv.block(kE+kbetaR, kE+kbetaR, kbetaF, kbetaF).data(), kbetaF * kbetaF);
            // Eigen::Map<const Eigen::VectorXd> vB(tmp.data(), tmp.size());
            // dlogmu_dbetaF_mat_t_vec(i) = vA.dot(vB);
            dlogmu_dbetaF_mat_t_vec(i) = (he_s_u_mat_inv.block(kE+kbetaR, kE+kbetaR, kbetaF, kbetaF).array() * tmp.array()).sum();
        }
        dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec.resize(n);
        for (int i = 0; i < n; i++) {
            Eigen::MatrixXd tmp = dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i);
            // Eigen::Map<const Eigen::VectorXd> vA(he_s_u_mat_inv.block(kE, kE+kbetaR, kbetaR, kbetaF).data(), kbetaR * kbetaF);
            // Eigen::Map<const Eigen::VectorXd> vB(tmp.data(), tmp.size());
            // dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec(i) = vA.dot(vB);
            dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec(i) = (he_s_u_mat_inv.block(kE, kE+kbetaR, kbetaR, kbetaF).array() * tmp.array()).sum();
        }
    }



    d2logdensity_dmudlog_theta_vec = d2logdensity_dmudlog_theta();
    d3logdensity_dmudmudlog_theta_vec = d3logdensity_dmudmudlog_theta();
    d3logdensity_dmudmudmu_vec = d3logdensity_dmudmudmu();

    Tri_alpha_f_vec = Tri_alpha_f(nthreads_derivH);
    Tri_alpha_f2_betaR_vec = Tri_alpha_f2_betaR(nthreads_derivH);
    Tri_alpha_f2_betaF_vec = Tri_alpha_f2_betaF(nthreads_derivH);
    Tri_alpha_f2_log_theta_scalar = Tri_alpha_f2_log_theta();

    Tri_betaR_vec = Tri_betaR(nthreads_derivH);
    Tri_betaR2_alpha_f_vec = Tri_betaR2_alpha_f(nthreads_derivH);
    Tri_betaR2_betaF_vec = Tri_betaR2_betaF(nthreads_derivH);
    Tri_betaR2_log_theta_scalar = Tri_betaR2_log_theta();

    Tri_betaF_vec = Tri_betaF(nthreads_derivH);
    Tri_betaF2_alpha_f_vec = Tri_betaF2_alpha_f(nthreads_derivH);
    Tri_betaF2_betaR_vec = Tri_betaF2_betaR(nthreads_derivH);
    Tri_betaF2_log_theta_scalar = Tri_betaF2_log_theta();


    Tri_alpha_f_betaR_log_theta_scalar = Tri_alpha_f_betaR_log_theta();
    Tri_alpha_f_betaR_alpha_f_vec = Tri_alpha_f_betaR_alpha_f(nthreads_derivH);
    Tri_alpha_f_betaR_betaR_vec = Tri_alpha_f_betaR_betaR(nthreads_derivH);
    Tri_alpha_f_betaR_betaF_vec = Tri_alpha_f_betaR_betaF(nthreads_derivH);

    Tri_alpha_f_betaF_log_theta_scalar = Tri_alpha_f_betaF_log_theta();
    Tri_alpha_f_betaF_alpha_f_vec = Tri_alpha_f_betaF_alpha_f(nthreads_derivH);
    Tri_alpha_f_betaF_betaR_vec = Tri_alpha_f_betaF_betaR(nthreads_derivH);
    Tri_alpha_f_betaF_betaF_vec = Tri_alpha_f_betaF_betaF(nthreads_derivH);

    Tri_betaR_betaF_log_theta_scalar = Tri_betaR_betaF_log_theta();
    Tri_betaR_betaF_alpha_f_vec = Tri_betaR_betaF_alpha_f(nthreads_derivH);
    Tri_betaR_betaF_betaR_vec = Tri_betaR_betaF_betaR(nthreads_derivH);
    Tri_betaR_betaF_betaF_vec = Tri_betaR_betaF_betaF(nthreads_derivH);
  }


  // functions for NegativeLogLikelihood
  void NegativeLogLikelihood() {

    double loglik = 0;
    for (int i = 0; i < n; i++) {
      loglik += lgamma(y(i) + theta) - lgamma(theta) - lgamma(y(i) + 1) -
                                    theta * log(1 + mu(i)/theta) +
                                    y(i)*( eta(i) + eta_remaining(i) + Xoffset(i) - log_theta - log(1 + mu(i)/theta) );
    }
    // part 1: DLNM
    for (int i = 0; i < mE; i++) {
        loglik += - 0.5 * smoothing_w(i) * alpha_f.segment(i*kEp, kEp).dot(Sw_list.at(i) * alpha_f.segment(i*kEp, kEp)) - 0.5 * smoothing_f(i) * alpha_f.segment(i*kEp, kEp).dot(Sf_list.at(i) * alpha_f.segment(i*kEp, kEp));
        loglik += logdetSplus_vec(i)/2.0;
    }


    // part 2: Remaining smooth terms
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      Eigen::VectorXd betaRi(ki);
      for (int j = 0; j < ki; j++) betaRi(j) = betaR(begin + j);
      loglik += -0.5 * smoothing(i) * betaRi.dot(betaRi); // smooth penalty
      loglik += ki/2.0 * logsmoothing(i); // scale

      begin += ki;
    }


    NegLogL = -1.0 * loglik; // NEGATIVE log-likelihood
  }

  // functions for NegativeLogLikelihood WITHOUT penalty for AIC
  void NegativeLogLikelihood_l() {

    double loglik = 0;
    for (int i = 0; i < n; i++) {
      loglik += lgamma(y(i) + theta) - lgamma(theta) - lgamma(y(i) + 1) -
                                    theta * log(1 + mu(i)/theta) +
                                    y(i)*( eta(i) + eta_remaining(i) + Xoffset(i) - log_theta - log(1 + mu(i)/theta) );
    }

    NegLogL_l = -1.0 * loglik; // NEGATIVE log-likelihood
  }

  void prepare_AIC () {
    NegativeLogLikelihood_l();
    // hessian of log likelihood without penalty
    I_alpha_f_mat = I_alpha_f();
    I_betaR_mat = I_betaR();

    I_mat.block(0,0,kE+kbetaR+kbetaF, kE+kbetaR+kbetaF) = he_s_u_mat;
    I_mat.block(0, 0, kE, kE)  = I_alpha_f_mat;
    I_mat.block(kE, kE, kbetaR, kbetaR) = I_betaR_mat;


    d2logdensity_dthetadtheta_scalar = d2logdensity_dthetadtheta();
    he_log_theta_scalar = he_log_theta();



    I_mat.row(kE+kbetaR+kbetaF) << he_s_par_u_mat.row(0), he_log_theta_scalar;
    I_mat.col(kE+kbetaR+kbetaF) = I_mat.row(kE+kbetaR+kbetaF).transpose();

    IS_mat.block(0,0,kE+kbetaR+kbetaF, kE+kbetaR+kbetaF) = he_s_u_mat;
    IS_mat.row(kE+kbetaR+kbetaF) << he_s_par_u_mat.row(0), he_log_theta_scalar;
    IS_mat.col(kE+kbetaR+kbetaF) = IS_mat.row(kE+kbetaR+kbetaF).transpose();


    // Matrix K
    K_alpha_f_mat = K_alpha_f();
    K_betaR_mat = K_betaR();
    K_betaF_mat = K_betaF();
    K_log_theta_vec = K_log_theta();

    Kleft.resize(kE+kbetaR+kbetaF+1, n);
    Kleft.setZero();

    Khat.resize(kE+kbetaR+kbetaF+1, kE+kbetaR+kbetaF+1);

    Kleft.block(0, 0, kE, n) = K_alpha_f_mat;
    Kleft.block(kE,0, kbetaR,n) = K_betaR_mat;
    Kleft.block(kE+kbetaR,0, kbetaF,n) = K_betaF_mat;
    Kleft.row(kE+kbetaR+kbetaF) = K_log_theta_vec.transpose();

    Khat = Kleft * Kleft.transpose();

  }
  void prepare_IS () {
    // IS matrix
    d2logdensity_dthetadtheta_scalar = d2logdensity_dthetadtheta();
    he_log_theta_scalar = he_log_theta();

    IS_mat.block(0,0,kE+kbetaR+kbetaF, kE+kbetaR+kbetaF) = he_s_u_mat;
    IS_mat.row(kE+kbetaR+kbetaF) << he_s_par_u_mat.row(0), he_log_theta_scalar;
    IS_mat.col(kE+kbetaR+kbetaF) = IS_mat.row(kE+kbetaR+kbetaF).transpose();



    // // Matrix K
    // K_alpha_f_mat = K_alpha_f();
    // K_betaR_mat = K_betaR();
    // K_betaF_mat = K_betaF();
    // K_log_theta_vec = K_log_theta();

    // Kleft.resize(kE+kbetaR+kbetaF+1, n);
    // Kleft.setZero();

    // Khat.resize(kE+kbetaR+kbetaF+1, kE+kbetaR+kbetaF+1);

    // Kleft.block(0, 0, kE, n) = K_alpha_f_mat;
    // Kleft.block(kE,0, kbetaR,n) = K_betaR_mat;
    // Kleft.block(kE+kbetaR,0, kbetaF,n) = K_betaF_mat;
    // Kleft.row(kE+kbetaR+kbetaF) = K_log_theta_vec.transpose();

    // Khat = Kleft * Kleft.transpose();
  }

  Eigen::MatrixXd prepare_NCV (Eigen::VectorXd nei_vec) {
    Eigen::MatrixXd Hunpen_nei(kE+kbetaR+kbetaF + 1, kE+kbetaR+kbetaF + 1);
    Hunpen_nei.setZero();
    Eigen::MatrixXd tmp(kE+kbetaR+kbetaF + 1, kE+kbetaR+kbetaF + 1);
    tmp.setZero();

    // compute Hunpen_nei
    for (size_t i = 0; i < nei_vec.size(); i++) {
      int index_int = static_cast<int>(nei_vec(i)) - 1; // start from 0 in C++
      I_alpha_f_i_mat = I_alpha_f_i(index_int);
      I_betaR_i_mat = I_betaR_i(index_int);
      he_betaF_i_mat = he_betaF_i(index_int);
      he_alpha_f_betaF_i_mat = he_alpha_f_betaF_i(index_int);
      he_alpha_f_betaR_i_mat = he_alpha_f_betaR_i(index_int);
      he_betaR_betaF_i_mat = he_betaR_betaF_i(index_int);

      he_alpha_f_log_theta_i_vec = he_alpha_f_log_theta_i(index_int);
      he_betaR_log_theta_i_vec = he_betaR_log_theta_i(index_int);
      he_betaF_log_theta_i_vec = he_betaF_log_theta_i(index_int);

      tmp.setZero();
      tmp.block(0, 0, kE, kE)  = I_alpha_f_i_mat;

      

      tmp.block(kE, kE, kbetaR, kbetaR) = I_betaR_i_mat;
      tmp.block(kE+kbetaR, kE+kbetaR, kbetaF, kbetaF) = he_betaF_i_mat;

      tmp.block(0,kE,kE,kbetaR) = he_alpha_f_betaR_i_mat;
      tmp.block(0,kE+kbetaR,kE,kbetaF) = he_alpha_f_betaF_i_mat;
      tmp.block(kE, kE+kbetaR, kbetaR, kbetaF) = he_betaR_betaF_i_mat;

      tmp.block(kE,0,kbetaR,kE) = he_alpha_f_betaR_i_mat.transpose();
      tmp.block(kE+kbetaR,0,kbetaF,kE) = he_alpha_f_betaF_i_mat.transpose();
      tmp.block(kE+kbetaR, kE, kbetaF, kbetaR) = he_betaR_betaF_i_mat.transpose();

      tmp.row(kE+kbetaR+kbetaF) << he_alpha_f_log_theta_i_vec.transpose(), he_betaR_log_theta_i_vec.transpose(), he_betaF_log_theta_i_vec.transpose(), he_log_theta_i(index_int);
      tmp.col(kE+kbetaR+kbetaF) = tmp.row(kE+kbetaR+kbetaF).transpose();

      Hunpen_nei += tmp;
    }

    return Hunpen_nei;

  }

  // ********* Derivatives *************

  // FUNCTIONS
  // 1. density function
  // d log(exponential family density) / d mu
  Eigen::VectorXd dlogdensity_dmu () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = y(i) / mu(i) - (theta + y(i)) / (theta + mu(i));
    }
    return out;
  }
  // d^2 log(exponential family density) / d mu d log_theta
  Eigen::VectorXd d2logdensity_dmudlog_theta () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = theta * ((theta + y(i)) / pow(theta + mu(i), 2) - 1.0 / (theta + mu(i)));
    }
    return out;
  }

  // d^2 log(exponential family density) / d mu^2
  Eigen::VectorXd d2logdensity_dmudmu () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = - y(i) / pow(mu(i), 2) + (theta + y(i)) / pow(theta + mu(i), 2);
    }
    return out;
  }
  // d^3 log(exponential family density) / d mu^2 d log_theta
  Eigen::VectorXd d3logdensity_dmudmudlog_theta () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = theta * (1 / pow(theta + mu(i), 2) - 2.0 * (theta + y(i)) / pow(theta + mu(i), 3));
    }
    return out;
  }

  Eigen::VectorXd d3logdensity_dmudmudmu () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = 2.0 * y(i) / pow(mu(i), 3) - 2.0 * (theta + y(i)) / pow(theta + mu(i), 3);
    }
    return out;
  }

  // d log(exponential family density) / d theta
  double dlogdensity_dtheta () {
    double out = 0.0;
    // std::cout << "x" << 3.5 << std::endl;

    // TO DO: optimize it. Use property of gamma function...
    for (int i = 0; i < n; i++) {
      out += log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta);
    }
    return out;
  }
  // d^2 log(exponential family density) / d theta^2
  double d2logdensity_dthetadtheta () {
    double out = 0.0;

    for (int i = 0; i < n; i++) {
      out += 1/theta - 1/(theta + mu(i)) - (mu(i) - y(i)) / ((theta + mu(i))*(theta + mu(i))) + lgamma2nd(y(i) + theta) - lgamma2nd(theta);
    }
    return out;
  }


  Eigen::VectorXd d2logdensity_dmudtheta () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = (y(i) - mu(i)) / pow(theta+mu(i), 2);
    }
    return out;
  }



  // // 2. mean model
  // d log(mu) / d alpha_f
  Eigen::MatrixXd dlogmu_df () {
    return Bf_matrix;
  }
  // d mu / d alpha_f
  Eigen::MatrixXd dmu_df () {
    Eigen::MatrixXd out(n, kE);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_df_mat.row(i) * mu(i);
    }
    return out;
  }

  // d log(mu) / d betaR
  Eigen::MatrixXd dlogmu_dbetaR () {
    return Xrand;
  }
  // d mu / d betaR
  Eigen::MatrixXd dmu_dbetaR () {
    Eigen::MatrixXd out(n, kbetaR);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dbetaR_mat.row(i) * mu(i);
    }
    return out;
  }
  // d log(mu) / d betaF
  Eigen::MatrixXd dlogmu_dbetaF () {
    return Xfix;
  }
  // d mu / d betaR
  Eigen::MatrixXd dmu_dbetaF () {
    Eigen::MatrixXd out(n, kbetaF);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dbetaF_mat.row(i) * mu(i);
    }
    return out;
  }





  // *** GRADIENT ***
  Eigen::VectorXd gr_alpha_f () {
    Eigen::VectorXd tmp(kE);
    for (int i = 0; i < mE; i++) {
        tmp.segment(i*kEp, kEp) = smoothing_w(i) * Sw_list.at(i) * alpha_f.segment(i*kEp, kEp) + smoothing_f(i) * Sf_list.at(i) * alpha_f.segment(i*kEp, kEp);
    }
    Eigen::VectorXd out = - dmu_df_mat.transpose() * dlogdensity_dmu_vec + tmp;
    return out;
  }
  Eigen::MatrixXd K_alpha_f () {
    // dmu_df_mat: n * kE
    // dlogdensity_dmu_vec: n * 1
    // out: kE * n
    Eigen::MatrixXd out(kE, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = dmu_df_mat.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }

  Eigen::VectorXd gr_betaR () {
    Eigen::VectorXd out = - dmu_dbetaR_mat.transpose() * dlogdensity_dmu_vec; // + smoothing * betaR;
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      // for (int j = 0; j < ki; j++) betaRi(j) = betaR(begin + j);
      // out += smoothing(i) * betaRi;
      for (int j = 0; j < ki; j++) out(begin + j) += smoothing(i) * betaR(begin + j);
      begin += ki;
    }
    return out;
  }

  Eigen::MatrixXd K_betaR () {
    // dmu_dbetaR_mat: n * kbetaR
    // dlogdensity_dmu_vec: n * 1
    // out: kbetaR * n
    Eigen::MatrixXd out(kbetaR, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = dmu_dbetaR_mat.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }

  Eigen::VectorXd gr_betaF () {
    Eigen::VectorXd out = - dmu_dbetaF_mat.transpose() * dlogdensity_dmu_vec;
    return out;
  }


  Eigen::MatrixXd K_betaF () {
    // dmu_dbetaF_mat: n * kbetaF
    // dlogdensity_dmu_vec: n * 1
    // out: kbetaF * n
    Eigen::MatrixXd out(kbetaF, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = dmu_dbetaF_mat.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }

  Eigen::VectorXd gr_log_smoothing_f () {
    Eigen::VectorXd out(mE);
    for (int i = 0; i < mE; i++) {
        out(i) = 0.5 * smoothing_f(i) * alpha_f.segment(i*kEp, kEp).dot(Sf_list.at(i) * alpha_f.segment(i*kEp, kEp)) - 0.5 * smoothing_f(i) * (Sinv_list.at(i) * SfR_list.at(i)).trace();
    }
    // return 0.5 * smoothing_f * alpha_f.dot(Sf * alpha_f) - 0.5 * smoothing_f * (Sinv * SfR).trace();
    return out;
  }

  Eigen::VectorXd gr_log_smoothing_w () {
    Eigen::VectorXd out(mE);
    for (int i = 0; i < mE; i++) {
        out(i) = 0.5 * smoothing_w(i) * alpha_f.segment(i*kEp, kEp).dot(Sw_list.at(i) * alpha_f.segment(i*kEp, kEp)) - 0.5 * smoothing_w(i) * (Sinv_list.at(i) * SwR_list.at(i)).trace();
    }
    // return 0.5 * smoothing_f * alpha_f.dot(Sf * alpha_f) - 0.5 * smoothing_f * (Sinv * SfR).trace();
    // return 0.5 * smoothing_w * alpha_f.dot(Sw * alpha_f) - 0.5 * smoothing_w * (Sinv * SwR).trace();
    return out;
  }

  double gr_log_theta () {
    return -1.0 * theta * dlogdensity_dtheta_scalar;
  }


  Eigen::VectorXd K_log_theta () {
    // out: n * 1. DO NOT FORGET TO TRANSPOSE WHEN USING
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = -1.0 * theta * (log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta));
    }
    return out;
  }


  Eigen::VectorXd gr_logsmoothing () {
    Eigen::VectorXd out(p);
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      Eigen::VectorXd betaRi(ki);
      for (int j = 0; j < ki; j++) betaRi(j) = betaR(begin + j);
      out(i) = 0.5 * smoothing(i) * betaRi.dot(betaRi) - 0.5*ki;
      begin += ki;
    }
    return out;
  }


  double he_log_theta () {
    return -1.0 * theta * dlogdensity_dtheta_scalar - theta * theta * d2logdensity_dthetadtheta_scalar;
  }
  double he_log_theta_i (int i) {
    double tmp1 = log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta);
    double tmp2 = 1/theta - 1/(theta + mu(i)) - (mu(i) - y(i)) / ((theta + mu(i))*(theta + mu(i))) + lgamma2nd(y(i) + theta) - lgamma2nd(theta);
    return -1.0 * theta * tmp1 - theta * theta * tmp2;
  }

  Eigen::MatrixXd he_log_smoothing_f () {
    Eigen::MatrixXd out(mE, mE);
    out.setZero();
    for (int i = 0; i < mE; i++) {
        out(i, i) = 0.5 * smoothing_f(i) * alpha_f.segment(i*kEp, kEp).dot(Sf_list.at(i) * alpha_f.segment(i*kEp, kEp)) - 0.5 * smoothing_f(i) * (Sinv_list.at(i) * SfR_list.at(i)).trace();
    }
    // return 0.5 * smoothing_f * alpha_f.dot(Sf * alpha_f) - 0.5 * smoothing_f * (Sinv * SfR).trace();
    return out;
  }
  Eigen::MatrixXd he_log_smoothing_w () {
    Eigen::MatrixXd out(mE, mE);
    out.setZero();
    for (int i = 0; i < mE; i++) {
        out(i, i) = 0.5 * smoothing_w(i) * alpha_f.segment(i*kEp, kEp).dot(Sw_list.at(i) * alpha_f.segment(i*kEp, kEp)) - 0.5 * smoothing_w(i) * (Sinv_list.at(i) * SwR_list.at(i)).trace();
    }
    // return 0.5 * smoothing_w * alpha_f.dot(Sw * alpha_f) - 0.5 * smoothing_w * (Sinv * SwR).trace();
    return out;
  }

  Eigen::MatrixXd he_logsmoothing () {
    Eigen::MatrixXd out(p, p);
    out.setZero();
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      Eigen::VectorXd betaRi(ki);
      for (int j = 0; j < ki; j++) betaRi(j) = betaR(begin + j);
      out(i, i) = 0.5 * smoothing(i) * betaRi.dot(betaRi);
      begin += ki;
    }
    return out;
  }


  // *** Hessian ***
  Eigen::MatrixXd he_alpha_f () {
    Eigen::MatrixXd out1(kE, kE);
    Eigen::MatrixXd out2(kE, kE);
    Eigen::MatrixXd outpen(kE, kE);
    out1.setZero();
    out2.setZero();
    outpen.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i));
    }

    for (int i = 0; i < mE; i++) outpen.block(i*kEp, i*kEp, kEp, kEp) = smoothing_f(i)*Sf_list.at(i) + smoothing_w(i)*Sw_list.at(i);

    return - out1 - out2 + outpen;
  }

  double Tri_alpha_f2_log_theta () {
    // cannot use const Eigen::Ref<const Eigen::MatrixXd>& mat, because of Eigen::Map later
    Eigen::MatrixXd out1(kE, kE);
    Eigen::MatrixXd out2(kE, kE);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d3logdensity_dmudmudlog_theta_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
      out2 += d2logdensity_dmudlog_theta_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i));
    }
    Eigen::MatrixXd out = - out1 - out2;
    // Eigen::Map<const Eigen::VectorXd> vA(mat.data(), mat.size());
    // Eigen::Map<const Eigen::VectorXd> vB(out.data(), out.size());
    return (out.array() * he_s_u_mat_inv.block(0, 0, kE, kE).array()).sum();
  }

  // third derivatives wrt alpha_f
  Eigen::VectorXd Tri_alpha_f (int nThreads) {
    // std::vector<Eigen::MatrixXd> out;
    Eigen::VectorXd out(kE);
    double out1;
    // Eigen::MatrixXd out2(kE, kE);
    if(nThreads > 1) {
        // std::cout << "Using multithreading with " << nThreads << " threads for third derivatives wrt alpha_f." << std::endl;
        for (int j = 0; j < kE; j++) {
            // std::cout << "Computing third derivative wrt alpha_f: " << j+1 << " out of " << kE << std::endl;
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = (d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                                    d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                                d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                                    dlogdensity_dmu_vec(i) * dmu_df_mat(i,j));
                    local += tmp * dlogmu_df_mat_t_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kE; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                // out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
                //   out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i));
                double tmp = (d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_df_mat(i,j));
                out1 += tmp * dlogmu_df_mat_t_vec(i);
            }
            out(j) = -1.0 * out1;
        }
    }
    return out;
  }


  // d^3 L / d alpha_f d alpha_f d betaR
  Eigen::VectorXd Tri_alpha_f2_betaR (int nThreads) {
    Eigen::VectorXd out(kbetaR);
    double out1;

    if(nThreads > 1) {
        for (int j = 0; j < kbetaR; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                                    d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                                d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                                    dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                    local += tmp * dlogmu_df_mat_t_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaR; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                out1 += tmp * dlogmu_df_mat_t_vec(i);
            }
            out(j) = -1.0 * out1;
        }
    }

    return out;
  }
  // d^3 L / d alpha_f d alpha_f d betaF
  Eigen::VectorXd Tri_alpha_f2_betaF (int nThreads) {
    Eigen::VectorXd out(kbetaF);
    double out1;
    if(nThreads > 1) {
        for (int j = 0; j < kbetaF; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                                    d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                                d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                                    dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
                    local += tmp * dlogmu_df_mat_t_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back( - out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaF; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
                out1 += tmp * dlogmu_df_mat_t_vec(i);
            }
            // out.push_back( - out1);
            out(j) = -1.0 * out1;
        }
    }
    return out;
  }

  Eigen::MatrixXd I_alpha_f () { // hessian of negative likelihood without penalty
    Eigen::MatrixXd out1(kE, kE);
    Eigen::MatrixXd out2(kE, kE);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i));
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd I_alpha_f_i (int i) {
    Eigen::MatrixXd out1(kE, kE);
    Eigen::MatrixXd out2(kE, kE);

    out1 = d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i));

    return - out1 - out2;
  }

  Eigen::MatrixXd he_betaR () {
    Eigen::MatrixXd out1(kbetaR, kbetaR);
    Eigen::MatrixXd out2(kbetaR, kbetaR);
    Eigen::MatrixXd Ones(kbetaR, kbetaR); // identity matrix with diagonal smoothing
    out1.setZero();
    out2.setZero();
    Ones.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i));
    }
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      for (int j = 0; j < ki; j++) Ones(begin + j, begin + j) = smoothing(i);
      begin += ki;
    }
    return - out1 - out2 + Ones;
  }

  double Tri_betaR2_log_theta () {
    Eigen::MatrixXd out1(kbetaR, kbetaR);
    Eigen::MatrixXd out2(kbetaR, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d3logdensity_dmudmudlog_theta_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
      out2 += d2logdensity_dmudlog_theta_vec(i) * (mu(i) * dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i));
    }
    Eigen::MatrixXd out = - out1 - out2;
    // Eigen::Map<const Eigen::VectorXd> vA(mat.data(), mat.size());
    // Eigen::Map<const Eigen::VectorXd> vB(out.data(), out.size());
    // return vA.dot(vB);
    return (out.array() * he_s_u_mat_inv.block(kE, kE, kbetaR, kbetaR).array()).sum();
  }

  Eigen::VectorXd Tri_betaR (int nThreads) {
    Eigen::VectorXd out(kbetaR);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kbetaR; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                                    d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                                d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                                    dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                    local += tmp * dlogmu_dbetaR_mat_t_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back( - out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaR; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                out1 += tmp * dlogmu_dbetaR_mat_t_vec(i);
            }
            // out.push_back( - out1);
            out(j) = -1.0 * out1;
        }
    }

    return out;
  }

  // d^3 L / d betaR d betaR d alpha_f
  Eigen::VectorXd Tri_betaR2_alpha_f (int nThreads) {
    Eigen::VectorXd out(kE);
    double out1;
    if(nThreads > 1) {
        for (int j = 0; j < kE; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                                d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_df_mat(i,j);
                    local += tmp * dlogmu_dbetaR_mat_t_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back( - out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kE; j++) {
        out1 = 0.0;
        for (int i = 0; i < n; i++) {
            double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                        d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                        d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                        dlogdensity_dmu_vec(i) * dmu_df_mat(i,j);
            out1 += tmp * dlogmu_dbetaR_mat_t_vec(i);
        }
        // out.push_back( - out1);
        out(j) = -1.0 * out1;
    }
    }


    return out;
  }

  // d^3 L / d betaR d betaR d betaF
  Eigen::VectorXd Tri_betaR2_betaF (int nThreads) {
    Eigen::VectorXd out(kbetaF);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kbetaF; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                        d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                        d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                        dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
                    local += tmp * dlogmu_dbetaR_mat_t_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back( - out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaF; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                        d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                        d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                        dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
                out1 += tmp * dlogmu_dbetaR_mat_t_vec(i);
            }
            // out.push_back( - out1);
            out(j) = -1.0 * out1;
        }
    }

    return out;
  }

  Eigen::MatrixXd I_betaR () {  // hessian of negative likelihood without penalty
    Eigen::MatrixXd out1(kbetaR, kbetaR);
    Eigen::MatrixXd out2(kbetaR, kbetaR);
    Eigen::MatrixXd Ones(kbetaR, kbetaR); // identity matrix with diagonal smoothing
    out1.setZero();
    out2.setZero();
    Ones.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i));
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd I_betaR_i (int i) {
    Eigen::MatrixXd out1(kbetaR, kbetaR);
    Eigen::MatrixXd out2(kbetaR, kbetaR);

    out1 = d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i));

    return - out1 - out2;
  }

  Eigen::MatrixXd he_betaF () {
    Eigen::MatrixXd out1(kbetaF, kbetaF);
    Eigen::MatrixXd out2(kbetaF, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaF_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i));
    }
    return - out1 - out2;
  }


  double Tri_betaF2_log_theta () {
    Eigen::MatrixXd out1(kbetaF, kbetaF);
    Eigen::MatrixXd out2(kbetaF, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d3logdensity_dmudmudlog_theta_vec(i) * dmu_dbetaF_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += d2logdensity_dmudlog_theta_vec(i) * (mu(i) * dlogmu_dbetaF_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i));
    }
    Eigen::MatrixXd out = - out1 - out2;
    // Eigen::Map<const Eigen::VectorXd> vA(mat.data(), mat.size());
    // Eigen::Map<const Eigen::VectorXd> vB(out.data(), out.size());
    // return vA.dot(vB);
    return (out.array() * he_s_u_mat_inv.block(kE+kbetaR, kE+kbetaR, kbetaF, kbetaF).array()).sum();
  }


  Eigen::VectorXd Tri_betaF (int nThreads) {
    Eigen::VectorXd out(kbetaF);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kbetaF; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                        d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                        d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                        dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
                    local += tmp * dlogmu_dbetaF_mat_t_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaF; j++) {
        out1 = 0.0;
        for (int i = 0; i < n; i++) {
            double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                        d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                        d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                        dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
            out1 += tmp * dlogmu_dbetaF_mat_t_vec(i);
        }
        // out.push_back( - out1);
        out(j) = -1.0 * out1;
    }
    }


    return out;
  }




  // d^3 L / d betaF d betaF d alpha_f
  Eigen::VectorXd Tri_betaF2_alpha_f (int nThreads) {
    Eigen::VectorXd out(kE);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kE; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                                d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_df_mat(i,j);
                    local += tmp * dlogmu_dbetaF_mat_t_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kE; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                            d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                            dlogdensity_dmu_vec(i) * dmu_df_mat(i,j);
                out1 += tmp * dlogmu_dbetaF_mat_t_vec(i);
            }
            out(j) = -1.0 * out1;
        }
    }

    return out;
  }




  // d^3 L / d betaF d betaF d betaR
  Eigen::VectorXd Tri_betaF2_betaR (int nThreads) {
    Eigen::VectorXd out(kbetaR);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kbetaR; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                        d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                                d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                    local += tmp * dlogmu_dbetaF_mat_t_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaR; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                            d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                out1 += tmp * dlogmu_dbetaF_mat_t_vec(i);
            }
            out(j) = -1.0 * out1;
        }
    }

    return out;
  }

  Eigen::MatrixXd he_betaF_i (int i) {
    Eigen::MatrixXd out1(kbetaF, kbetaF);
    Eigen::MatrixXd out2(kbetaF, kbetaF);

    out1 = d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaF_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i));

    return - out1 - out2;
  }


  Eigen::MatrixXd he_alpha_f_betaF () {
    Eigen::MatrixXd out1(kE, kbetaF);
    Eigen::MatrixXd out2(kE, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    }
    return - out1 - out2;
  }

  double Tri_alpha_f_betaF_log_theta () {
    Eigen::MatrixXd out1(kE, kbetaF);
    Eigen::MatrixXd out2(kE, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d3logdensity_dmudmudlog_theta_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += d2logdensity_dmudlog_theta_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    }
    Eigen::MatrixXd out = - out1 - out2;
    // Eigen::Map<const Eigen::VectorXd> vA(mat.data(), mat.size());
    // Eigen::Map<const Eigen::VectorXd> vB(out.data(), out.size());
    // return vA.dot(vB);
    return (out.array() * he_s_u_mat_inv.block(0,kE+kbetaR,kE,kbetaF).array()).sum();
  }

  Eigen::VectorXd Tri_alpha_f_betaF_alpha_f (int nThreads) {
    Eigen::VectorXd out(kE);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kE; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                           d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                         d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                           dlogdensity_dmu_vec(i) * dmu_df_mat(i,j);
                    local += tmp * dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kE; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                            d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                            dlogdensity_dmu_vec(i) * dmu_df_mat(i,j);
                out1 += tmp * dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec(i);
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    }
    return out;
  }

  Eigen::VectorXd Tri_alpha_f_betaF_betaR (int nThreads) {
    Eigen::VectorXd out(kbetaR);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kbetaR; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                           d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                         d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                           dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                    local += tmp * dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaR; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                            d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                            dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                out1 += tmp * dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec(i);
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    }
    return out;
  }


  Eigen::VectorXd Tri_alpha_f_betaF_betaF (int nThreads) {
    Eigen::VectorXd out(kbetaF);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kbetaF; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                           d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                         d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                           dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
                    local += tmp * dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaF; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                            d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                            dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
                out1 += tmp * dlogmu_df_mat_t_dlogmu_dbetaF_mat_vec(i);
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    }
    return out;
  }

  Eigen::MatrixXd he_alpha_f_betaF_i (int i) {
    Eigen::MatrixXd out1(kE, kbetaF);
    Eigen::MatrixXd out2(kE, kbetaF);
    out1 = d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);

    return - out1 - out2;
  }

  Eigen::MatrixXd he_alpha_f_betaR () {
    Eigen::MatrixXd out1(kE, kbetaR);
    Eigen::MatrixXd out2(kE, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
    }
    return - out1 - out2;
  }


  double Tri_alpha_f_betaR_log_theta () {
    Eigen::MatrixXd out1(kE, kbetaR);
    Eigen::MatrixXd out2(kE, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d3logdensity_dmudmudlog_theta_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
      out2 += d2logdensity_dmudlog_theta_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
    }
    Eigen::MatrixXd out = - out1 - out2;
    // Eigen::Map<const Eigen::VectorXd> vA(mat.data(), mat.size());
    // Eigen::Map<const Eigen::VectorXd> vB(out.data(), out.size());
    // return vA.dot(vB);
    return (out.array() * he_s_u_mat_inv.block(0,kE,kE,kbetaR).array()).sum();
  }

  Eigen::VectorXd Tri_alpha_f_betaR_alpha_f (int nThreads) {
    Eigen::VectorXd out(kE);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kE; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                           d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                                 d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                                 dlogdensity_dmu_vec(i) * dmu_df_mat(i,j);
                    local += tmp * dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kE; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_df_mat(i,j);
                out1 += tmp * dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec(i);
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    }
    return out;
  }

  Eigen::VectorXd Tri_alpha_f_betaR_betaR (int nThreads) {
    Eigen::VectorXd out(kbetaR);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kbetaR; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                           d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                         d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                           dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                    local += tmp * dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaR; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                out1 += tmp * dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec(i);
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    }
    return out;
  }

  Eigen::VectorXd Tri_alpha_f_betaR_betaF (int nThreads) {
    Eigen::VectorXd out(kbetaF);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kbetaF; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                           d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                         d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                           dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
                    local += tmp * dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaF; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
                out1 += tmp * dlogmu_df_mat_t_dlogmu_dbetaR_mat_vec(i);
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    }
    return out;
  }


  Eigen::MatrixXd he_alpha_f_betaR_i (int i) {
    Eigen::MatrixXd out1(kE, kbetaR);
    Eigen::MatrixXd out2(kE, kbetaR);

    out1 = d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);

    return - out1 - out2;
  }


  Eigen::MatrixXd he_betaR_betaF () {
    Eigen::MatrixXd out1(kbetaR, kbetaF);
    Eigen::MatrixXd out2(kbetaR, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * dlogmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    }
    return - out1 - out2;
  }

  double Tri_betaR_betaF_log_theta () {
    Eigen::MatrixXd out1(kbetaR, kbetaF);
    Eigen::MatrixXd out2(kbetaR, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d3logdensity_dmudmudlog_theta_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += d2logdensity_dmudlog_theta_vec(i) * dlogmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    }
    Eigen::MatrixXd out = - out1 - out2;
    // Eigen::Map<const Eigen::VectorXd> vA(mat.data(), mat.size());
    // Eigen::Map<const Eigen::VectorXd> vB(out.data(), out.size());
    // return vA.dot(vB);
    return (out.array() * he_s_u_mat_inv.block(kE, kE+kbetaR, kbetaR, kbetaF).array()).sum();
  }


  Eigen::VectorXd Tri_betaR_betaF_alpha_f (int nThreads) {
    Eigen::VectorXd out(kE);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kE; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                           d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                         d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                           dlogdensity_dmu_vec(i) * dmu_df_mat(i,j);
                    local += tmp * dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kE; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_df_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_df_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_df_mat(i,j);
                out1 += tmp * dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec(i);
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    }
    return out;
  }

  Eigen::VectorXd Tri_betaR_betaF_betaR (int nThreads) {
    Eigen::VectorXd out(kbetaR);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kbetaR; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                           d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                         d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                           dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                    local += tmp * dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaR; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaR_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_dbetaR_mat(i,j);
                out1 += tmp * dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec(i);
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    }
    return out;
  }

  Eigen::VectorXd Tri_betaR_betaF_betaF (int nThreads) {
    Eigen::VectorXd out(kbetaF);
    double out1;
    if (nThreads > 1) {
        for (int j = 0; j < kbetaF; j++) {
            std::vector<double> partial(
                nThreads, 0.0
            );
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                double &local = partial[tid];
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                           d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                         d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                           dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
                    local += tmp * dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec(i);
                }
            }
            out1 = 0.0;
            for (int t = 0; t < nThreads; ++t) {
                out1 += partial[t];
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    } else {
        for (int j = 0; j < kbetaF; j++) {
            out1 = 0.0;
            for (int i = 0; i < n; i++) {
                double tmp = d3logdensity_dmudmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) * mu(i) +
                                d2logdensity_dmudmu_vec(i) * 2 * mu(i) * dmu_dbetaF_mat(i,j) +
                            d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat(i,j) * mu(i) +
                                dlogdensity_dmu_vec(i) * dmu_dbetaF_mat(i,j);
                out1 += tmp * dlogmu_dbetaR_mat_t_dlogmu_dbetaF_mat_vec(i);
            }
            // out.push_back(- out1);
            out(j) = -1.0 * out1;
        }
    }
    return out;
  }

  Eigen::MatrixXd he_betaR_betaF_i (int i) {
    Eigen::MatrixXd out1(kbetaR, kbetaF);
    Eigen::MatrixXd out2(kbetaR, kbetaF);
    out1 = d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * dlogmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);

    return - out1 - out2;
  }



  Eigen::MatrixXd he_alpha_f_log_smoothing_f () {
    Eigen::MatrixXd out(kE, mE);
    out.setZero();
    Eigen::VectorXd tmp(kE);
    for (int i = 0; i < mE; i++) {
      tmp.setZero();
      tmp.segment(i*kEp, kEp) = smoothing_f(i) * Sf_list.at(i) * alpha_f.segment(i*kEp, kEp);
      out.col(i) = tmp;
    }

    return out;
    // return smoothing_f * Sf * alpha_f;
  }

  Eigen::MatrixXd he_alpha_f_log_smoothing_w () {
    Eigen::MatrixXd out(kE, mE);
    out.setZero();
    Eigen::VectorXd tmp(kE);
    for (int i = 0; i < mE; i++) {
      tmp.setZero();
      tmp.segment(i*kEp, kEp) = smoothing_w(i) * Sw_list.at(i) * alpha_f.segment(i*kEp, kEp);
      out.col(i) = tmp;
    }

    return out;
    // return smoothing_w * Sw * alpha_f;
  }



  Eigen::MatrixXd he_betaR_logsmoothing () {
    Eigen::MatrixXd out(kbetaR, p);
    out.setZero();
    int begin = 0;
    for (int i = 0; i < p; i++) {
      int ki = static_cast<int>(r(i));
      for (int j = 0; j < ki; j++) out(begin + j, i) = smoothing(i) * betaR(begin + j);
      begin += ki;
    }
    return out;
  }
  Eigen::VectorXd he_alpha_f_log_theta () {
    // he_alpha_f_theta = dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec;
    return -1.0*dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }
  Eigen::VectorXd he_alpha_f_log_theta_i (int i) {
    // he_alpha_f_theta = dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec;
    return -1.0*dmu_df_mat.row(i).transpose() * d2logdensity_dmudtheta_vec(i) * theta;
  }
  Eigen::VectorXd he_betaR_log_theta () {
    return -1.0*dmu_dbetaR_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }
  Eigen::VectorXd he_betaR_log_theta_i (int i) {
    return -1.0*dmu_dbetaR_mat.row(i).transpose() * d2logdensity_dmudtheta_vec(i) * theta;
  }
  Eigen::VectorXd he_betaF_log_theta () {
    return -1.0*dmu_dbetaF_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }
  Eigen::VectorXd he_betaF_log_theta_i (int i) {
    return -1.0*dmu_dbetaF_mat.row(i).transpose() * d2logdensity_dmudtheta_vec(i) * theta;
  }





  // *********** LAML ***********
  double logdetH05() {
    double logdetH05 = 0.0;
    // Eigen::PartialPivLU<Eigen::MatrixXd> lu(he_s_u_mat);
    // Eigen::MatrixXd LU = lu.matrixLU();
    // // double c = lu.permutationP().determinant(); // -1 or 1
    // double lii;
    // for (int i = 0; i < (kE+kbetaR+kbetaF); i++) {
    //   lii = LU(i,i);
    //   logdetH05 += log(abs(lii));
    //   // logdetH05 += log(abs(lii) + 1e-8); // add small value to avoid log(0)
    //   // logdetH05 += log(std::max(lii, 1e-8));
    // }
    // logdetH05 += log(c);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(he_s_u_mat);
    Eigen::VectorXd eigvals = eigvec.eigenvalues().array();
    for (int i = 0; i < (kE+kbetaR+kbetaF); i++) {
        logdetH05 += log(abs(eigvals(i)));
    }
    return logdetH05/2.0;
  }

  // copy constructor
  Model(const Model&) = default;

};




// optimize alpha_f and betaF for a given log_smoothing_f, log_smoothing_w, and log_theta
// Use newton method with eigenvalue modification
void PL(Model& modelobj, bool verbose){
    int maxitr = 20, itr = 0;
    const double eps = 1e-05;
    double mineig = 1e-03; // minimum eigenvalue of Hessian, to ensure it is PD
    int maxstephalve = 20, stephalve = 0;
    int resetitr = 0, maxreset = 1; // if step is nan, reset coefficients as 0. only once.
    int additr = 20; // allow further iterations after resetting coefficients as 0

    Eigen::VectorXd alpha_f = modelobj.alpha_f;
    Eigen::VectorXd betaR = modelobj.betaR;
    Eigen::VectorXd betaF = modelobj.betaF;

    int kE = modelobj.kE;
    int mE = modelobj.mE;
    int kbetaR = modelobj.kbetaR;
    int kbetaF = modelobj.kbetaF;
    int converge = 0;

    int paraSize = kE+kbetaR+kbetaF;
    // Optimize ALPHA_F
    double u;
    double u_tmp;



    // update steps
    Eigen::VectorXd step(paraSize);
    step.setZero();

    Eigen::MatrixXd H(paraSize, paraSize);
    Eigen::VectorXd g(paraSize);

    g.setZero();
    H.setZero();

    // START DEFINE lanczos algorithm for smallest eigenvalue
    // Code from https://github.com/mrcdr/lambda-lanczos/blob/master/src/samples/sample4_use_Eigen_library.cpp
    // the matrix-vector multiplication routine
    auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
      auto eigen_in = Eigen::Map<const Eigen::VectorXd>(&in[0], in.size());
      auto eigen_out = Eigen::Map<Eigen::VectorXd>(&out[0], out.size());

      // eigen_out = H * eigen_in; // Easy version
      eigen_out.noalias() += H * eigen_in; // Efficient version
    };

    LambdaLanczos<double> engine(mv_mul, paraSize, false, 1); // Find 1 minimum eigenvalue
    std::vector<double> smallest_eigenvalues;
    std::vector<std::vector<double>> smallest_eigenvectors;
    double smallest_eigval; // smallest eigenvalue
    // END DEFINE lanczos for smallest eigenvalue

    // eigen decomposition
    Eigen::VectorXd eigvals(paraSize);
    eigvals.setZero();
    Eigen::VectorXd invabseigvals(paraSize);
    invabseigvals.setZero();
    // double eigval;
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(H,false); // Only values, not vectors
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(H,true); // Both values and vectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec; // Both values and vectors

    double delta = 1.0; // for New Q-Newton method
    double g_norm;

    for (int i = 0; i < paraSize; i++) g(i) = 1. + eps;

    while (itr < maxitr) {
      itr++;
      modelobj.derivative_coef();
      modelobj.derivative_he();

      g = modelobj.gr_s_u_vec;
      g_norm = g.norm();
    //   g_norm = g.lpNorm<Eigen::Infinity>();

      if(verbose) std::cout << "-- AlphaF Gradient Max: " << g.maxCoeff() << std::endl;

      if (g_norm < eps) break;
      modelobj.NegativeLogLikelihood();
      u = modelobj.NegLogL;

      H = modelobj.he_s_u_mat;
      // std::cout << "H has NaN: " << containsNaN(H) << std::endl;

      engine.run(smallest_eigenvalues, smallest_eigenvectors);
      smallest_eigval = smallest_eigenvalues[0]; // the smallest eigenvalue
      if ((smallest_eigval < 1e-2) || std::isnan(smallest_eigval)) {
        // Do Q-Newton's Step
        eigvec.compute(H, Eigen::ComputeEigenvectors); // Compute eigenvalues and vectors
        eigvals = eigvec.eigenvalues().array();

        if (abs(eigvals.prod()) < 1e-3) {
          for (int iii = 0; iii < paraSize; iii++) eigvals(iii) += delta*g_norm;
        }
        // std::cout << "eigvals" << eigvals.transpose() << std::endl;
        // for (int i = 0; i < paraSize; i++) invabseigvals(i) = 1. / max(abs(eigvals(i)), mineig); // flip signs
        for (int i = 0; i < paraSize; i++) invabseigvals(i) = 1. / abs(eigvals(i)); // flip signs
        step = eigvec.eigenvectors() * (invabseigvals.asDiagonal()) * (eigvec.eigenvectors().transpose()) * g;
      } else {
        // smallest eigenvalue > 1e-3
        // regular Newton's step
        // step = H.llt().solve(g);
        step = H.ldlt().solve(g);
      }

      // check nan in step
      // Really needed
      if(hasNaN(step)){
        if (resetitr < maxreset){
          resetitr++;
          alpha_f.setZero(); // reset alpha_f
          betaR.setZero();
          betaF.setZero();
          // alpha_f.setOnes(); // reset alpha_f
          // betaR.setOnes();
          // betaF.setOnes();
          modelobj.setAlphaF(alpha_f);
          modelobj.setBetaR(betaR);
          modelobj.setBetaF(betaF);
          if(verbose) std::cout << "reset alpha_f and betaF as 0" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          continue; // do next iteration
        } else {
          converge = 99;
          break;
        }
      }

      alpha_f -= step.segment(0, kE);
      betaR -= step.segment(kE, kbetaR);
      betaF -= step.segment(kE+kbetaR, kbetaF);

      modelobj.setAlphaF(alpha_f);
      modelobj.setBetaR(betaR);
      modelobj.setBetaF(betaF);

      modelobj.NegativeLogLikelihood();

      u_tmp = modelobj.NegLogL;

      // halving if objective function increase.
      stephalve = 0;
      while ((u_tmp > u + 1e-8) & (stephalve < maxstephalve)){
        stephalve++;
        step /= 2.;

        alpha_f += step.segment(0, kE);
        betaR += step.segment(kE, kbetaR);
        betaF += step.segment(kE+kbetaR, kbetaF);
        modelobj.setAlphaF(alpha_f);
        modelobj.setBetaR(betaR);
        modelobj.setBetaF(betaF);

        modelobj.NegativeLogLikelihood();
        u_tmp = modelobj.NegLogL;
      }


      stephalve = 0;
      // Check feasibility of step. If u is nan then we went too far;
      // halve the step and try again
      while (std::isnan(u_tmp) & (stephalve < maxstephalve)) {
        stephalve++;
        step /= 2.; // This is still the step from the previous iteration

        alpha_f += step.segment(0, kE);
        betaR += step.segment(kE, kbetaR);
        betaF += step.segment(kE+kbetaR, kbetaF);
        modelobj.setAlphaF(alpha_f);
        modelobj.setBetaR(betaR);
        modelobj.setBetaF(betaF);
        modelobj.NegativeLogLikelihood();
        u_tmp = modelobj.NegLogL;
      }

      stephalve = 0;
      // if (stephalve > 0) std::cout << "Performed " << stephalve << " iterations of step-halving." << std::endl;
      if (std::isnan(u_tmp)) {
        // Step-halving didn't work
        // std::cout << "AlphaF: Step-halving failed with nan function value. Returning failure." << std::endl;
        converge = 99;
        break;
      }
    }
    if(itr == maxitr){
      // std::cout << "AlphaF: Newton method for updating alpha fails" << std::endl;
      converge = 99;
    }


    // // check Hessian
    // H = modelobj.he_s_u_mat;
    // engine.run(smallest_eigenvalues, smallest_eigenvectors);
    // smallest_eigval = smallest_eigenvalues[0]; // the smallest eigenvalue
    // std::cout << "smallest eigenvalue of Hessian: " << smallest_eigval << std::endl;

    modelobj.PLg = g.maxCoeff();
    modelobj.converge = converge; // 0: converge. 99: not converge
}




#endif

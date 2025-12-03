/** Include **/
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

#include <random> // for generating samples from standard normal distribution

#ifdef _OPENMP
  #include <omp.h>
#endif


#include "drfDLNMadditiveeigenClass.hpp"

// #include <LBFGSB.h>
// using namespace LBFGSpp;
// TODO: use https://github.com/yixuan/LBFGSpp



std::vector<Eigen::MatrixXd> convertListToVectorEigenMatrixXd(List rList) {

  int n = rList.size();

  std::vector<Eigen::MatrixXd> matrixVector;


  for (int i = 0; i < n; ++i) {
    NumericMatrix rMatrix = as<NumericMatrix>(rList[i]);
    Eigen::Map<Eigen::MatrixXd> eigenMatrix(as<Eigen::Map<Eigen::MatrixXd>>(rMatrix));
    matrixVector.push_back(eigenMatrix);
  }

  return matrixVector;
}


std::vector<Eigen::VectorXd> convertListToVectorEigenVectorXd(List rList) {

  int n = rList.size();
  std::vector<Eigen::VectorXd> vectorList;

  for (int i = 0; i < n; ++i) {
    NumericVector rVector = as<NumericVector>(rList[i]);
    Eigen::Map<Eigen::VectorXd> eigenVector(as<Eigen::Map<Eigen::VectorXd>>(rVector));
    vectorList.push_back(eigenVector);
  }

  return vectorList;
}



struct LAMLResult {
    double fn;
    Eigen::VectorXd gradient;
};


// [[Rcpp::export]]
List drfDLNMadditivebuild(const Eigen::VectorXd R_y,
                   const List R_B_inner_list,
                   const List R_knots_f_list,
                   const List R_knots_w_list,
                   const List R_Sw_list,
                   const List R_SwR_list,
                   const List R_Sf_list,
                   const List R_SfR_list,
                   const Eigen::MatrixXd R_Xrand,
                   const Eigen::MatrixXd R_Xfix,
                   const List R_Zf_list,
                   const Eigen::VectorXd R_Xoffset,
                   const Eigen::VectorXd R_r,
                   Eigen::VectorXd R_alpha_f,
                   double R_log_theta,
                   Eigen::VectorXd R_log_smoothing_f,
                   Eigen::VectorXd R_log_smoothing_w,
                   Eigen::VectorXd R_betaR,
                   Eigen::VectorXd R_betaF,
                   Eigen::VectorXd R_logsmoothing) {
    // convert
    std::vector<Eigen::MatrixXd> R_B_inner_Eigen_list = convertListToVectorEigenMatrixXd(R_B_inner_list);
    std::vector<Eigen::VectorXd> R_knots_f_Eigen_list = convertListToVectorEigenVectorXd(R_knots_f_list);
    std::vector<Eigen::VectorXd> R_knots_w_Eigen_list = convertListToVectorEigenVectorXd(R_knots_w_list);
    std::vector<Eigen::MatrixXd> R_Sw_Eigen_list = convertListToVectorEigenMatrixXd(R_Sw_list);
    std::vector<Eigen::MatrixXd> R_SwR_Eigen_list = convertListToVectorEigenMatrixXd(R_SwR_list);
    std::vector<Eigen::MatrixXd> R_Sf_Eigen_list = convertListToVectorEigenMatrixXd(R_Sf_list);
    std::vector<Eigen::MatrixXd> R_SfR_Eigen_list = convertListToVectorEigenMatrixXd(R_SfR_list);
    std::vector<Eigen::MatrixXd> R_Zf_Eigen_list = convertListToVectorEigenMatrixXd(R_Zf_list);


    Model* modelobj_ptr = new Model(R_y, R_B_inner_Eigen_list, R_knots_f_Eigen_list, R_knots_w_Eigen_list,
                                    R_Sw_Eigen_list, R_SwR_Eigen_list, R_Sf_Eigen_list, R_SfR_Eigen_list,
                                    R_Xrand, R_Xfix, R_Zf_Eigen_list, R_Xoffset, R_r,
                                    R_alpha_f, R_log_theta, R_log_smoothing_f, R_log_smoothing_w,
                                    R_betaR, R_betaF, R_logsmoothing);
    Rcpp::XPtr<Model> ptr(modelobj_ptr);


    return List::create(Named("address.eigen") = ptr);
}




LAMLResult LAML(Model& modelobj, bool verbose, int nthreads_derivH = 1) {
    LAMLResult result;
    double u_LAML;
    int kE = modelobj.kE;
    int kEp = modelobj.kEp;
    int kbetaR = modelobj.kbetaR;
    int kbetaF = modelobj.kbetaF;
    int mE = modelobj.mE;
    int p = modelobj.p;
    Eigen::VectorXd r = modelobj.r;
    modelobj.derivative_coef();
    modelobj.derivative_he();
    modelobj.derivative_full();
    modelobj.NegativeLogLikelihood();

    Eigen::VectorXd alpha_f = modelobj.alpha_f;
    Eigen::VectorXd betaR = modelobj.betaR;
    Eigen::VectorXd betaF = modelobj.betaF;

    std::vector<Eigen::MatrixXd> Sf_list = modelobj.Sf_list;
    std::vector<Eigen::MatrixXd> Sw_list = modelobj.Sw_list;

    Eigen::VectorXd gr_s_u_vec = modelobj.gr_s_u_vec;
    Eigen::VectorXd gr_s_par_vec = modelobj.gr_s_par_vec;
    Eigen::MatrixXd he_s_u_mat = modelobj.he_s_u_mat;
    Eigen::MatrixXd he_s_par_u_mat = modelobj.he_s_par_u_mat;

    

    double log_theta = modelobj.log_theta;
    Eigen::VectorXd log_smoothing_f = modelobj.log_smoothing_f;
    Eigen::VectorXd log_smoothing_w = modelobj.log_smoothing_w;
    Eigen::VectorXd logsmoothing = modelobj.logsmoothing;


    u_LAML = modelobj.logdetH05() + modelobj.NegLogL - modelobj.n/2.0 * log(2*3.141592653589793238462643383279);

    // Eigen::PartialPivLU<Eigen::MatrixXd> lu(he_s_u_mat);
    // Eigen::MatrixXd LU = lu.matrixLU();
    // double lii = 1.0;
    // for (int i = 0; i < (kE+kbetaR+kbetaF); i++) {
    //     if(LU(i,i) > 0) {
    //         lii *= 1.0;
    //     } else {
    //         lii *= -1.0;
    //     }
    // }
    // std::cout << "Prod Diagonal element " << lii << std::endl;

    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(he_s_u_mat);
    // Eigen::VectorXd eigvals = eigvec.eigenvalues().array();
    // double lii = 1.0;
    // for (int i = 0; i < (kE+kbetaR+kbetaF); i++) {
    //     if(eigvals(i) > 0) {
    //         lii *= 1.0;
    //     } else {
    //         lii *= -1.0;
    //     }
    // }
    // std::cout << "Prod Eigenvalues sign:" << lii << std::endl;
    // std::cout << "modelobj.logdetH05()" << modelobj.logdetH05() << std::endl;
    // std::cout << "modelobj.NegLogL" << modelobj.NegLogL << std::endl;
    // std::cout << "u_LAML: " << u_LAML << std::endl;

    Eigen::VectorXd gr(1 + 2*mE + p);

    if(verbose) std::cout << "Calculating gradient of LAML..." << std::endl;
    // derivative of LAML
    modelobj.derivative_third(nthreads_derivH);
    // Eigen::MatrixXd he_s_u_mat_inv = modelobj.he_s_u_mat_inv;
    Eigen::MatrixXd he_s_u_mat_inv = he_s_u_mat.ldlt().solve(Eigen::MatrixXd::Identity(kE + kbetaR + kbetaF, kE + kbetaR + kbetaF));

    Eigen::VectorXd g_LAML(kE+kbetaR+kbetaF + 1 + 2*mE + p); // d log det H / d (alpha_f, betaR, betaF, log_theta, log_smoothing_f, log_smoothing_w, logsmoothing)
    int begin = 0; // for logsmoothing
    Eigen::MatrixXd Ones(kbetaR, kbetaR); // identity matrix with diagonal smoothing for logsmoothing
    Eigen::MatrixXd mat_3rd_tmp(kE+kbetaR+kbetaF, kE+kbetaR+kbetaF); // for smoothing parameters
    for (int j = 0; j < kE+kbetaR+kbetaF + 1 + 2*mE + p; j++) {
        mat_3rd_tmp.setZero();
        if(j < kE) {
            // alpha_f
            g_LAML(j) = 0.5 * (modelobj.Tri_alpha_f_vec(j) + 
                                modelobj.Tri_betaR2_alpha_f_vec(j) + 
                                modelobj.Tri_betaF2_alpha_f_vec(j) +
                                2.0*modelobj.Tri_alpha_f_betaR_alpha_f_vec(j) +
                                2.0*modelobj.Tri_alpha_f_betaF_alpha_f_vec(j) +
                                2.0*modelobj.Tri_betaR_betaF_alpha_f_vec(j)
                               );
        }
        if (j >= kE && j < kE + kbetaR) {
            // betaR
            int idx = j - kE;
            g_LAML(j) = 0.5 * (modelobj.Tri_alpha_f2_betaR_vec(idx) + 
                                modelobj.Tri_betaR_vec(idx) + 
                                modelobj.Tri_betaF2_betaR_vec(idx) +
                                2.0*modelobj.Tri_alpha_f_betaR_betaR_vec(idx) +
                                2.0*modelobj.Tri_alpha_f_betaF_betaR_vec(idx) +
                                2.0*modelobj.Tri_betaR_betaF_betaR_vec(idx)
                               );
            // mat_3rd_tmp.block(0, 0, kE, kE)  = modelobj.Tri_alpha_f2_betaR_list.at(idx);
            // mat_3rd_tmp.block(kE, kE, kbetaR, kbetaR) = modelobj.Tri_betaR_list.at(idx);
            // mat_3rd_tmp.block(kE+kbetaR, kE+kbetaR, kbetaF, kbetaF) = modelobj.Tri_betaF2_betaR_list.at(idx);
            // mat_3rd_tmp.block(0,kE,kE,kbetaR) = modelobj.Tri_alpha_f_betaR_betaR_list.at(idx);
            // mat_3rd_tmp.block(0,kE+kbetaR,kE,kbetaF) = modelobj.Tri_alpha_f_betaF_betaR_list.at(idx);
            // mat_3rd_tmp.block(kE, kE+kbetaR, kbetaR, kbetaF) = modelobj.Tri_betaR_betaF_betaR_list.at(idx);
            // mat_3rd_tmp.block(kE,0,kbetaR,kE) = modelobj.Tri_alpha_f_betaR_betaR_list.at(idx).transpose();
            // mat_3rd_tmp.block(kE+kbetaR,0,kbetaF,kE) = modelobj.Tri_alpha_f_betaF_betaR_list.at(idx).transpose();
            // mat_3rd_tmp.block(kE+kbetaR, kE, kbetaF, kbetaR) = modelobj.Tri_betaR_betaF_betaR_list.at(idx).transpose();
        }
        if (j >= kE + kbetaR && j < kE + kbetaR + kbetaF) {
            // betaF
            int idx = j - kE - kbetaR;
            g_LAML(j) = 0.5 * (modelobj.Tri_alpha_f2_betaF_vec(idx) + 
                                modelobj.Tri_betaR2_betaF_vec(idx) + 
                                modelobj.Tri_betaF_vec(idx) +
                                2.0*modelobj.Tri_alpha_f_betaR_betaF_vec(idx) +
                                2.0*modelobj.Tri_alpha_f_betaF_betaF_vec(idx) +
                                2.0*modelobj.Tri_betaR_betaF_betaF_vec(idx)
                               );
            // mat_3rd_tmp.block(0, 0, kE, kE) = modelobj.Tri_alpha_f2_betaF_list.at(idx);
            // mat_3rd_tmp.block(kE, kE, kbetaR, kbetaR) = modelobj.Tri_betaR2_betaF_list.at(idx);
            // mat_3rd_tmp.block(kE+kbetaR, kE+kbetaR, kbetaF, kbetaF) = modelobj.Tri_betaF_list.at(idx);
            // mat_3rd_tmp.block(0,kE,kE,kbetaR) = modelobj.Tri_alpha_f_betaR_betaF_list.at(idx);
            // mat_3rd_tmp.block(0,kE+kbetaR,kE,kbetaF) = modelobj.Tri_alpha_f_betaF_betaF_list.at(idx);
            // mat_3rd_tmp.block(kE, kE+kbetaR, kbetaR, kbetaF) = modelobj.Tri_betaR_betaF_betaF_list.at(idx);
            // mat_3rd_tmp.block(kE,0,kbetaR,kE) = modelobj.Tri_alpha_f_betaR_betaF_list.at(idx).transpose();
            // mat_3rd_tmp.block(kE+kbetaR,0,kbetaF,kE) = modelobj.Tri_alpha_f_betaF_betaF_list.at(idx).transpose();
            // mat_3rd_tmp.block(kE+kbetaR, kE, kbetaF, kbetaR) = modelobj.Tri_betaR_betaF_betaF_list.at(idx).transpose();
        }
        if (j == kE + kbetaR + kbetaF) {
            // log_theta
            g_LAML(j) = 0.5 * (modelobj.Tri_alpha_f2_log_theta_scalar + 
                                modelobj.Tri_betaR2_log_theta_scalar + 
                                modelobj.Tri_betaF2_log_theta_scalar +
                                2.0*modelobj.Tri_alpha_f_betaR_log_theta_scalar +
                                2.0*modelobj.Tri_alpha_f_betaF_log_theta_scalar +
                                2.0*modelobj.Tri_betaR_betaF_log_theta_scalar
                               );
            // mat_3rd_tmp.block(0, 0, kE, kE) = modelobj.Tri_alpha_f2_log_theta_mat;
            // mat_3rd_tmp.block(kE, kE, kbetaR, kbetaR) = modelobj.Tri_betaR2_log_theta_mat;
            // mat_3rd_tmp.block(kE+kbetaR, kE+kbetaR, kbetaF, kbetaF) = modelobj.Tri_betaF2_log_theta_mat;
            // mat_3rd_tmp.block(0,kE,kE,kbetaR) = modelobj.Tri_alpha_f_betaR_log_theta_mat;
            // mat_3rd_tmp.block(0,kE+kbetaR,kE,kbetaF) = modelobj.Tri_alpha_f_betaF_log_theta_mat;
            // mat_3rd_tmp.block(kE, kE+kbetaR, kbetaR, kbetaF) = modelobj.Tri_betaR_betaF_log_theta_mat;
            // mat_3rd_tmp.block(kE,0,kbetaR,kE) = modelobj.Tri_alpha_f_betaR_log_theta_mat.transpose();
            // mat_3rd_tmp.block(kE+kbetaR,0,kbetaF,kE) = modelobj.Tri_alpha_f_betaF_log_theta_mat.transpose();
            // mat_3rd_tmp.block(kE+kbetaR, kE, kbetaF, kbetaR) = modelobj.Tri_betaR_betaF_log_theta_mat.transpose();
        }
        if (j >= kE + kbetaR + kbetaF + 1 && j < kE + kbetaR + kbetaF + 1 + mE) {
            // log_smoothing_f
            int idx = j - (kE + kbetaR + kbetaF + 1);
            mat_3rd_tmp.block(idx*kEp, idx*kEp, kEp, kEp) = exp(log_smoothing_f(idx)) * Sf_list.at(idx);
            g_LAML(j) = 0.5 * (he_s_u_mat_inv * mat_3rd_tmp).trace();
        }
        if (j >= kE + kbetaR + kbetaF + 1 + mE && j < kE + kbetaR + kbetaF + 1 + 2*mE) {
            // log_smoothing_w
            int idx = j - (kE + kbetaR + kbetaF + 1 + mE);
            mat_3rd_tmp.block(idx*kEp, idx*kEp, kEp, kEp) = exp(log_smoothing_w(idx)) * Sw_list.at(idx);
            g_LAML(j) = 0.5 * (he_s_u_mat_inv * mat_3rd_tmp).trace();
        }
        if (j >= kE + kbetaR + kbetaF + 1 + 2*mE) {
            // logsmoothing
            Ones.setZero();
            int idx = j - (kE + kbetaR + kbetaF + 1 + 2*mE);
            // Smooth Penalty
            int ki = static_cast<int>(r(idx));
            Eigen::VectorXd betaRi(ki);
            for (int jj = 0; jj < ki; jj++) Ones(begin + jj, begin + jj) = exp(logsmoothing(idx));
            begin += ki;

            mat_3rd_tmp.block(kE, kE, kbetaR, kbetaR) = Ones;
            g_LAML(j) = 0.5 * (he_s_u_mat_inv * mat_3rd_tmp).trace();
        }
    }




    // In R: grad[-(1:(kE+kw-1))] - H.full[-(1:(kE+kw-1)),(1:(kE+kw-1))] %*% as.vector(solve(H.alpha, grad[(1:(kE+kw-1))]))
    Eigen::VectorXd g1 = g_LAML.segment(0, kE + kbetaR+kbetaF) + gr_s_u_vec;
    Eigen::VectorXd g2 = g_LAML.segment(kE+kbetaR+kbetaF, 1 + 2*mE + p) + gr_s_par_vec;
    // gr = g2 - he_s_par_u_mat * he_s_u_mat.ldlt().solve(g1);
    gr = g2 - he_s_par_u_mat * he_s_u_mat_inv * g1;
    
    if(verbose) std::cout << "Gradient calculation completed." <<  std::endl;

    result.fn = u_LAML;
    result.gradient = gr;
    return result;
}




// [[Rcpp::export]]
List drfDLNMadditiveopt(SEXP ptr,
                        Eigen::VectorXd R_alpha_f,
                        double R_log_theta,
                        Eigen::VectorXd R_log_smoothing_f,
                        Eigen::VectorXd R_log_smoothing_w,
                        Eigen::VectorXd R_betaR,
                        Eigen::VectorXd R_betaF,
                        Eigen::VectorXd R_logsmoothing,
                        int nthreads_derivH,
                        bool verbose) {


    Rcpp::XPtr<Model> modelobj_ptr(ptr);

    Model& modelobj = *modelobj_ptr;

    modelobj.setAlphaF(R_alpha_f);
    modelobj.setBetaR(R_betaR);
    modelobj.setBetaF(R_betaF);
    modelobj.setLogTheta(R_log_theta);
    modelobj.setLogSmoothingF(R_log_smoothing_f);
    modelobj.setLogSmoothingW(R_log_smoothing_w);
    modelobj.setLogsmoothing(R_logsmoothing);


    PL(modelobj, verbose);

    LAMLResult LAMLresult;
    // // Inner opt.
    // Inner(modelobj, verbose);
    // // get gr of LAML
    LAMLresult = LAML(modelobj, verbose, nthreads_derivH);
    return List::create(Named("LAML.fn") = LAMLresult.fn,
                        Named("LAML.gradient") = LAMLresult.gradient,
                        Named("alpha_f.mod") = modelobj.alpha_f,
                        Named("betaR.mod") = modelobj.betaR,
                        Named("betaF.mod") = modelobj.betaF,
                        Named("PLg") = modelobj.PLg,
                        Named("address") = modelobj_ptr,
                        Named("converge") = modelobj.converge
                        );
}



// [[Rcpp::export]]
List drfDLNMadditiveCI(SEXP ptr,
                const int Rci,
                const int rseed,
                bool ifeta,
                bool verbose,
                Rcpp::Nullable<Rcpp::NumericMatrix> R_he_input = R_NilValue) {

  Rcpp::XPtr<Model> modelobj_ptr(ptr);
  Model& modelobj = *modelobj_ptr;


  Eigen::VectorXd R_alpha_f = modelobj.alpha_f;
  Eigen::VectorXd R_betaR = modelobj.betaR;
  Eigen::VectorXd R_betaF = modelobj.betaF;

  int n = modelobj.n;
  int kE = modelobj.kE;
  int mE = modelobj.mE;
  int kbetaR = modelobj.kbetaR;
  int kbetaF = modelobj.kbetaF;
  Eigen::MatrixXd Bf_matrix = modelobj.Bf_matrix;

  int paraSize = kE+kbetaR+kbetaF;
  int paraSizefull;

  // hessian
  Eigen::MatrixXd R_he;

  // Vectors for sampling
  Eigen::VectorXd R_alpha_f_sample(kE);
  Eigen::VectorXd R_betaR_sample(kbetaR);
  Eigen::VectorXd R_betaF_sample(kbetaF);


  // Matrices to save results
  Eigen::MatrixXd alpha_f_sample_mat(Rci, kE);
  Eigen::MatrixXd betaR_sample_mat(Rci, kbetaR);
  Eigen::MatrixXd betaF_sample_mat(Rci, kbetaF);

  // components for eta
  Eigen::VectorXd R_E;
  Eigen::VectorXd eta_sample;
  Eigen::MatrixXd eta_sample_mat;
  Eigen::VectorXd eta_E_sample;
  Eigen::VectorXd eta_other_sample;
  Eigen::MatrixXd eta_E_sample_mat;
  Eigen::MatrixXd eta_other_sample_mat;


  Eigen::MatrixXd R_Xfix = modelobj.getXfix();
  Eigen::MatrixXd R_Xrand = modelobj.getXrand();
  Eigen::VectorXd R_Xoffset = modelobj.getXoffset();

  if(ifeta) {
    R_E.resize(n);
    eta_sample.resize(n);
    eta_sample_mat.resize(Rci, n);
    eta_E_sample.resize(n);
    eta_other_sample.resize(n);
    eta_E_sample_mat.resize(Rci, n);
    eta_other_sample_mat.resize(Rci, n);
  }





  paraSizefull = paraSize;

  // Joint
  if (R_he_input.isNotNull()) {
    R_he = Rcpp::as<Eigen::MatrixXd>(R_he_input.get());
  } else {
    R_he = modelobj.he_s_u_mat;
  }

  Eigen::VectorXd R_u_mod(paraSizefull);
  // Hessian
  // cholesky of inverse Hessian
  Eigen::MatrixXd R_he_u_L(paraSize, paraSize);
  Eigen::MatrixXd R_he_u_L_inv(paraSizefull, paraSize);
  Eigen::VectorXd zjoint(paraSize);
  Eigen::VectorXd samplejoint(paraSizefull);


  R_u_mod << R_alpha_f, R_betaR, R_betaF;

  // cholesky of inverse Hessian
  R_he_u_L = R_he.llt().matrixL();

  R_he_u_L_inv = (invertL(R_he_u_L)).transpose();


  // std::random_device rd;
  // std::mt19937 gen(rd());
  std::mt19937 gen(rseed);
  std::normal_distribution<> dist(0, 1);


  for(int i = 0; i < Rci; i++)
  {
    // Jointly sample
    for (int j = 0; j < paraSize; j++) {
      zjoint(j) = dist(gen);
    }
    samplejoint = R_u_mod + R_he_u_L_inv * zjoint;
    // get alpha_f
    R_alpha_f_sample = samplejoint.segment(0, kE);

    // get betaR
    R_betaR_sample = samplejoint.segment(kE, kbetaR);
    // get betaF
    R_betaF_sample = samplejoint.segment(kE+kbetaR, kbetaF);
    
    if(ifeta) {
      
      eta_E_sample = Bf_matrix * R_alpha_f_sample;
      eta_other_sample = R_Xfix * R_betaF_sample + R_Xrand * R_betaR_sample; // + R_Xoffset;
      eta_sample = eta_E_sample + eta_other_sample + R_Xoffset;
      

      eta_E_sample = eta_E_sample + eta_other_sample.mean()*Eigen::VectorXd::Ones(n);
      eta_other_sample = eta_other_sample - eta_other_sample.mean()*Eigen::VectorXd::Ones(n);

      eta_E_sample_mat.row(i) = eta_E_sample.transpose();
      eta_other_sample_mat.row(i) = eta_other_sample.transpose();
      eta_sample_mat.row(i) = eta_sample.transpose();
    }

    // save
    alpha_f_sample_mat.row(i) = R_alpha_f_sample.transpose();
    betaR_sample_mat.row(i) = R_betaR_sample.transpose();
    betaF_sample_mat.row(i) = R_betaF_sample.transpose();
  }

  // point estimate for eta
  Eigen::VectorXd eta_point(n);
  Eigen::VectorXd eta_E_point(n);
  Eigen::VectorXd eta_other_point(n);

  if(ifeta) {
    

    eta_E_point = Bf_matrix * R_alpha_f;
    eta_other_point = R_Xfix * R_betaF + R_Xrand * R_betaR;
    eta_point = eta_E_point + eta_other_point + R_Xoffset;


    eta_E_point = eta_E_point + eta_other_point.mean()*Eigen::VectorXd::Ones(n);
    eta_other_point = eta_other_point - eta_other_point.mean()*Eigen::VectorXd::Ones(n);
  }


  return List::create(
                    Named("alpha_f_sample") = alpha_f_sample_mat,
                    Named("betaR_sample") = betaR_sample_mat,
                    Named("betaF_sample") = betaF_sample_mat,
                    Named("eta_sample_mat") = eta_sample_mat,
                    Named("eta_E_sample_mat") = eta_E_sample_mat,
                    Named("eta_other_sample_mat") = eta_other_sample_mat,
                    Named("eta_point") = eta_point,
                    Named("eta_E_point") = eta_E_point,
                    Named("eta_other_point") = eta_other_point,
                    Named("Hessian_inner") = R_he);
}



// [[Rcpp::export]]
List ConditionalAICdrfDLNMadditive(SEXP ptr) {

  Rcpp::XPtr<Model> modelobj_ptr(ptr);
  Model& modelobj = *modelobj_ptr;

  modelobj.prepare_AIC();
  
  Eigen::MatrixXd R_I;
  R_I = modelobj.I_mat;

  Eigen::MatrixXd IS_mat = modelobj.IS_mat;

  Eigen::MatrixXd Vbeta = IS_mat.ldlt().solve(Eigen::MatrixXd::Identity(IS_mat.rows(), IS_mat.cols())); // solve(IS_mat)
  Eigen::MatrixXd Khat = modelobj.Khat;
  

  double l = modelobj.NegLogL_l;
 
  // the widely used version of conditional AIC proposed by Hastie and Tibshirani (1990).
  // See Wood et al. 2016 JASA
  Eigen::MatrixXd mat_AIC = Vbeta * R_I;
  // double edf1 = (2 * mat_AIC - mat_AIC * mat_AIC).trace();
  double edf_conventional = mat_AIC.trace();
  double AIC_conventional = 2.0*l + 2.0*edf_conventional;

  // proposed conditional AIC
  Eigen::MatrixXd mat_cAIC = Khat * Vbeta;
  double edf_cAIC = mat_cAIC.trace();
  double AIC_cAIC = 2.0*l + 2.0*edf_cAIC;

  return List::create(Named("l") = -1.0*l,
                        Named("edf_conventional") = edf_conventional,
                        Named("AIC_conventional") = AIC_conventional,
                        Named("edf_cAIC") = edf_cAIC,
                        Named("AIC_cAIC") = AIC_cAIC
                        );
}





// [[Rcpp::export]]
List NCVdrfDLNMadditive(SEXP ptr, const List nei_list, bool verbose = false, int nthreads = 1) {
  Rcpp::XPtr<Model> modelobj_ptr(ptr);
  Model& modelobj = *modelobj_ptr;


  int kE = modelobj.kE;
  int kbetaR = modelobj.kbetaR;
  int kbetaF = modelobj.kbetaF;
  int mE = modelobj.mE;
  int p = modelobj.p;
  int n = modelobj.n;
  Eigen::VectorXd r = modelobj.r;

  Eigen::VectorXd alpha_f = modelobj.alpha_f;
  Eigen::VectorXd betaR = modelobj.betaR;
  Eigen::VectorXd betaF = modelobj.betaF;
  double log_theta = modelobj.log_theta;




  Eigen::VectorXd beta_mod(kE+kbetaR+kbetaF+1);
  beta_mod << alpha_f, betaR, betaF, log_theta;
  Eigen::MatrixXd beta_nei(beta_mod.size(), nei_list.size());

  modelobj.prepare_AIC(); // for gunpen_nei
  Eigen::MatrixXd Hpen = modelobj.IS_mat;

  Eigen::MatrixXd Kleft = modelobj.Kleft;
  Eigen::VectorXd gnei(Kleft.rows());


  Eigen::VectorXd nei_vec;
  Eigen::MatrixXd Hunpen_nei;
  Eigen::MatrixXd Hlambdanei(kE+kbetaR+kbetaF+1, kE+kbetaR+kbetaF+1); // Hpen - Hunpen_nei

  Eigen::VectorXd Dnei(nei_list.size());
  Eigen::VectorXd Pnei(nei_list.size());


  // Eigen::VectorXd Dfull(nei_list.size());
  for (size_t i = 0; i < nei_list.size(); i++) {
    // Dfull(i) = modelobj.NegativeLogLikelihood_l_i(i);

    
    nei_vec = as<Eigen::VectorXd>(nei_list[i]);
    // std::cout << "nei_vec" << nei_vec.transpose() << std::endl;

    // compute Hunpen_nei
    Hunpen_nei = modelobj.prepare_NCV(nei_vec);
    Hlambdanei = Hpen - Hunpen_nei;

    // compute gunpen_nei from K
    gnei.setZero();
    for (size_t j = 0; j < nei_vec.size(); j++) {
      int j_int = static_cast<int>(nei_vec(j)) - 1; // R to C++ index
      gnei += Kleft.col(j_int); // col(j_int);
    }

    beta_nei.col(i) = beta_mod - Hlambdanei.completeOrthogonalDecomposition().solve(gnei);


    Eigen::VectorXd beta_i = beta_nei.col(i);

    // for NCV loss
    Dnei(i) = modelobj.get_D_i(beta_i.segment(0, kE), beta_i.segment(kE, kbetaR),
                               beta_i.segment(kE+kbetaR, kbetaF), 
                               beta_i(kE+kbetaR+kbetaF),
                               static_cast<int>(i));
    if(verbose) {
        if ((i + 1) % 100 == 0) {
            std::cout << "NCV for neighbor " << i+1 << " / " << nei_list.size() << std::endl;
        }
    }
  }


  int MCR = 1000;
  Eigen::VectorXd p_i_vec(MCR);


  // single core version
  if (nthreads == 1) {

    std::mt19937 gen(123);
    std::normal_distribution<> dist(0, 1);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(Hlambdanei,true); // Both values and vectors
    Eigen::LLT<Eigen::MatrixXd> cholSolver(Hlambdanei);
    p_i_vec.setZero();
    Eigen::MatrixXd R_he_u_L_inv;
    for (size_t i = 0; i < nei_list.size(); i++) {
      
      Eigen::VectorXd beta_i = beta_nei.col(i);
      modelobj.setAlphaF(beta_i.segment(0, kE));
      modelobj.setBetaR(beta_i.segment(kE, kbetaR));
      modelobj.setBetaF(beta_i.segment(kE+kbetaR, kbetaF));
      modelobj.setLogTheta(beta_i(kE+kbetaR+kbetaF));

      
      modelobj.derivative_coef();
      modelobj.derivative_he();
      modelobj.derivative_full();
      modelobj.prepare_IS_K();

      

      nei_vec = as<Eigen::VectorXd>(nei_list[i]);
      Hlambdanei = modelobj.IS_mat - modelobj.prepare_NCV(nei_vec);

      cholSolver.compute(Hlambdanei);
      if(cholSolver.info()!=Eigen::Success) {
        eigvec.compute(Hlambdanei); // Compute eigenvalues and vectors
        Eigen::VectorXd eigvals = eigvec.eigenvalues().array();
        Eigen::VectorXd invabseigvals(eigvals.size());
        for (int ii = 0; ii < eigvals.size(); ii++) invabseigvals(ii) = 1. / max(abs(eigvals(ii)), 1e-3);
        R_he_u_L_inv = eigvec.eigenvectors() * (invabseigvals.cwiseSqrt().asDiagonal());
        // if(verbose) {
        //   std::cout << "Warning: HLambdanei is not positive definite for neighbor " << i + 1 << " . Using eigen decomposition to compute the inverse of Cholesky factor." << std::endl;
        // }
      } else {
        Eigen::MatrixXd chol_L = cholSolver.matrixL();
        R_he_u_L_inv = invertL(chol_L).transpose();
      }

      Eigen::VectorXd zjoint(kE+kbetaR+kbetaF + 1);


      p_i_vec.setZero();
      for(int r = 0; r < MCR; r++) {
        // Jointly sample
        for (int j = 0; j < kE+kbetaR+kbetaF+1; j++) {
          zjoint(j) = dist(gen);
        }
        Eigen::VectorXd samplejoint = beta_i + R_he_u_L_inv * zjoint;
        // get alpha_f
        Eigen::VectorXd R_alpha_f_sample = samplejoint.segment(0, kE);

        // get betaR
        Eigen::VectorXd R_betaR_sample = samplejoint.segment(kE, kbetaR);
        // get betaF
        Eigen::VectorXd R_betaF_sample = samplejoint.segment(kE+kbetaR, kbetaF);

        double R_log_theta_sample = samplejoint(kE+kbetaR+kbetaF);

        double p_i = modelobj.get_p_i(R_alpha_f_sample, R_betaR_sample, R_betaF_sample, R_log_theta_sample, static_cast<int>(i));

        p_i_vec(r) = p_i;
        
      }
      // PneiMat.col(i) = p_i_vec;
      Pnei(i) = p_i_vec.mean();
      // print if is multiple of 100
      if(verbose) {
        if ((i + 1) % 100 == 0) {
          std::cout << "calculate NCV loss and predictive density: for neighbor " << i + 1 << " / " << nei_list.size() << std::endl;
        }
      }
    }
  }
  // end single core version

  // multi-thread version

  if(nthreads > 1) {

    // openMP version
    const size_t N = nei_list.size();
    std::vector<Eigen::VectorXd> nei_vecs(N);
    for (size_t i = 0; i < N; ++i) nei_vecs[i] = as<Eigen::VectorXd>(nei_list[i]);

    #ifdef _OPENMP
      omp_set_num_threads(nthreads);
    #endif

    #ifdef _OPENMP
      if(verbose) {
        std::cout << "OpenMP is ON" << std::endl;
        std::cout << "Using " << nthreads << " threads for NCV loss and predictive density calculation." << std::endl;
      }
    #else
      std::cout << "Warnings: OpenMP is OFF. Install OpenMP to enable parallel processing." << std::endl;
    #endif

    Eigen::setNbThreads(1);

    #pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(N); ++i) {
      Model local_model(modelobj);


      std::mt19937 gen(123);
      std::normal_distribution<> dist(0, 1);



      const Eigen::Ref<const Eigen::VectorXd> beta_i = beta_nei.col(i);

      local_model.setAlphaF(beta_i.segment(0, kE));
      local_model.setBetaR(beta_i.segment(kE, kbetaR));
      local_model.setBetaF(beta_i.segment(kE+kbetaR, kbetaF));
      local_model.setLogTheta(beta_i(kE+kbetaR+kbetaF));

      local_model.derivative_coef();
      local_model.derivative_he();
      local_model.derivative_full();
      local_model.prepare_IS_K();

      Eigen::MatrixXd local_Hlambdanei = local_model.IS_mat - local_model.prepare_NCV(nei_vecs[i]);
      Eigen::MatrixXd R_he_u_L_inv;

      Eigen::LLT<Eigen::MatrixXd> cholSolver(local_Hlambdanei);
      if(cholSolver.info()!=Eigen::Success) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(local_Hlambdanei,true); // Both values and vectors
        Eigen::VectorXd eigvals = eigvec.eigenvalues().array();
        Eigen::VectorXd invabseigvals(eigvals.size());
        for (int ii = 0; ii < eigvals.size(); ii++) invabseigvals(ii) = 1. / max(abs(eigvals(ii)), 1e-3);
        R_he_u_L_inv = eigvec.eigenvectors() * (invabseigvals.cwiseSqrt().asDiagonal());
        // if(verbose) {
        //   std::cout << "Warning: local_Hlambdanei is not positive definite for neighbor " << i + 1 << " . Using eigen decomposition to compute the inverse of Cholesky factor." << std::endl;
        // }
      } else {
        Eigen::MatrixXd chol_L = cholSolver.matrixL();
        R_he_u_L_inv = invertL(chol_L).transpose();
      }

      Eigen::VectorXd zjoint(kE+kbetaR+kbetaF + 1);
      Eigen::VectorXd local_p_i_vec(MCR);
      local_p_i_vec.setZero();

      for(int r = 0; r < MCR; r++) {
        // Jointly sample
        for (int j = 0; j < kE+kbetaR+kbetaF + 1; j++) {
          zjoint(j) = dist(gen);
        }
        Eigen::VectorXd samplejoint = beta_i + R_he_u_L_inv * zjoint;
        // get alpha_f
        Eigen::VectorXd R_alpha_f_sample = samplejoint.segment(0, kE);
        // get betaR
        Eigen::VectorXd R_betaR_sample = samplejoint.segment(kE, kbetaR);
        // get betaF
        Eigen::VectorXd R_betaF_sample = samplejoint.segment(kE+kbetaR, kbetaF);

        double R_log_theta_sample = samplejoint(kE+kbetaR+kbetaF);

        double p_i = local_model.get_p_i(R_alpha_f_sample, R_betaR_sample, R_betaF_sample, R_log_theta_sample, static_cast<int>(i));

        local_p_i_vec(r) = p_i;
      }
      Pnei(i) = local_p_i_vec.mean();

      // print if is multiple of 100
      if(verbose) {
        if ((i + 1) % 100 == 0) {
          std::cout << "calculate NCV loss and predictive density: for neighbor " << i + 1 << " / " << nei_list.size() << std::endl;
        }
      }
    }

  }

  if(verbose) {
      if ((nei_list.size() + 1) % 100 != 0) {
        std::cout << "calculate NCV loss and predictive density: for neighbor " << nei_list.size() << " / " << nei_list.size() << std::endl;
      }
    }


  return List::create(Named("beta_nei") = beta_nei,
                      Named("beta_mod") = beta_mod,
                      Named("Dnei") = Dnei,
                      Named("Pnei") = Pnei
                      );

}

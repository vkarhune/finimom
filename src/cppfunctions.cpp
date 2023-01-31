// #include <cmath>  // std::pow

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <RcppArmadilloExtensions/sample.h>
#include <limits>

#define pi        3.141592653589793238462643383280    /* pi */

constexpr double minusinf = std::numeric_limits<double>::lowest();

// https://github.com/oldregan/RcppNotes
// [[Rcpp::export]]
inline arma::vec arma_setdiff(arma::vec x, arma::vec y){
  arma::vec x1 = arma::unique(x);
  arma::vec y1 = arma::unique(y);

  arma::vec x_ret = x1;
  arma::uvec q1;

  if(x.is_empty() && y.is_empty()){
    return(x);
  }else{
    for (size_t j = 0; j < y1.n_elem; j++) {
      if(x_ret.is_empty()){
        break;
      }else{
        q1 = arma::find(x_ret == y1(j));
        if(q1.is_empty()){
          break;
        }
      }
      x_ret.shed_row(q1(0));
    }
    return(x_ret);
  }
}

// [[Rcpp::export]]
inline arma::vec subset_vector(arma::vec x, arma::uvec pos) {
  return x.elem(pos);
}

// [[Rcpp::export]]
inline arma::vec set_vector_vals(arma::vec x, arma::uvec pos, arma::vec vals) {
  x.elem(pos) = vals;
  return x;
}




// [[Rcpp::export]]
double LMarlik(arma::vec beta, arma::mat sematinv, arma::vec z,
                       arma::vec tau, double psi, double r, int d, arma::mat LDmat, double gval) {
  arma::mat hessian = abs(sematinv*LDmat*sematinv + arma::diagmat(6*tau/pow(beta, 4) - (r + 1)/square(beta)));

  double logdeth;
  bool success = log_det_sympd(logdeth, hessian);

  if(success){
    arma::vec grad = (r + 1)/beta - 2*tau/pow(beta, 3) + sematinv*LDmat*sematinv*beta - sematinv*z;

    // return -gval + 0.5*d*log(2*pi) - 0.5*log_det_sympd(hessian) + 0.5*arma::as_scalar(grad.t() * inv(hessian) * grad);
    // return -gval + 0.5*d*log(2*pi) - 0.5*logdeth + 0.5*arma::as_scalar(grad.t() * inv(hessian) * grad);
    return -gval + 0.5*d*log(2*pi) - 0.5*logdeth + 0.5*arma::as_scalar(grad.t() * solve(hessian, grad));
  } else {
    return minusinf;
  }
}

// [[Rcpp::export]]
inline double gfunc(arma::vec x, arma::vec z, arma::mat sematinv, arma::mat LDmat,
                    arma::vec tau, double r){
  Function lg("lgamma");

  return -0.5*r*arma::accu(log(tau)) + z.size()*(*REAL(lg(0.5*r))) +
    0.5*(r + 1)*arma::accu(log(arma::square(x))) + arma::accu(tau/(arma::square(x))) +
    arma::as_scalar(0.5*x.t()*sematinv*LDmat*sematinv*x) - arma::as_scalar(z.t()*sematinv*x);
}



// [[Rcpp::export]]
Rcpp::List posterior(Rcpp::List dat, arma::vec tau, int maxsize, double r,
                     int p, int niter, arma::vec lpriorval){

  Function lg("lgamma");

  Rcpp::List out;

  arma::vec val(niter);
  arma::vec outmodsize(niter);

  arma::vec indexvec = arma::linspace(0, p-1, p);
  arma::mat betavecmat(p, niter);

  double thresh;



  arma::vec input;

  arma::vec beta;
  arma::vec se;
  arma::vec z;

  beta = as<arma::vec>(dat["beta"]);
  se = as<arma::vec>(dat["se"]);

  z = beta/se;

  arma::uvec inds;
  arma::uvec indsplus;
  inds = index_max(abs(z));

  // arma::mat LDmatprop;



  double gval;
  arma::vec gvalpar;

  bool issympd;
  arma::mat LDmatinv;

  arma::mat sematinv;
  sematinv = arma::diagmat(1/se);

  arma::mat LDmat;
  LDmat = as<arma::mat>(dat["LDmat"]);


  input = beta.elem(inds);

  int modelsize;
  modelsize = inds.size();

  arma::mat LDmatinds = LDmat(inds, inds);
  arma::mat sematinvinds = sematinv(inds, inds);


  gvalpar = arma::inv_sympd(LDmatinds)*subset_vector(beta, inds);

  gval = gfunc(gvalpar, subset_vector(z, inds), sematinvinds, LDmatinds,
               subset_vector(tau, inds), r);

  double lm;
  double lp;

  lm = LMarlik(gvalpar, sematinvinds, subset_vector(z, inds),
               subset_vector(tau, inds), 1, r, modelsize,
               LDmatinds, gval);
  lp = lm + lpriorval[modelsize - 1];
  val[0] = lm;

  arma::vec betavec(p);

  betavec = set_vector_vals(betavec, inds, gvalpar);
  betavecmat.col(0) = betavec;

  arma::uvec indsprop;
  arma::vec indsprop1;
  arma::vec indsprop2;
  int modelsizeprop;
  int add;
  arma::vec seqvec = {0, 1, 2};
  arma::vec seqvecmin = {1, 2};
  arma::vec seqvecmax = {0, 2};

  arma::vec betaprop(p);

  double lmlnew;
  double lpnew;
  double barker;

  arma::vec xtr;
  xtr = beta - LDmat*betavec;

  arma::vec pr;
  arma::vec probs(p);
  arma::vec swapindex(1);
  arma::vec addindex(1);

  outmodsize[0] = 1;

  arma::vec zerovec;

  std::ostringstream ss;

  Rcpp::StringVector modindices(niter);
  indsplus = inds + 1;

  ss.str(std::string());
  indsplus.st().raw_print(ss);
  // get string version of vector with an end-of-line character at end
  std::string s1 = ss.str();
  // remove the end-of-line character at end
  std::string s2 = s1.substr(0, (s1.size() > 0) ? (s1.size()-1) : 0);
  modindices[0] = s2;

for(int i = 1; i < niter; ++i){

  // randomly pick add = 1, delete = 0, swap = 2
  if(modelsize == 1){
    add = as_scalar(Rcpp::RcppArmadillo::sample(seqvecmin, 1, false));
  } else {
    if(modelsize != maxsize){
      add = as_scalar(Rcpp::RcppArmadillo::sample(seqvec, 1, false));
    } else {
      add = as_scalar(Rcpp::RcppArmadillo::sample(seqvecmax, 1, false));
    }
  }



  if(add == 1){

    // xtr = beta - LDmat*betavec;
    zerovec = arma::zeros(inds.size());

    probs = set_vector_vals(arma::conv_to<arma::vec>::from(xtr), inds, zerovec);
    probs = arma::square(probs);
    probs = probs/arma::accu(probs);

    swapindex = Rcpp::RcppArmadillo::sample(indexvec, 1, false, probs);
    indsprop1 = arma::join_cols(arma::conv_to<arma::vec>::from(inds), swapindex);

    indsprop = arma::conv_to<arma::uvec>::from(indsprop1(sort_index(indsprop1)));

  } else{
    if(add == 2){

      swapindex = Rcpp::RcppArmadillo::sample(arma::conv_to<arma::vec>::from(inds), 1, false);

      pr = LDmat.col(arma::conv_to<arma::uword>::from(swapindex));

      zerovec = arma::zeros(inds.size());

      probs = set_vector_vals(xtr, inds, zerovec);
      probs = arma::square(probs);
      probs = probs/arma::accu(probs);

      addindex = Rcpp::RcppArmadillo::sample(indexvec, 1, false, probs);

      indsprop1 = arma::join_cols(arma_setdiff(arma::conv_to<arma::vec>::from(inds), swapindex), addindex);

      indsprop = arma::conv_to<arma::uvec>::from(indsprop1(sort_index(indsprop1)));


    } else {

      swapindex = Rcpp::RcppArmadillo::sample(arma::conv_to<arma::vec>::from(inds), 1, false);

      indsprop = arma::conv_to<arma::uvec>::from(arma_setdiff(arma::conv_to<arma::vec>::from(inds), swapindex));

    }
  }

  modelsizeprop = indsprop.size();

  input = beta.elem(indsprop);


  arma::mat LDmatprop = LDmat(indsprop, indsprop);

  arma::mat sematinvindsprop = sematinv(indsprop, indsprop);

  // bool issympd;
  // arma::mat LDmatinv;
  issympd = arma::inv_sympd(LDmatinv, LDmatprop);


  if(issympd){
  //if(arma::rank(LDmatprop) == modelsizeprop){
  // see script used in: https://github.com/RcppCore/RcppArmadillo/issues/257

  // gvalpar = LDmatinv*beta.elem(indsprop);
  gvalpar = solve(LDmatprop, beta.elem(indsprop));

  gval = gfunc(gvalpar, subset_vector(z, indsprop), sematinvindsprop, LDmatprop,
               subset_vector(tau, indsprop), r);

  lmlnew = LMarlik(gvalpar, sematinvindsprop, subset_vector(z, indsprop),
                   subset_vector(tau, indsprop), 1, r, modelsizeprop,
                   LDmatprop, gval);
  lpnew = lmlnew + lpriorval[modelsizeprop - 1];

  betaprop = arma::zeros(p);
  betaprop = set_vector_vals(betaprop, indsprop, gvalpar);

  barker = 1/(1+exp(lp - lpnew));
  } else {
    barker = -1;
  }

  thresh = arma::randu(1)[0];

  if(thresh < barker){
    lm = lmlnew;
    lp = lpnew;
    modelsize = modelsizeprop;
    inds = indsprop;
    betavec = betaprop;
    xtr = beta - LDmat*betavec;
  }


  // https://stackoverflow.com/questions/63594363/is-there-a-way-to-convert-armadillo-vector-to-a-string-in-c
  // .st() to transpose column vector into row vector
  // https://stackoverflow.com/questions/20731/how-do-you-clear-a-stringstream-variable
  indsplus = inds + 1;
  ss.str(std::string());
  indsplus.st().raw_print(ss);
  // get string version of vector with an end-of-line character at end
  std::string s1 = ss.str();
  // remove the end-of-line character at end
  std::string s2 = s1.substr(0, (s1.size() > 0) ? (s1.size()-1) : 0);

  // val[i] = barker;
  // val[i] = lpnew;
  val[i] = lm;
  outmodsize[i] = modelsize;
  modindices[i] = s2;
  betavecmat.col(i) = betavec;

}

 out["betavecmat"] = betavecmat;
 out["value"] = val;
 out["modsize"] = outmodsize;
 out["modindices"] = modindices;

 return out;
}

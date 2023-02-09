// #include <cmath>  // std::pow

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <RcppArmadilloExtensions/sample.h>

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

using namespace roptim;

#define pi        3.141592653589793238462643383280    /* pi */

constexpr double minusinf = std::numeric_limits<double>::lowest();

// https://github.com/oldregan/RcppNotes
// [[Rcpp::export]]
inline arma::vec arma_setdiff(arma::vec x, arma::vec y){
  arma::vec x1 = arma::unique(x);
  arma::vec y1 = arma::unique(y);
  //Rcout<<"x1:"<<x1<<endl;
  //Rcout<<"y1"<<y1<<endl;
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
double LMarlik(arma::vec beta, arma::mat sematinv, arma::vec tau,
               double psi, double r, int d, arma::mat LDmat, double gval) {
  arma::mat hessian = sematinv*LDmat*sematinv + arma::diagmat(6*tau/pow(beta, 4) - (r + 1)/square(beta));
  return -gval + 0.5*d*log(2*pi) - 0.5*log_det_sympd(hessian);
}



// [[Rcpp::export]]
double gf(arma::vec x, arma::vec z, arma::mat sematinv, arma::mat LDmat, arma::vec tau, int k, int r){
  Function lg("lgamma");

  double out;

  out = -0.5*r*arma::accu(log(tau)) + k*(*REAL(lg(0.5*r))) +
    0.5*(r + 1)*arma::accu(log(arma::square(x))) + arma::accu(tau/(arma::square(x))) +
    arma::as_scalar(0.5*x.t()*sematinv*LDmat*sematinv*x) - arma::as_scalar(z.t()*sematinv*x);

  return out;
}



class Gfunc : public Functor {
  public:
    Gfunc(const Rcpp::List pars) : pars_(pars){

    }
    double operator()(const arma::vec &x) override {
      Function lg("lgamma");
      double out;

      arma::vec z;
      arma::mat sematinv;
      arma::mat LDmat;
      arma::vec seinv;
      arma::vec tau;
      double r;

      tau = as<arma::vec>(pars_["tau"]);
      z = as<arma::vec>(pars_["z"]);
      // seinv = as<arma::vec>(pars_["seinv"]);
      // sematinv = arma::diagmat(seinv);
      sematinv = as<arma::mat>(pars_["sematinv"]);
      LDmat = as<arma::mat>(pars_["LDmat"]);
      r = pars_["r"];

      int k = z.size();

      out =  -0.5*r*arma::accu(log(tau)) + k*(*REAL(lg(0.5*r))) +
        0.5*(r + 1)*arma::accu(log(arma::square(x))) + arma::accu(tau/(arma::square(x))) +
        // 0.5*arma::as_scalar(-2*z.t()*sematinv*x + x.t()*sematinv*LDmat*sematinv*x);
        arma::as_scalar(0.5*x.t()*sematinv*LDmat*sematinv*x) - arma::as_scalar(z.t()*sematinv*x);

      return out;

    }

private:
  Rcpp::List pars_;
};


// [[Rcpp::export]]
Rcpp::List opt_nm(arma::vec &initval, Rcpp::List &pars) {
  // Rosen rb;
  Gfunc gf(pars);
  Roptim<Gfunc> opt("Nelder-Mead");
  opt.control.trace = 0;
  opt.control.warn_1d_NelderMead = false;
  opt.set_hessian(false);

  Rcpp::List out;

  opt.minimize(gf, initval);

  out["par"] = opt.par();
  out["value"] = opt.value();

  return out;

}



// [[Rcpp::export]]
double LMarlikApprox(arma::vec beta, arma::mat sematinv, arma::vec z, arma::vec tau,
                     double psi, double r, int d, arma::mat LDmat, double gval) {

  arma::mat hessian = sematinv*LDmat*sematinv + arma::diagmat(6*tau/pow(beta, 4) - (r + 1)/square(beta));

  double logdeth;
  double sign;
  // bool success = log_det_sympd(logdeth, hessian);
  bool success = log_det(logdeth, sign, hessian);

  if(success){
    //if(arma::rank(hessian) == beta.size()){
    arma::vec gr = (r + 1)/beta - 2*tau/pow(beta, 3) + sematinv*LDmat*sematinv*beta - sematinv*z;

    return -gval + 0.5*d*log(2*pi) - 0.5*sign*logdeth + 0.5*arma::as_scalar(gr.t()*inv(hessian)*gr);

  } else {

    Rcpp::List parsinit;
    Rcpp::List funlist;

    parsinit["tau"] = tau;
    parsinit["z"] = z;
    parsinit["LDmat"] = LDmat;
    parsinit["sematinv"] = sematinv;
    parsinit["r"] = r;

    funlist = opt_nm(beta, parsinit);

    double gval;
    gval = funlist["value"];

    arma::vec gvalpar;
    gvalpar = as<arma::vec>(funlist["par"]);

    double out;
    out = LMarlik(gvalpar, parsinit["sematinv"], parsinit["tau"], 1, r, gvalpar.size(),
                  parsinit["LDmat"], gval);

    return out;
    // return minusinf;
  }
}



// [[Rcpp::export]]
Rcpp::List posterior(Rcpp::List dat, arma::vec tau, int maxsize, double r,
                     int p, int niter, arma::vec lpriorval, int approx){

  Rcpp::List out;
  // arma::vec theta(niter);
  arma::vec val(niter);
  arma::vec outmodsize(niter);
  // arma::vec outtemp(niter);
  arma::vec indexvec = arma::linspace(0, p-1, p);
  arma::mat betavecmat(p, niter);
  // double current;
  double thresh;

  Rcpp::List funlist;
  Rcpp::List pars;

  // val[0] = 0;

  arma::vec input;

  arma::vec beta;
  arma::vec se;
  arma::vec z;

  beta = as<arma::vec>(dat["beta"]);
  se = as<arma::vec>(dat["se"]);

  z = beta/se;

  arma::uvec inds;
  arma::uvec indsplus;
  // inds = {182, 233};
  // inds = {0, 1}; // REMEMBER INDEXING STARTS AT ZERO!
  inds = index_max(abs(z));
  // inds = {7};


  double gval;
  arma::vec gvalpar;

  // z = {10.724232, 3.703226};

  // arma::vec seinv;
  arma::mat sematinv;
  // seinv = {1/0.01854270, 1/0.01883774};
  // sematinv = arma::diagmat(seinv);
  sematinv = arma::diagmat(1/se);

  arma::mat LDmat;
  // arma::mat LDmat = {{1.00000000, -0.08843057}, {-0.08843057, 1.00000000}};
  LDmat = as<arma::mat>(dat["LDmat"]);

  // input[0] = 0;

  // tau = {0.0083, 0.0083};

  pars["tau"] = tau;
  // pars["z"] = {10.724232, 3.703226};
  // pars["sematinv"] = arma::diagmat({1/0.01854270, 1/0.01883774});
  // pars["LDmat"] = {{1.00000000, -0.08843057}, {-0.08843057, 1.00000000}};

  pars["z"] = z;
  pars["sematinv"] = sematinv;
  pars["LDmat"] = LDmat;

  input = beta.elem(inds);
  // input = {0.19885617, 0.06976041};

  // https://arma.sourceforge.net/docs.html#submat
  // matrices: X.submat()

  Rcpp::List parsinit;

  parsinit["tau"] = subset_vector(tau, inds);
  parsinit["z"] = subset_vector(z, inds);

  arma::mat tmpmatl = LDmat(inds, inds);
  parsinit["LDmat"] = tmpmatl;

  arma::mat tmpmats = sematinv(inds, inds);
  parsinit["sematinv"] = tmpmats;

  parsinit["r"] = r;
  // arma::vec tmpvec = {1/0.01854270, 1/0.01883774};
  // arma::mat tmp = arma::diagmat(tmpvec);
  // parsinit["sematinv"] = tmp;
  // arma::mat LDmattmp = {{1.00000000, -0.08843057}, {-0.08843057, 1.00000000}};
  // parsinit["LDmat"] = LDmattmp;

  int modelsize;
  modelsize = inds.size();

  double lm;
  double lp;

  if(approx == 0){
  // optimise using proposed indices
  funlist = opt_nm(input, parsinit);

  gval = funlist["value"];
  gvalpar = as<arma::vec>(funlist["par"]);

  // theta[0] = gval;


  lm = LMarlik(gvalpar, parsinit["sematinv"], parsinit["tau"], 1, r, modelsize,
                   parsinit["LDmat"], gval);
  } else {

    gvalpar = subset_vector(beta, inds);
    gval = gf(gvalpar, subset_vector(z, inds),
              parsinit["sematinv"], parsinit["LDmat"], parsinit["tau"],
                                                               modelsize, r);

    lm = LMarlikApprox(gvalpar, parsinit["sematinv"], subset_vector(z, inds),
                       parsinit["tau"], 1, r, modelsize, parsinit["LDmat"], gval);
  }

  lp = lm + lpriorval[modelsize - 1];

  val[0] = lm;

  // out["gval"] = gval;
  // out["gvalpar"] = gvalpar;

  // int modelsize;
  // modelsize = inds.size();

  // arma::vec betavec(p, arma::fill::zeros);
  arma::vec betavec(p);

  // betavec[inds] <- gvalpar
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

  // testing for now:
  // indsprop = inds;

  arma::vec betaprop(p);

  double lmlnew;
  double lpnew;
  double barker;

  arma::vec xtr = beta - LDmat*betavec;
  arma::vec pr;
  // arma::mat xtr;
  arma::vec probs(p);
  // arma::vec probs2;
  // arma::uvec swapindex(1);
  arma::vec swapindex(1);
  arma::vec addindex(1);
  // arma::mat outbeta; ?!

  outmodsize[0] = 1;
  // outtemp[0] = 0;

  // arma::mat A;
  // arma::vec B;

  arma::vec zerovec;

  std::ostringstream ss;

  Rcpp::StringVector modindices(niter);
  indsplus = inds + 1;
  // modindices[0] = inds;
  ss.str(std::string());
  indsplus.st().raw_print(ss);
  // get string version of vector with an end-of-line character at end
  std::string s1 = ss.str();
  // remove the end-of-line character at end
  std::string s2 = s1.substr(0, (s1.size() > 0) ? (s1.size()-1) : 0);
  modindices[0] = s2;

for(int i = 1; i < niter; ++i){

  // current = theta[i-1] + rnorm(1,0, 1)[0];
  // theta[i] = current;

  // input[i] = current;
  // tau = tau + 1;
  // pars["tau"] = tau;

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

    // works:
    // A = beta - LDmat*betavec;
    // xtr = A.col(0);
    // zerovec = arma::zeros(inds.size());
    // xtr = beta - LDmat*betavec; // this moved to end
    zerovec = arma::zeros(inds.size());

    probs = set_vector_vals(arma::conv_to<arma::vec>::from(xtr), inds, zerovec);
    // probs = set_vector_vals(xtr, inds, zerovec);
    probs = arma::square(probs);
    probs = probs/arma::accu(probs);

    swapindex = Rcpp::RcppArmadillo::sample(indexvec, 1, false, probs);
    indsprop1 = arma::join_cols(arma::conv_to<arma::vec>::from(inds), swapindex);
    // indsprop1 = indsprop1(sort_index(indsprop1));
    // indsprop = arma::conv_to<arma::uvec>::from(indsprop1);

    indsprop = arma::conv_to<arma::uvec>::from(indsprop1(sort_index(indsprop1)));

  } else{
    if(add == 2){
      // outtemp[i] = add + 1;

      swapindex = Rcpp::RcppArmadillo::sample(arma::conv_to<arma::vec>::from(inds), 1, false);

      // probs = LDmat.col(arma::conv_to<arma::uword>::from(swapindex));
      pr = LDmat.col(arma::conv_to<arma::uword>::from(swapindex));

      zerovec = arma::zeros(inds.size());

      // probs = set_vector_vals(probs, inds, zerovec);
      // probs = set_vector_vals(xtr, inds, zerovec);
      probs = set_vector_vals(pr, inds, zerovec);
      probs = arma::square(probs);
      probs = probs/arma::accu(probs);

      addindex = Rcpp::RcppArmadillo::sample(indexvec, 1, false, probs);


      // this works:
      // indsprop1 = arma::join_cols(arma::conv_to<arma::vec>::from(inds), addindex);
      // indsprop1 = arma_setdiff(indsprop1, swapindex);
      // indsprop1 = indsprop1(sort_index(indsprop1));
      // indsprop = arma::conv_to<arma::uvec>::from(indsprop1);

      indsprop1 = arma::join_cols(arma_setdiff(arma::conv_to<arma::vec>::from(inds), swapindex), addindex);

      // indsprop1 = indsprop1(sort_index(indsprop1));
      // indsprop = arma::conv_to<arma::uvec>::from(indsprop1);

      indsprop = arma::conv_to<arma::uvec>::from(indsprop1(sort_index(indsprop1)));

    } else {
      // outtemp[i] = add;

      swapindex = Rcpp::RcppArmadillo::sample(arma::conv_to<arma::vec>::from(inds), 1, false);

      indsprop = arma::conv_to<arma::uvec>::from(arma_setdiff(arma::conv_to<arma::vec>::from(inds), swapindex));
      // indsprop1 = arma_setdiff(arma::conv_to<arma::vec>::from(inds), swapindex);
      // indsprop = arma::conv_to<arma::uvec>::from(indsprop1);

      // indsprop = arma::conv_to<arma::uvec>::from(arma_setdiff(arma::conv_to<arma::vec>::from(inds), swapindex));
      // swapindex <- sample2(which(betavec != 0), 1)
      // indsprop <- setdiff(which(betavec != 0), swapindex)
    }
  }

  modelsizeprop = indsprop.size();

  bool useala;

  if(approx == 0){
    parsinit["tau"] = subset_vector(tau, indsprop);
    parsinit["z"] = subset_vector(z, indsprop);

    arma::mat tmpmat1 = sematinv(indsprop, indsprop);
    parsinit["sematinv"] = tmpmat1;

    arma::mat tmpmat2 = LDmat(indsprop, indsprop);
    parsinit["LDmat"] = tmpmat2;
    // r is the same all the time

    input = beta.elem(indsprop);

    funlist = opt_nm(input, parsinit);
    gval = funlist["value"];
    gvalpar = as<arma::vec>(funlist["par"]);

    // NOTE: makes more sense to pre-calculate Lprior, and then give as input
    lmlnew = LMarlik(gvalpar, parsinit["sematinv"], parsinit["tau"], 1, r, modelsizeprop,
                     parsinit["LDmat"], gval);
    lpnew = lmlnew + lpriorval[modelsizeprop - 1];

    betaprop = arma::zeros(p);
    betaprop = set_vector_vals(betaprop, indsprop, gvalpar);

    barker = 1/(1+exp(lp - lpnew));

  } else {

    arma::mat LDmatprop = LDmat(indsprop, indsprop);

    arma::mat sematinvindsprop = sematinv(indsprop, indsprop);

    // boolean useala;

    if(modelsizeprop == 1){
      useala = true;
    } else{

      arma::mat LDmatupper = arma::abs(arma::trimatu(LDmatprop, 1));
      double mval = LDmatupper.max();

      if(mval <= 0.99){
        useala = true;
      } else {
        useala = false;
      }

    }
    // arma::mat LDmatinv;
    // bool invertible = inv_sympd(LDmatinv, LDmatprop, arma::inv_opts::allow_approx);
    // bool invertible = arma::pinv(LDmatinv, LDmatprop);

    if(useala){
    // if(invertible){
    //if(arma::rank(LDmatprop) == modelsizeprop){
      // see script used in: https://github.com/RcppCore/RcppArmadillo/issues/257

      // gvalpar = LDmatinv*beta.elem(indsprop);
      // gvalpar = solve(LDmatprop, beta.elem(indsprop));
      gvalpar = solve(LDmatprop, beta.elem(indsprop), arma::solve_opts::no_approx);

      gval = gf(gvalpar, subset_vector(z, indsprop), sematinvindsprop, LDmatprop,
                   subset_vector(tau, indsprop), modelsizeprop, r);

      lmlnew = LMarlikApprox(gvalpar, sematinvindsprop, subset_vector(z, indsprop),
                       subset_vector(tau, indsprop), 1, r, modelsizeprop,
                       LDmatprop, gval);
      lpnew = lmlnew + lpriorval[modelsizeprop - 1];

      betaprop = arma::zeros(p);
      betaprop = set_vector_vals(betaprop, indsprop, gvalpar);

      barker = 1/(1+exp(lp - lpnew));
    } else {
      // hybrid sampling!
      // NOTE: does not really work as well as hoped...
      // barker = -1;

      parsinit["tau"] = subset_vector(tau, indsprop);
      parsinit["z"] = subset_vector(z, indsprop);

      arma::mat tmpmat1 = sematinv(indsprop, indsprop);
      parsinit["sematinv"] = tmpmat1;

      arma::mat tmpmat2 = LDmat(indsprop, indsprop);
      parsinit["LDmat"] = tmpmat2;
      // r is the same all the time

      input = beta.elem(indsprop);

      funlist = opt_nm(input, parsinit);
      gval = funlist["value"];
      gvalpar = as<arma::vec>(funlist["par"]);

      // NOTE: makes more sense to pre-calculate Lprior, and then give as input
      lmlnew = LMarlik(gvalpar, parsinit["sematinv"], parsinit["tau"], 1, r, modelsizeprop,
                       parsinit["LDmat"], gval);
      lpnew = lmlnew + lpriorval[modelsizeprop - 1];

      betaprop = arma::zeros(p);
      betaprop = set_vector_vals(betaprop, indsprop, gvalpar);

      barker = 1/(1+exp(lp - lpnew));

    }
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

 // out["iter"] = niter;
 out["betavecmat"] = betavecmat;
 // out["theta"] = theta;
 out["value"] = val;
 out["modsize"] = outmodsize;
 // out["outtemp"] = outtemp;
 // out["probs"] = probs;
 // out["swapindex"] = swapindex;
 // out["ip"] = addindex;
 out["modindices"] = modindices;
 // out["betavec"] = betavec;
 // out["betaprop"] = betaprop;
 // out["betavec"] = betavec;
 // out["sematinv"] = parsinit["sematinv"];
 // out["LDmat"] = parsinit["LDmat"];



 return out;
}

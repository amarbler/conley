#include <iostream>
//#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]

#define ARMA_64BIT_WORD 1

//#include <cmath>
//#include <algorithm>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// ------------------------------------------

// geo constants and safe acos() to avoid rounding error induced NaNs

namespace geoc {
const double DE2RA             = 0.01745329252;
const double RA2DE             = 57.2957795129;
const double ERAD              = 6378.135;
const double ERADM             = 6378135.0;
const double AVG_ERAD          = 6371.0;
const double FLATTENING        = 1.0/298.257223563;
// Earth flattening (WGS '84)
const double EPS               = 0.000000000005;
const double KM2MI             = 0.621371;
const double GEOSTATIONARY_ALT = 35786.0;    // km
}

double SafeAcos(double x) {
  if(x < -1.0) x = -1.0;
  else if(x > 1.0) x = 1.0;
  return acos(x);
}


// euclidean distance on vector, use for chord distance
template <typename InputIterator1, typename InputIterator2>
inline double vectorDistance(InputIterator1 begin1, InputIterator1 end1,
                            InputIterator2 begin2) {

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;

  double ret = 0.0;
  while (it1 != end1) {
    double dist = (*it1++) - (*it2++);
    ret += dist * dist;
  }
  return ret > 0.0 ? sqrt(ret) : 0.0;
}


// just for testing
// [[Rcpp::export]]
double vecdist(NumericVector a, NumericVector b){

  double res;
  res = vectorDistance(a.begin(), a.end(),b.begin());
  return res;

}
// ------------------------------------------

// Chord distance:
double chord_cpp(double lat1, double lon1,
                 double lat2, double lon2){

  lat1 *= geoc::DE2RA;
  lon1 *= geoc::DE2RA;
  lat2 *= geoc::DE2RA;
  lon2 *= geoc::DE2RA;

  double  R = geoc::AVG_ERAD;
  std::vector<double> u1 = {R * cos(lat1) * cos(lon1), R * cos(lat1) * sin(lon1), R * sin(lat1)};
  std::vector<double> u2 = {R * cos(lat2) * cos(lon2), R * cos(lat2) * sin(lon2), R * sin(lat2)};

  double dist = vectorDistance(u1.begin(), u1.end(),u2.begin());
  return dist;

}


// Spherical law of cosines:
double spherical_cpp(double lat1, double lon1,
                     double lat2, double lon2){

    lat1 *= geoc::DE2RA;
    lon1 *= geoc::DE2RA;
    lat2 *= geoc::DE2RA;
    lon2 *= geoc::DE2RA;

    double d = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon1 - lon2);
    return (geoc::AVG_ERAD * SafeAcos(d));
}

// Haversine distance:
double haversine_cpp(double lat1, double lon1,
                     double lat2, double lon2) {

  lat1 *= geoc::DE2RA;
  lon1 *= geoc::DE2RA;
  lat2 *= geoc::DE2RA;
  lon2 *= geoc::DE2RA;

  double delta_phi = lat2 - lat1;
  double delta_lambda = lon2 - lon1;
  double term1 = pow(sin(delta_phi / 2), 2);
  double term2 = cos(lat1) * cos(lat2) * pow(sin(delta_lambda/2), 2);
  double the_terms = term1 + term2;
  double delta_sigma = 2 * atan2(sqrt(the_terms), sqrt(1-the_terms));

  return (geoc::AVG_ERAD * delta_sigma);

}


// Distance Formula used by SH:
double sh_cpp(double lat1, double long1,
              double lat2, double long2) {
  double distance = pow(pow(111 * (lat1 - lat2), 2) + pow(cos(lat1 * M_PI / 180) * 111 * (long1 - long2), 2), .5);
  return(distance);
}

// ------------------------------------------

// Easy to write in par, just call distmat par, then parallel sandwich from cpp

// [[Rcpp::export]]
arma::mat XeeXhC(arma::mat M, double cutoff,
                 arma::mat X, arma::vec e, int n1, int k,
                 std::string kernel="bartlett",
                 std::string dist_fn="haversine"){

  long long int nrow = M.n_rows;
  arma::mat dmat(nrow, nrow, fill::zeros);

  for( long long int i = 0; i < nrow; i++ ){
    dmat(i, i) = 1;

    // removed i+1 same issue and below, think about this later
    for( long long int j = i; j < nrow; j++ ){
      double d;
      // Distance functions:
      if(dist_fn == "haversine") d = haversine_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else if(dist_fn == "spherical") d = spherical_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else if(dist_fn == "chord") d = chord_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else if(dist_fn == "flatearth") d = sh_cpp(M(i,0), M(i,1), M(j,0), M(j,1));

      // Kernel:
      int v = d <= cutoff;
      if(kernel != "bartlett") dmat(i,j) = dmat(j,i) = v;
      else dmat(i,j) = dmat(j,i) = (1 - d / cutoff) * v;
    }
  }

  arma::mat XeeXh(k, k, fill::zeros);
  for( long long int i = 0; i < nrow; i++ ){
    arma::mat e_mat(1, n1, fill::zeros);
    e_mat.fill(e[i]);

    arma::mat k_mat(k, 1, fill::ones);

    arma::mat d_row(1, n1, fill::ones);
    d_row %= dmat.row(i); d_row %= e.t();

    arma::mat X_row(k, 1, fill::ones);
    X_row %= X.row(i).t();

    XeeXh += (X_row * e_mat % (k_mat * d_row)) * X;
  }
  return XeeXh;
}


// ------------------------------------------

// [[Rcpp::export]]
arma::mat Bal_XeeXhC(arma::mat dmat,
                     arma::mat X, arma::vec e, int n1, int k){

  long long int nrow = dmat.n_rows;

  arma::mat XeeXh(k, k, fill::zeros);
  for( long long int i = 0; i < nrow; i++ ){
    arma::mat e_mat(1, n1, fill::zeros);
    e_mat.fill(e[i]);

    arma::mat k_mat(k, 1, fill::ones);

    arma::mat d_row(1, n1, fill::ones);
    d_row %= dmat.row(i); d_row %= e.t();

    arma::mat X_row(k, 1, fill::ones);
    X_row %= X.row(i).t();

    XeeXh += (X_row * e_mat % (k_mat * d_row)) * X;
  }
  return XeeXh;
}


// ------------------------------------------

// [[Rcpp::export]]
arma::mat XeeXhC_Lg(arma::mat& M, double& cutoff,
                    arma::mat& X, arma::vec& e, int&  n1, int& k,
                    std::string& kernel,
                    std::string& dist_fn){

  long long int nrow = M.n_rows;
  arma::mat XeeXh(k, k, fill::zeros);

  for( long long int i = 0; i < nrow; i++ ){
   // arma::vec d_row(nrow);
    arma::mat d_row(1, nrow, fill::ones);

    for( long long int j = 0; j < nrow; j++ ){
      double d;

      // Distance functions:
      if(dist_fn == "haversine") d = haversine_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else if(dist_fn == "spherical") d = spherical_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else if(dist_fn == "chord") d = chord_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else if(dist_fn == "flatearth") d = sh_cpp(M(i,0), M(i,1), M(j,0), M(j,1));


      // Kernel:
      int v = d <= cutoff;
      if(kernel != "bartlett") d_row[j] = v;
      else d_row[j] = (1 - d / cutoff) * v;
    }

    arma::mat e_mat(1, nrow, fill::zeros);
    e_mat.fill(e[i]);

    arma::mat k_mat(k, 1, fill::ones);

    d_row %= e.t();

    arma::mat X_row(k, 1, fill::ones);
    X_row %= X.row(i).t();

    XeeXh += (X_row * e_mat % (k_mat * d_row)) * X;
  }
  return XeeXh;
}


// ------------------------------------------


// [[Rcpp::export]]
arma::mat TimeDist(arma::vec times, double cutoff,
                   arma::mat X, arma::vec e, int n1, int k){

  long long int nrow = times.n_elem;
  arma::mat dmat(nrow, nrow, fill::ones);
  // dmat.each_row() %= times.t();

  for( long long int i = 0; i < nrow; i++ ){
    arma::vec t_diff(nrow);
    t_diff = times;

    t_diff -= times[i];
    t_diff = abs(t_diff);

    NumericVector v1(nrow); NumericVector v2(nrow);
    for( long long int j = 0; j < nrow; j++ ) {
      v1[j] = t_diff[j] <= cutoff;
      v2[j] = t_diff[j] != t_diff[i];
      t_diff[j] = v1[j] * v2[j] * (1 - t_diff[j] / (cutoff + 1));
    }

    dmat.row(i) %= t_diff.t();
  }

  arma::mat XeeXh(k, k, fill::zeros);
  for( long long int i = 0; i < nrow; i++ ){
    arma::mat e_mat(1, n1, fill::zeros);
    e_mat.fill(e[i]);

    arma::mat k_mat(k, 1, fill::ones);

    arma::mat d_row(1, n1, fill::ones);
    d_row %= dmat.row(i); d_row %= e.t();

    arma::mat X_row(k, 1, fill::ones);
    X_row %= X.row(i).t();

    XeeXh += (X_row * e_mat % (k_mat * d_row)) * X;
  }

  return XeeXh;
}


// [[Rcpp::export]]
arma::mat DistMat(arma::mat M, double cutoff,
                  std::string kernel="bartlett",
                  std::string dist_fn="haversine"){

  long long int nrow = M.n_rows;
  arma::mat dmat(nrow, nrow, fill::zeros);

  for( long long int i = 0; i < nrow; i++ ){
    dmat(i, i) = 1;

    for( long long int j = i+1; j < nrow; j++ ){
      double d;

      // Distance functions:
      if(dist_fn == "haversine") d = haversine_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else if(dist_fn == "spherical") d = spherical_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else if(dist_fn == "chord") d = chord_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else d = sh_cpp(M(i,0), M(i,1), M(j,0), M(j,1));

      // Kernel:
      int v = d <= cutoff;

      if(kernel != "bartlett") dmat(i,j) = dmat(j,i) = v;
      else  dmat(i,j) = dmat(j,i) =  (1 - d / cutoff) * v;
    }
  }
  return dmat;
}


arma::mat DistMatPure(arma::mat M, double cutoff,
                  std::string dist_fn="haversine"){

  long long int nrow = M.n_rows;
  arma::mat dmat(nrow, nrow, fill::zeros);

  for( long long int i = 0; i < nrow; i++ ){
    dmat(i, i) = 1;

    for( long long int j = i+1; j < nrow; j++ ){
      double d;

      // Distance functions:
      if(dist_fn == "haversine") d = haversine_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else if(dist_fn == "spherical") d = spherical_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else if(dist_fn == "chord") d = chord_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      else d = sh_cpp(M(i,0), M(i,1), M(j,0), M(j,1));

      dmat(i,j) = dmat(j,i) =  d;
    }
  }
  return dmat;
}

struct ParDistance: public Worker {

  // input matrix to read from
  const arma::mat& mat;
  const double& cutoff;
  std::string& kernel, dist_fn;

  // output matrix to write to
  arma::mat& dmat;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  ParDistance(const arma::mat& mat, arma::mat& dmat,  const double& cutoff,
             std::string& kernel, std::string& dist_fn)
    : mat(mat), dmat(dmat), cutoff(cutoff), kernel(kernel), dist_fn(dist_fn) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    // needs to be declared here to avoid strange numerical issues
    double d = 0;
    double v = 1;

    for (std::size_t i = begin; i < end; i++) {
   //   dmat(i,i) = 1;
   // used to be (i+1) to not calculate the diagonal, but that somehow gives wrong results with the spherical dists
   // think about this later
      for (std::size_t j = i; j < mat.n_rows; j++) {

        if(dist_fn == "haversine") d = haversine_cpp(mat(i,0), mat(i,1), mat(j,0), mat(j,1));
        else if(dist_fn == "spherical") d = spherical_cpp(mat(i,0), mat(i,1), mat(j,0), mat(j,1));
        else if(dist_fn == "chord") d = chord_cpp(mat(i,0), mat(i,1), mat(j,0), mat(j,1));
        else d = sh_cpp(mat(i,0), mat(i,1), mat(j,0), mat(j,1));

        // truncated kernel
        v = (d <=  cutoff);
        // optional: apply bartlett
        if( kernel == "bartlett") v = (1 - d / cutoff) * v;

        // write to output matrix
        dmat(i,j) = dmat(j,i) =  v;
      }
    }
  }
};

// [[Rcpp::export]]
arma::mat DistMatPar(arma::mat& mat, double& cutoff,
                     std::string kernel,
                     std::string dist_fn) {

  // allocate the matrix we will return
  arma::mat dmat(mat.n_rows, mat.n_rows, fill::zeros);

  // create the worker
  ParDistance parDistance(mat, dmat, cutoff, kernel, dist_fn);

  // call it with parallelFor
  parallelFor(0, mat.n_rows, parDistance);

  return dmat;
}

//---------- try again now for XeeX_LG

struct XeeXhC_Lg_Distance: public Worker {

  // input matrices and so on to read from
  const arma::mat&  mat, X;
  const arma::vec&  e;
  const double& cutoff;
  std::size_t k, n;
  std::string kernel, dist_fn;

  // output mat
  arma::mat XeeXh;

  // two constructors for parallelReduce
  // initialize from Rcpp input and output matrixes (mostly arma types)
  XeeXhC_Lg_Distance(const arma::mat& mat, const double& cutoff,
                     const arma::mat& X, const arma::vec& e, std::size_t k, std::size_t n,
                     std::string& kernel, std::string& dist_fn)
    : mat(mat), cutoff(cutoff), X(X), e(e), k(k), n(n),
      kernel(kernel), dist_fn(dist_fn), XeeXh(arma::mat(k,k)) { XeeXh.zeros(); }
  XeeXhC_Lg_Distance(const XeeXhC_Lg_Distance& billy, Split)
    : mat(billy.mat), cutoff(billy.cutoff),
      X(billy.X), e(billy.e), k(billy.k), n(billy.n),
      kernel(billy.kernel), dist_fn(billy.dist_fn), XeeXh(arma::mat(k,k)) { XeeXh.zeros(); }

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    // have to initialize this here! not above.
    double d = 0;
    double v = 1;
    for (std::size_t i = begin; i < end; i++) {
      arma::mat d_row(1, n, fill::ones);

      for (std::size_t j = 0; j < mat.n_rows; j++) {

        // calculate distances only off the lower diagonal
        if(dist_fn == "haversine") d = haversine_cpp(mat(i,0), mat(i,1), mat(j,0), mat(j,1));
        else if(dist_fn == "spherical") d = spherical_cpp(mat(i,0), mat(i,1), mat(j,0), mat(j,1));
        else if(dist_fn == "chord") d = chord_cpp(mat(i,0), mat(i,1), mat(j,0), mat(j,1));
        else d = sh_cpp(mat(i,0), mat(i,1), mat(j,0), mat(j,1));

        // truncated kernel
        v = (d <=  cutoff);
        // apply bartlett if option is passed
        if( kernel == "bartlett" )  v = (1 - d / cutoff) * v;

        // write to output vector
        d_row[j] = v;
      }


      arma::mat e_mat(1, n, fill::zeros);
      e_mat.fill(e[i]);

      arma::mat k_mat(k, 1, fill::ones);
      d_row %= e.t();

      arma::mat X_row(k, 1, fill::ones);
      X_row = X.row(i).t();

      XeeXh += (X_row * e_mat % (k_mat * d_row)) * X;

    }
  }
  // the join
  void join(const XeeXhC_Lg_Distance& rhs) {
    XeeXh += rhs.XeeXh;
  }
};

// [[Rcpp::export]]
arma::mat XeeXhC_Lg_Par(arma::mat & mat,
                        double& cutoff,
                        arma::mat & X,
                        arma::vec & e,
                        std::string kernel,
                        std::string dist_fn) {

  // dim sizes
  std::size_t k = X.n_cols;
  std::size_t n = X.n_rows;

  // create the worker
  XeeXhC_Lg_Distance XeeXhC_lg_Distance(mat, cutoff, X, e, k, n, kernel, dist_fn);

  // call it with parallelReduce
  parallelReduce(0, mat.n_rows, XeeXhC_lg_Distance);

  return XeeXhC_lg_Distance.XeeXh;
}


// test this function with the correct matrices
struct ParSandwich: public Worker {

  // input matrices and so on to read from
  const arma::mat&  dmat, X;
  const arma::vec&  e;
  std::size_t k, n;

  // output mat
  arma::mat output;

  // two constructors for parallelReduce
  // initialize from Rcpp input and output matrixes (mostly arma types)
  ParSandwich(const arma::mat& dmat, const arma::mat& X, const arma::vec& e,
                  std::size_t k, std::size_t n) : dmat(dmat), X(X), e(e), k(k), n(n),
                  output(arma::mat(k,k)) { output.zeros(); }
  ParSandwich(const ParSandwich& billy, Split) : dmat(billy.dmat), X(billy.X), e(billy.e),
                  k(billy.k), n(billy.n),
                  output(arma::mat(k,k)) { output.zeros(); }

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t i = begin; i < end; i++) {

      arma::mat e_mat(1, n, fill::zeros);
      e_mat.fill(e[i]);

      arma::mat k_mat(k, 1, fill::ones);

      arma::mat d_row(1, n, fill::ones);
      d_row %= dmat.row(i); d_row %= e.t();

      arma::mat X_row(k, 1, fill::ones);
      X_row %= X.row(i).t();

      output += (X_row * e_mat % (k_mat * d_row)) * X;

    }
  }

  // the join
  void join(const ParSandwich& rhs) {
    output += rhs.output;
  }
};


// [[Rcpp::export]]
arma::mat Bal_XeeXhC_Par(arma::mat & dmat,
                         arma::mat & X,
                         arma::vec & e) {

  std::size_t k = X.n_cols;
  std::size_t n = X.n_rows;

  // create the worker
  ParSandwich parSandwich(dmat, X, e, k, n);

  // call it with parallelReduce
  parallelReduce(0, dmat.n_rows, parSandwich);

  return parSandwich.output;
}

// easy both components already exist as parallel code

// [[Rcpp::export]]
arma::mat XeeXhC_Par(arma::mat& mat, double& cutoff,
                     arma::mat& X, arma::vec& e,
                     std::string& kernel,
                     std::string& dist_fn){

  std::size_t n;
  std::size_t k;

  n  = mat.n_rows;
  k  = mat.n_cols;

  arma::mat dmat(n, n, fill::zeros);
  arma::mat XeeXh(k, k, fill::zeros);

  dmat = DistMatPar(mat, cutoff, kernel, dist_fn);
  XeeXh = Bal_XeeXhC_Par(dmat, X, e);

  return XeeXh;

}


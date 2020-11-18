#include <RcppArmadillo.h>
#include <iostream>
//#include <htool/include/htool/htool.hpp>





using namespace arma;
using namespace Rcpp;

////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat tps(int ep, arma::mat r){
  
  arma::mat rbf(r.n_rows, r.n_cols);
  arma::uvec rGt0 = arma::find(r > 0);
  rbf.elem(rGt0) = arma::pow(ep*r.elem(rGt0),2) % log(ep*r.elem(rGt0));
  
  return rbf;
}
//_____________________________________________________________



// SINGULAR VALUE DESCOMPOSITION (SVD)
/////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::vec baseSVD(const arma::mat & X) {
  arma::mat U; 
  arma::mat V;
  arma::vec S;
  arma::svd(U, S, V, X, "standard");
  return S;
}

// [[Rcpp::export]]
arma::vec dcSVD(const arma::mat & X) {
  arma::mat U;
  arma::mat V;
  arma::vec S;
  arma::svd(U, S, V, X, "dc");
  return S;
}
//_____________________________________________________________


////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat f1(arma::mat x, arma::mat y){
  
  arma::mat value(size(x));
  
  value = 0.75 * exp(- (pow((9*x - 2), 2) + pow((9*y - 2), 2)) / 4);
  
  return value;
}

////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat f2(arma::mat x, arma::mat y){
  
  arma::mat value(size(x));
  
  value = 0.75 * exp(- (pow((9*x + 1), 2) / 49 + pow((9*y + 1), 2) / 10));
  
  return value;
}

////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat f3(arma::mat x, arma::mat y){
  
  arma::mat value(size(x));
  
  value = 0.5 * exp(- (pow((9*x - 7), 2) + pow((9*y - 3), 2)) / 4);
  
  return value;
}

////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat f4(arma::mat x, arma::mat y){
  
  arma::mat value(size(x));
  
  value = 0.2 * exp(- (pow((9*x - 4), 2) + pow((9*y - 7), 2)));
  
  return value;
}

////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat testfunction(arma::mat x, arma::mat y){
  
  arma::mat value(size(x));
  
  value = f1(x, y) + f2(x, y) + f3(x, y) - f4(x, y);
  
  return value;
}
//______________________________________________________________________________________




////////////////////////////////////////////////////////////
// [[Rcpp::export]]
void meshgrid(arma::mat & x, arma::mat & y, arma::vec & xv, arma::vec &yv){
  
  /* Copia los valores de xv en cada fila de x, también copia los valores de 
   yv en cada columna de y. 
   x es una matriz de tamaño (y.n_elem, x.n_elem)
   y es una matriz de tamaño (y_n_elem, x.n_elem)
   */
  
  int i;
  
  const int xv_l = xv.n_elem;
  const int yv_l = yv.n_elem;
  
  /* Copia de los valores de x */
  
  for(i=0; i<xv_l; i++){
    x.col(i).fill(as_scalar(xv(i)));
  }
  
  /* Copia de los valores de y */
  for(i=0; i<yv_l; i++){
    y.row(i).fill(as_scalar(yv(i)));
  }
  
  return;
}
//_____________________________________________________________________________



////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat radialFunction(arma::mat & r, const int RBFtype, const double R){
  
  r = r / R;
  
  arma::mat phi(r.n_rows, r.n_cols);
  
  if(RBFtype == 1)               // linear spline (R1)
    phi = r;
  
  if(RBFtype == 2){               // Thin Plate Spline (TPS2??)
    arma::uvec rGt0 = arma::find(r > 0);
    phi.elem(rGt0) = arma::pow(r.elem(rGt0),2) % log(r.elem(rGt0));
  }
  
  if(RBFtype == 3)               // Gaussian (GS)
    phi = exp(-arma::pow(r,2));
  
  if(RBFtype == 4)               // Multiquadric (MQ)
    phi = sqrt(1 + arma::pow(r,2));
  
  if(RBFtype == 5)               // Inverse multiquadric (IMQ)
    phi = 1 / sqrt(1 + arma::pow(r,2));     
  
  // if(RBFtype == 6 & r < 1)      // Compact support
  //   phi = arma::pow((1 - r), 2);
  // 
  
  return(phi);
}
//_____________________________________________________________________________

////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat DistanceMatrix(const arma::mat dsites, const arma::mat ctrs){

  int M  = dsites.n_rows;
  int dim = dsites.n_cols;
  int N   = ctrs.n_rows;
  
  arma::mat DM_data(M, N, arma::fill::zeros);
  
  //if ((x.n_cols == dim) && (fPar.n_rows == Ns + dim +1))
  
  for (int i = 0; i < M; i++){
    for (int j = 0; j < N; j++){
      DM_data(i, j) = norm(dsites.row(i) - ctrs.row(j));
    }
  }
  return(DM_data);
}
//_______________________________________________________________________________




////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat DistanceMatrix2(const arma::mat dsites, const arma::mat ctrs){
  
  int d;
  int i;
  
  const int M = dsites.n_rows;
  const int N = ctrs.n_rows;
  
  /* Se asume que dsites y ctrs tienen el mismo numero de columnas */
  const int s = dsites.n_cols;
  
  arma::mat DM_data(M, N, arma::fill::zeros);
  
  arma::mat dr(M, N);
  arma::mat cc(M, N);
  
  
  for(d=0; d<s; d++){
    
    for(i=0; i<N; i++){
      dr.col(i) = dsites.col(d);
    }
    
    for(i=0; i<M; i++){
      cc.row(i) = arma::trans(ctrs.col(d));
    }
    
    DM_data += arma::pow(dr - cc, 2);
  }
  
  DM_data = arma::sqrt(DM_data);
  
  return(DM_data);
}



/////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List RBF_LS(arma::mat dsites, arma::mat ctrs, const int RBFtype, const double R, const int neval){

  int N   = dsites.n_rows;
  int dim = dsites.n_cols;
  int M   = ctrs.n_rows;
  int ncol = 1 + dim;
  
arma::mat DM_data = DistanceMatrix2(dsites, ctrs);

// CM de linea 14 de programa 19.2
arma::mat CM = radialFunction(DM_data, RBFtype, R);          // Matriz E del libro Green and Silverman

// PM de linea 15 programa 19.2 (mantiene el mismo nombre)
arma::mat PM = join_rows(ones(N,1), dsites);                 // termino polinomial (Matriz T) 
arma::mat PtM = trans(join_rows(ones(M,1), ctrs));           // matriz transpuesta termino polinomial (T transpuesta)

// CM de linea 16 programa 19.2
arma::mat A = join_cols(join_rows(CM,PM), join_rows(PtM, zeros(ncol, ncol)));
//arma::mat A = join_cols(join_rows(CM,PM), join_rows(PtM, zeros(3,3))); 


// Simulacion de una variable respuesta cins distribucion Normal
// linea 17 de programa 19.2
arma::colvec rhs = testfunction(dsites.col(0), dsites.col(1));

// linea 18 de programa 19.2 (lo nombro como rhs2 ya que no se puede renombrar en Rcpparmadillo)
arma::colvec rhs2 = rhs + 0.03*randn(rhs.n_rows, rhs.n_cols);

// linea 19 del programa 19.2 
// Vector de variable respuesta con fila 1 = rhs3 y fila 2 matriz 3x3 de ceros
arma::colvec rhs3 = join_cols(rhs2,zeros(3, rhs.n_cols));

// Creo una neval x neval igualmente espaciados sitios de evaluacion
// linea 20
arma::vec grid = linspace(0, 1, neval);

arma::mat x(neval, neval);
arma::mat y(neval, neval);
arma::mat epoints(neval*neval, neval*neval);

meshgrid(x, y, grid, grid);

epoints = join_rows(vectorise(x), vectorise(y));

// Calcular matriz de distancia entre los puntos de evaluacion y los centros
// linea 22 del programa 19.2
arma::mat DM_eval = DistanceMatrix2(epoints, ctrs);

// Matriz de evaluacion para datos de evaluacion
// EM de linea 23 programa 19.2
   arma::mat CMe = radialFunction(DM_eval, RBFtype, R); 
//arma::mat CMe = tps(ep, DM_eval); 


// PM de linea 24 del programa 19.2
arma::mat PMe = join_rows(ones(neval*neval, 1), epoints);

// EM de linea 24 lado derecho programa 19.2 (aparece como EM = [EM PM])
arma::mat CMe2 = join_rows(CMe, PMe);

// Calculo de RBF por minimos cuadrados 
// linea 25 del programa 19.2
arma::colvec Pf = CMe2 * solve(A, rhs3);

// Calcular solucíon exacta, por ejemplo: evaluar 'testfunction' en los puntos de evaluacion
// linea 26 programa 19.2
arma::colvec exact = testfunction(epoints.col(0), epoints.col(1));

// Calcular el error maximo en la grilla
// linea 27 del programa 19.2
double maxerr = norm(Pf - exact, "inf");

return Rcpp::List::create(Rcpp::Named("phi_coefficients") = CM,
                          Rcpp::Named("Pf")       = Pf,
                          Rcpp::Named("maxerr")  = maxerr,
                          Rcpp::Named("Distances")  = DM_data);

}







/////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List RBF_LS2(arma::vec y, arma::mat dsites, arma::mat ctrs, int RBFtype, const double R, const int neval){
  
  int N   = dsites.n_rows;
  int dim = dsites.n_cols;
  int M   = ctrs.n_rows;
  int ncol = 1 + dim;
  
  arma::mat DM_data = DistanceMatrix2(dsites, ctrs);
  
  // CM de linea 14 de programa 19.2
  arma::mat CM = radialFunction(DM_data, RBFtype, R);          // Matriz E del libro Green and Silverman
  
  //arma::mat CM = tps(ep, DM_data);          // Matriz E del libro Green and Silverman
  
  
  // PM de linea 15 programa 19.2 (mantiene el mismo nombre)
  arma::mat PM  = join_rows(ones(N,1), dsites);                  // termino polinomial (Matriz T) 
  arma::mat PtM = trans(join_rows(ones(M,1), ctrs));             // matriz transpuesta termino polinomial (T transpuesta)
  
  // CM de linea 16 programa 19.2
  arma::mat aux  = join_rows(CM,PM);
  arma::mat aux2 = join_rows(PtM, zeros(ncol, ncol));
  
  arma::mat A = join_cols(aux, aux2);  
  
  
  arma::vec grid = linspace(0, 1, neval);
  
  arma::mat xs(neval, neval);
  arma::mat ys(neval, neval);
  arma::mat epoints(neval*neval, neval*neval);
  
  meshgrid(xs, ys, grid, grid);
  
  epoints = join_rows(vectorise(xs), vectorise(ys));
  

  // Calcular matriz de distancia entre los puntos de evaluacion y los centros
  // linea 22 del programa 19.2
  arma::mat DM_eval = DistanceMatrix2(epoints, ctrs);
  
  // Matriz de evaluacion para datos de evaluacion
  // EM de linea 23 programa 19.2
  //arma::mat CMe = tps(ep, DM_eval); 
    
  arma::mat CMe = radialFunction(DM_eval, RBFtype, R); 
  
  // PM de linea 24 del programa 19.2
  arma::mat PMe = join_rows(ones(neval*neval, 1), epoints);
  
  // EM de linea 24 lado derecho programa 19.2 (aparece como EM = [EM PM])
  arma::mat CMe2 = join_rows(CMe, PMe);
  
  arma::colvec b = join_cols(y, zeros(3, 1));
  
  arma::colvec Pf = CMe2 * solve(A, b);
  
  // Calcular solucíon exacta, por ejemplo: evaluar 'testfunction' en los puntos de evaluacion
  // linea 26 programa 19.2
  arma::colvec exact = testfunction(epoints.col(0), epoints.col(1));
  
  // Calcular el error maximo en la grilla
  // linea 27 del programa 19.2
  double maxerr = norm(Pf - exact, "inf");
  
  return Rcpp::List::create(Rcpp::Named("phi_coefficients") = CM,
                            Rcpp::Named("Pf")       = Pf,
                            Rcpp::Named("maxerr")  = maxerr,
                            Rcpp::Named("Distances")  = DM_data);
  
  
}







/////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List RBF_LSP(arma::mat dsites, arma::mat ctrs, int RBFtype, const double R, const int neval,
                  const int omega){
  
  int N   = dsites.n_rows;
  int dim = dsites.n_cols;
  int M   = ctrs.n_rows;
  int ncol = 1 + dim;
  
  arma::mat DM_data = DistanceMatrix2(dsites, ctrs);
  
  // CM de linea 14 de programa 19.2
  arma::mat IM_0 = radialFunction(DM_data, RBFtype, R);          // Matriz E del libro Green and Silverman
  
  arma::mat IM = IM_0 + eye(IM_0.n_rows, IM_0.n_cols) / (2*omega); 
  
  
  // PM de linea 15 programa 19.2 (mantiene el mismo nombre)
  arma::mat PM  = join_rows(ones(N,1), dsites);                  // termino polinomial (Matriz T) 
  arma::mat PtM = trans(PM);          // matriz transpuesta termino polinomial (T transpuesta)
  
  // CM de linea 16 programa 19.2
  arma::mat aux  = join_rows(IM,PM);
  arma::mat aux2 = join_rows(PtM, zeros(ncol, ncol));
  
  arma::mat A = join_cols(aux, aux2);  
  
  
  // Simulacion de una variable respuesta cins distribucion Normal
  // linea 17 de programa 19.2
  arma::colvec rhs = testfunction(dsites.col(0), dsites.col(1));
  
  // linea 18 de programa 19.2 (lo nombro como rhs2 ya que no se puede renombrar en Rcpparmadillo)
  arma::colvec rhs2 = rhs + 0.03*randn(rhs.n_rows, rhs.n_cols);
  
  // linea 19 del programa 19.2 
  // Vector de variable respuesta con fila 1 = rhs3 y fila 2 matriz 3x3 de ceros
  arma::colvec rhs3 = join_cols(rhs2,zeros(3, rhs.n_cols));
  
  
  arma::vec grid = linspace(0, 1, neval);
  
  arma::mat xs(neval, neval);
  arma::mat ys(neval, neval);
  arma::mat epoints(neval*neval, neval*neval);
  
  meshgrid(xs, ys, grid, grid);
  
  epoints = join_rows(vectorise(xs), vectorise(ys));
  
  
  // Calcular matriz de distancia entre los puntos de evaluacion y los centros
  // linea 22 del programa 19.2
  arma::mat DM_eval = DistanceMatrix2(epoints, ctrs);
  
  // Matriz de evaluacion para datos de evaluacion
  // EM de linea 23 programa 19.2
  //arma::mat CMe = tps(ep, DM_eval); 
  
  arma::mat IMe = radialFunction(DM_eval, RBFtype, R); 
  
  // PM de linea 24 del programa 19.2
  arma::mat PMe = join_rows(ones(neval*neval, 1), epoints);
  
  // EM de linea 24 lado derecho programa 19.2 (aparece como EM = [EM PM])
  arma::mat IMe2 = join_rows(IMe, PMe);
  
  arma::colvec Pf = IMe2 * solve(A, rhs3);
  
  // Calcular solucíon exacta, por ejemplo: evaluar 'testfunction' en los puntos de evaluacion
  // linea 26 programa 19.2
  arma::colvec exact = testfunction(epoints.col(0), epoints.col(1));
  
  // Calcular el error maximo en la grilla
  // linea 27 del programa 19.2
  double maxerr = norm(Pf - exact, "inf");
  
  return Rcpp::List::create(Rcpp::Named("phi_coefficients") = IM,
                            Rcpp::Named("Pf")       = Pf,
                            Rcpp::Named("maxerr")  = maxerr,
                            Rcpp::Named("epoints")  = epoints,
                            Rcpp::Named("Distances")  = DM_data);

}





/////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List RBF_LSP2(arma::vec y, arma::mat dsites, arma::mat ctrs, int RBFtype, const double R, const int neval,
                  const int omega){
  
  int N   = dsites.n_rows;
  int dim = dsites.n_cols;
  int M   = ctrs.n_rows;
  int ncol = 1 + dim;
  
  arma::mat DM_data = DistanceMatrix2(dsites, ctrs);
  
  // CM de linea 14 de programa 19.2
  arma::mat IM_0 = radialFunction(DM_data, RBFtype, R);          // Matriz E del libro Green and Silverman
  
  //arma::mat IM_0 = tps(1, DM_data);          // Matriz E del libro Green and Silverman
  
  
  arma::mat IM = IM_0 + eye(IM_0.n_rows, IM_0.n_cols) / (2*omega); 
  
  
  // PM de linea 15 programa 19.2 (mantiene el mismo nombre)
  arma::mat PM  = join_rows(ones(N,1), dsites);                  // termino polinomial (Matriz T) 
  arma::mat PtM = trans(PM);          // matriz transpuesta termino polinomial (T transpuesta)
  
  // CM de linea 16 programa 19.2
  arma::mat aux  = join_rows(IM,PM);
  arma::mat aux2 = join_rows(PtM, zeros(ncol, ncol));
  
  arma::mat A = join_cols(aux, aux2);  
  
  
  arma::vec grid = linspace(0, 1, neval);
  
  arma::mat xs(neval, neval);
  arma::mat ys(neval, neval);
  arma::mat epoints(neval*neval, neval*neval);
  
  meshgrid(xs, ys, grid, grid);
  
  epoints = join_rows(vectorise(xs), vectorise(ys));
  
  
  // Calcular matriz de distancia entre los puntos de evaluacion y los centros
  // linea 22 del programa 19.2
  arma::mat DM_eval = DistanceMatrix2(epoints, ctrs);
  
  // Matriz de evaluacion para datos de evaluacion
  // EM de linea 23 programa 19.2
  //arma::mat CMe = tps(ep, DM_eval); 
  
  arma::mat IMe = radialFunction(DM_eval, RBFtype, R); 
  
  //arma::mat IMe = tps(1, DM_eval); 
  
  // PM de linea 24 del programa 19.2
  arma::mat PMe = join_rows(ones(neval*neval, 1), epoints);
  
  // EM de linea 24 lado derecho programa 19.2 (aparece como EM = [EM PM])
  arma::mat IMe2 = join_rows(IMe, PMe);
  
  arma::colvec b = join_cols(y, zeros(3, 1));
  
  arma::colvec Pf = IMe2 * solve(A, b);
  
  // Calcular solucíon exacta, por ejemplo: evaluar 'testfunction' en los puntos de evaluacion
  // linea 26 programa 19.2
  arma::colvec exact = testfunction(epoints.col(0), epoints.col(1));
  
  // Calcular el error maximo en la grilla
  // linea 27 del programa 19.2
  double maxerr = norm(Pf - exact, "inf");
  
  return Rcpp::List::create(Rcpp::Named("coefficients") = IM,
                            Rcpp::Named("Pf")       = Pf,
                            Rcpp::Named("maxerr")  = maxerr,
                            Rcpp::Named("epoints")  = epoints,
                            Rcpp::Named("Distances")  = DM_data);
                            
  
}


Rcpp::cppFunction('double logDiffExp(double x, double y)
{return x > y ? 
Rf_logspace_sub(x, y) :
Rf_logspace_sub(y, x);}') 
##
log1m <- function(x) {
  return(log1p(-x))
}
##
log1m_exp <- function(a) {
  logDiffExp(0, a)
}
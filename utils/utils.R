library(Rmpfr)

logsumexp <- function(a, prec = 64L) {

  if (!length(a)) {
    return(Rmpfr::mpfr(-Inf, prec))
  }

  # a is expected to be an mpfr vector
  a_max <- max(a)

  if (is.infinite(a_max) && a_max < 0) {
    return(Rmpfr::mpfr(-Inf, prec))
  }

  a_max + log(sum(exp(a - a_max)))
}


fma <- function(a, b, c, prec = 64) {
  a <- mpfr(a, prec)
  b <- mpfr(b, prec)
  c <- mpfr(c, prec)

  a * b + c
}


logdiffexp <- function(a, b, prec = 64) {
  a <- mpfr(a, prec)
  b <- mpfr(b, prec)

  if (b > a) {
    tmp <- a
    a <- b
    b <- tmp
  }

  if (is.infinite(a) && is.infinite(b)) {
    return(mpfr(-Inf, prec))
  }

  if (a == b) {
    return(mpfr(-Inf, prec))
  }

  a + log1p(-exp(b - a))
}
log_diff_exp <- function(x, y){
  M <- max(x, y)
  m <- min(x, y)
  A <- M
  B <- m-M
  return(
    A + Rmpfr::log1mexp(-B)
  )
}
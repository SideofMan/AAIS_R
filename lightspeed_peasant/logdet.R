logdet <- function(A){
  A = as.matrix(A)
  U = chol(A)
  y = 2*sum(log(diag(U)))
  return(y)
}
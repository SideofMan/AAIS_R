cholproj <- function(A){
  # CHOLPROJ  Projected Cholesky factorization.
  # cholproj(A) returns an upper triangular matrix U so that U'*U = A,
  # provided A is symmetric positive semidefinite (sps).
  #
  # If A is not sps, then U will approximately satisfy U'*U = A.   
  # This is useful when dealing with matrices that are affected
  # by roundoff errors.  By multiplying U'*U you effectively round A to the 
  # nearest sps matrix.
  #
  # [U,isPosDef] = cholproj(A) also returns whether A is positive definite.
  
  A = as.matrix(A)

  U = matrix(0, dim(A)[1], dim(A)[2])
  isPosDef = T
  for (i in 1:ncol(A)){
    for (j in 1:nrow(A)){
      k = 1:(i-1)
      s = A[i,j] - t(U[k,i])%*%U[k,j]
      if (i == j){
        if (s <= 0){
          isPosDef = F
          U[i,i] = 0
        } else {
          U[i,i] = sqrt(s)
        }
      } else {
        if (U[i,i] > 0){
          U[i,j] = s/U[i,i]
        } else {
          U[i,j] = 0
        }
      }
    }
  }
  return(U)
}
inv_posdef <- function(A){
  # INV_POSDEF        Invert positive definite matrix.
  # INV_POSDEF(A) is like INV(A) but faster and more numerically stable.
  # See test_inv_posdef for a timing test.
  
  # Written by Tom Minka
  
  U = cholproj(A)
  iU = inv_triu(U)
  x = iU%*%t(iU)
  return(x)
}

inv_triu <- function(U){
  # INV_TRIU     Invert upper triangular matrix.
  
  # Singularity test: 
  # inv_triu([1 1; 0 0])
  
  x = solve_triu(U,diag(dim(U)[1]))
  return(x)
}
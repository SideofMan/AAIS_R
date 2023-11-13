solve_triu <- function(T,b){
  # SOLVE_TRIU      Left division by upper triangular matrix.
  # SOLVE_TRIU(T,b) is the same as T\b but requires T to be upper triangular 
  # and runs faster.
  
  # note that this is Josh's bad version that isn't correct but should give the
  # sameish result
  
  return(solve(T)%*%b)
}
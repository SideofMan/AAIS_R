randnorm <- function(n,m,S,V){
  # RANDNORM      Sample from multivariate normal.
  # RANDNORM(n,m) returns a matrix of n columns where each column is a sample 
  # from a multivariate normal with mean m (a column vector) and unit variance.
  # RANDNORM(n,m,S) specifies the standard deviation, or more generally an 
  # upper triangular Cholesky factor of the covariance matrix.  
  # This is the most efficient option.
  # RANDNORM(n,m,[],V) specifies the covariance matrix.
  #
  # Example:
  #   x = randnorm(5, matrix(0,3,1), [], diag(3));

  `[`<-base::`[`
  
  if(nargs()==1){
    x=matrix(rnorm(n), nrow=1)
    return(x)
  }
  d=dim(m)[1]; nm=dim(m)[2]
  x=matrix(rnorm(d*n),nrow=d,ncol=n)
  if(nargs()>2){
    if(nargs()==4){
      if(d==1){
        S=sqrt(V)
      } else{
        S=chol(V)
      }
    }
    if(d==1){
      x=S*as.numeric(x)
    } else{
      x=t(S)%*%x
    }
  }
  if(nm==1){
    # x=x+repmat(m,1,n)
    x=x+matrix(m,1*nrow(m),n*ncol(m),byrow = T) # this might be much faster
  } else{
    x=x+m
  }
}
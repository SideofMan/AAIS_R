sqdist <- function(p,q,A){
  # SQDIST      Squared Euclidean or Mahalanobis distance.
  # SQDIST(p,q)   returns m(i,j) = (p(:,i) - q(:,j))'*(p(:,i) - q(:,j)).
  # SQDIST(p,q,A) returns m(i,j) = (p(:,i) - q(:,j))'*A*(p(:,i) - q(:,j)).
  
  # Written by Tom Minka
  
  d = dim(p)[1]; pn = dim(p)[2]
  d = dim(q)[1]; qn = dim(q)[2]
  
  if (pn == 0 | qn == 0){
    m = matrix(0, pn, qn)
    return(m)
  }
  
  if (nargs() == 2){
    pmag = matrix(colSums(p*p), nrow=1)
    qmag = matrix(colSums(q*q), nrow=1)
    m = repmat(qmag, pn, 1) + repmat(t(pmag), 1, qn) - 2*t(p)%*%q
  } else {
    Ap = A%*%p
    Aq = A%*%q
    pmag = matrix(colSums(p*Ap), nrow=1)
    qmag = matrix(colSums(q*Aq), nrow=1)
    m = repmat(qmag, pn, 1) + repmat(t(pmag), 1, qn) - 2*t(p)%*%Aq
  }
  
 return(m)
}
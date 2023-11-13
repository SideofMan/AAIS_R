library(hypergeo)
library(BayesVarSel)
library(pracma)

X = as.matrix(unname(Hald[,1:4]))
n=dim(X)[1]
X=cbind(ones(n,1),X)
n=dim(X)[1]; k=dim(X)[2]
Y = as.matrix(unname(Hald[,-(1:4)]))

testX = X[1:12,]
testX
testY = Y[1:12,]
n=12;

# X2 = X[,-1]
# BigSigma22 = solve(t(X2)%*%(diag(n) - 1/n*ones(n,n))%*%X2)

LinearModelML_exact <- function(X,Y){
  S2=var(Y)*(n-1)
  Q=1/S2*t(Y)%*%(diag(n) - X%*%solve(t(X)%*%X)%*%t(X))%*%Y; Q=as.numeric(Q)
  
  output = gamma((n-1)/2)/(k*sqrt(n)*(pi*S2*Q)^((n-1)/2))*(k/(n+1))^((k-1)/2)*hypergeo(k/2,(n-1)/2,(k+2)/2,k*(1-1/Q)/(n+1))
  
  return(as.numeric(Re(output)))
}

# LinearModelML_exact(X,Y)
LinearModelML_exact(testX,testY)

# run through all possible combinations of covariates and the one with the biggest ML is the "best"

library(hypergeo)
library(BayesVarSel)
library(pracma)
library(rje)

X = as.matrix(unname(Hald[,1:4]))
n=dim(X)[1]
X=cbind(ones(n,1),X)
n=dim(X)[1]; k=dim(X)[2]
Y = as.matrix(unname(Hald[,-(1:4)]))

testX = X[,]
testX
testY = Y[,]
n=dim(testX)[1];

LinearModelML_exact <- function(X,Y){
  S2=var(Y)*(n-1)
  Q=1/S2*t(Y)%*%(diag(n) - X%*%solve(t(X)%*%X)%*%t(X))%*%Y; Q=as.numeric(Q)
  
  output = gamma((n-1)/2)/(k*sqrt(n)*(pi*S2*Q)^((n-1)/2))*(k/(n+1))^((k-1)/2)*hypergeo(k/2,(n-1)/2,(k+2)/2,k*(1-1/Q)/(n+1))
  
  return(as.numeric(Re(output)))
}

LinearModelML_exact(testX,testY)

# run through all possible combinations of covariates and the one with the biggest ML is the "best"
data=matrix(Y,ncol=1)
LinearModel_likelihood()

#---------------just keeping track of MarLik vs exact--------------#
my_df <- read.csv("C:/Users/jdseidma/Dropbox/Research/SU23/AAIS/GitHub/AAIS_R/linear_model_ML.csv")

df <- my_df
# df = data.frame(matrix(ncol = 4, nrow = 2^4-1)); colnames(df) <- c("algorithm_run","MarLik","Exact_MarLik","Ratio")
# columns <- 1:(2^4-1)
# df$algorithm_run<-columns
# columns_in_run = powerSet(1:4)

#-------add a new data frame row here--------#
# columns_in_run = powerSet(1:4)
#---which columns of Hald did you run?---#
# your_columns = c(1,2,3)
# ROW = which(unlist(lapply(columns_in_run,function(e) identical(as.numeric(e),c(1,2,3)))))
# df[ROW,2]=MarLik; df[ROW,3]=LinearModelML_exact(linearX,data); df[ROW,4]=df[ROW,3]/df[ROW,]
# df[15,2]=6.670089018839259923516e-19; df[15,3]=3.822703833994052497434e-18; df[15,4]=df[15,2]/df[15,3]

write.csv(df, "C:/Users/jdseidma/Dropbox/Research/SU23/AAIS/GitHub/AAIS_R/linear_model_ML.csv", row.names=FALSE)


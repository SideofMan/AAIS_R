library(hypergeo)
library(BayesVarSel)
library(pracma)
library(rje)

#---------Source this script to use readline() properly-------#

# The Hald data cleaned up ------------------------------------------------

X = as.matrix(unname(Hald[,1:4]))
n=dim(X)[1]
X=cbind(ones(n,1),X)
n=dim(X)[1]; k=dim(X)[2]
Y = as.matrix(unname(Hald[,5]))

LinearModelML_exact(X,Y)

# run through all possible combinations of covariates and the one with the biggest ML is the "best"


# Only run this section if you need to restart the data frame from --------
# 
df = data.frame(matrix(ncol = 6, nrow = 2^4-1)); colnames(df) <- c("algorithm_run","M","ESS","MarLik","Exact_MarLik","Ratio")
columns <- 1:(2^4-1)
df$algorithm_run<-columns

# Read the data frame of stored AAIS runs ---------------------------------

my_df <- read.csv("C:/Users/jdseidma/Dropbox/Research/SU23/AAIS/GitHub/AAIS_R/linear_model_ML_long_run.csv")
df <- my_df

# add a new run of the AAIS to the data frame
columns_in_run = powerSet(1:4)[-1] # all possible linear model runs
# your_columns = readline(prompt = "Enter the vector of columns you ran AAIS with: ") # the linear model run you did
# your_columns = eval(parse(text=your_columns))

your_columns = columns # the linear model run you did
ROW = which(unlist(lapply(columns_in_run,function(e) identical(as.numeric(e),as.numeric(your_columns))))) # row in the df for your run

df[ROW,2]=ESS_final; df[ROW,3]=MarLik; df[ROW,4]=LinearModelML_exact(linearX,data); df[ROW,5]=df[ROW,3]/df[ROW,4]

print(df)

write.csv(df, "C:/Users/jdseidma/Dropbox/Research/SU23/AAIS/GitHub/AAIS_R/linear_model_ML_sigma1e2.csv", row.names=FALSE)


#install.packages("hypergeo")
#install.packages("BayesVarSel")
#install.packages("mvtnorm")
#install.packages("abind")
#install.packages("foreach")
#install.packages("doParallel")
#install.packages("abind")
#install.packages("pracma")
#install.packages("NPflow") # not necessary to include in final script

library("hypergeo")
library("BayesVarSel")
library("mvtnorm")
library("abind")
library("profvis")
library("foreach")
library("doParallel")
library("abind")
library("pracma")
library("ggplot2")
library("plotly")
tic()
rm(list = ls())
options(error=browser)

#-----------alter some base functions-----------------------------#
# we need to alter some base functions in R to make sure that the algorithm
# works for all dimensions

`[` <- function(...) base::`[`(...,drop=FALSE)
drop <- function(S){
  if(all(dim(S)==c(1,1,1))){  # ensures that the algorithm works if your target function is 1D
    return(matrix(S)) 
  }else{
    return(base::drop(S))
  }
}

#----------source the function scripts needed here----------------#

# setwd("C:/Users/jdseidma/Dropbox/Research/SU23/AAIS/Seidman R code/lightspeed_peasant") # sourcing lightspeed scripts
# files.sources = list.files()
# sapply(files.sources, source)
#
# setwd('..')
# files.sources = list.files(pattern = '*.R') # sourcing helper functions
# files.sources = files.sources[grepl("\\.R$", files.sources)]
# files.sources = files.sources[grep("main", files.sources, invert = T)]
# files.sources = files.sources[grep("Plot", files.sources, invert = T)]
# sapply(files.sources, source)

setwd("C:/Users/jdseidma/Dropbox/Research/SU23/AAIS/GitHub/AAIS_R") # sourcing all scripts
source("all_functions.R")

#-------------set up/initialize empty lists------------#

X=list()
X2=list()
X_try=list()
Proposal=list()
Proposal_temp=list()
Proposal_try=list()
ESS=c()

#------Define the name of your Target Function here--------#

Target_function='LinearModelHald'

#-----Simulating Annealing Importance Sampling-------------#

if (Target_function=='outerproduct'){
  data=1
  dim=7 # dimension of the target likelihood function
  N=2e4;X2=list();X2$N=2e4 # N: particle size of importance sampling; X2$N: sample size for adding a new mixture component
  Proposal=list();Proposal$M=10; # Number of mixing components in the initial proposal
  df=5  # degree of freedom of student's t distributions involved in all proposals
  # Proposal$W=rep(1/Proposal$M, Proposal$M) # Initial weight of each mixing component
  Proposal$W=matrix(1/Proposal$M,1,Proposal$M) # Initial weight of each mixing component
  Proposal$Mu=matrix(runif(Proposal$M*dim, min = -10, max = 10), nrow = Proposal$M, ncol = dim) # Initial mean of mixing components
  Proposal$Sigma=array(0, dim = c(dim, dim, Proposal$M)) # Initial covariance of mixing components
  
  Sigma_initial=1e3*diag(dim) # Initial covariance of added components
  for (i in 1:Proposal$M){
    Proposal$Sigma[,,i]=1e3*diag(dim); # Initial covariance of mixing components
  }
  
  gama=c(0.0001, 0.001, 0.01, seq(0.05,1,by=0.05)) # the annealing schedule
  print('outer product toy function')
  
  Proposal0=Proposal
  ESS_Arr=c()
  
  Cor_thr=0.8
  ESS_thr=0.6
  str0='-----Made by Bin Liu-----'
  w_thr=1e-5
  
  # tracking_length corresponds to how long it will take to break the add component while loop
  # I'd suggest something like 10 (or 20 if you'd like it to run longer)
  tracking_length=10
  tracking_thr=0.1 # something like 0.1 for faster runtimes, 0.001 for slower but more added components
  
  # determine the number of components you would like to add when adding components: 1 or 2
  # please note that adding 2 components is faster but may lead to poor and approximations by the algorithm
  add_component_number=2
} else if(Target_function=='Rastrigin'){
  data=1
  dim=2 # dimension of the target likelihood function
  N=4e4;X2=list();X2$N=4e3 # N: particle size of importance sampling; X2$N: sample size for adding a new mixture component
  Proposal=list();Proposal$M=1; # Number of mixing components in the initial proposal
  df=300  # degree of freedom of student's t distributions involved in all proposals
  # Proposal$W=rep(1/Proposal$M, Proposal$M) # Initial weight of each mixing component
  Proposal$W=matrix(1/Proposal$M,1,Proposal$M) # Initial weight of each mixing component
  Proposal$Mu=zeros(Proposal$M,dim) # Initial mean of mixing components
  Proposal$Sigma=array(0, dim = c(dim, dim, Proposal$M)) # Initial covariance of mixing components
  
  Sigma_initial=1.6e-1*diag(dim) # changed from 1e1 to 1e0 for cutoff
  for (i in 1:Proposal$M){
    Proposal$Sigma[,,i]=2e0*diag(dim);
  }
  
  gama=c(0.001, 0.01, 0.1, 1) # the annealing schedule
  print('Rastrigin 2d toy function')
  
  Proposal0=Proposal
  ESS_Arr=c()
  
  Cor_thr=0.9
  ESS_thr=0.9
  str0='-----Made by Bin Liu-----'
  plot_id=F
  axis_lmt=c(-10,10,-10,10) # FLAG
  w_thr=1e-6 # If the weights of a mixture component is smaller than this value, delete it
  
  # tracking_length corresponds to how long it will take to break the add component while loop
  # I'd suggest something like 10 (or 20 if you'd like it to run longer)
  tracking_length=20
  tracking_thr=0.001 # something like 0.1 for faster runtimes, 0.001 for slower but more added components
  
  # determine the number of components you would like to add when adding components: 1 or 2
  # please note that adding 2 components is faster but may lead to poor and approximations by the algorithm
  add_component_number=1
}else if(Target_function=='helix'){
  data=1
  dim=3 # dimension of the target likelihood function
  N=2e4;X2=list();X2$N=2e3 # N: particle size of importance sampling; X2$N: sample size for adding a new mixture component
  Proposal=list();Proposal$M=10; # Number of mixing components in the initial proposal
  df=5  # degree of freedom of student's t distributions involved in all proposals
  # Proposal$W=rep(1/Proposal$M, Proposal$M) # Initial weight of each mixing component
  Proposal$W=matrix(1/Proposal$M,1,Proposal$M) # Initial weight of each mixing component
  # Proposal$Mu=matrix(runif(Proposal$M*dim, min = -10, max = 10), nrow = Proposal$M, ncol = dim) # Initial mean of mixing components
  Proposal$Mu=cbind(matrix(runif(Proposal$M*2, min = -50, max = 50), nrow = Proposal$M, ncol = 2), matrix(runif(Proposal$M*1, min = 0, max = 40), nrow = Proposal$M, ncol = 1))
  Proposal$Sigma=array(0, dim = c(dim, dim, Proposal$M)) # Initial covariance of mixing components
  
  Sigma_initial=1e3*diag(dim)
  for (i in 1:Proposal$M){
    Proposal$Sigma[,,i]=1e3*diag(dim);
  }
  
  gama=c(0.0001, 0.001, 0.01, seq(0.05,1,by=0.05)) # the annealing schedule
  gama=seq(0.1,1,by=0.1)
  print('Helix 3d toy function')
  
  Proposal0=Proposal
  ESS_Arr=c()
  
  Cor_thr=0.7
  ESS_thr=0.9
  str0='-----Made by Bin Liu-----'
  plot_id=F
  axis_lmt=c(-100,100,-100,100) # FLAG
  w_thr=1e-4
  
  # tracking_length corresponds to how long it will take to break the add component while loop
  # I'd suggest something like 10 (or 20 if you'd like it to run longer)
  tracking_length=20 # change back to 20 and 0.001
  tracking_thr=0.001 # something like 0.1 for faster runtimes, 0.001 for slower but more added components
  
  # determine the number of components you would like to add when adding components: 1 or 2
  # please note that adding 2 components is faster but may lead to poor approximations by the algorithm in some cases
  add_component_number=2
}else if(Target_function == "LinearModelHald"){
  # specify the number of columns to test (between 1 to 4 for Hald)
  columns=c(4)
  X = as.matrix(unname(Hald[,columns]))
  n=dim(X)[1]
  X=cbind(ones(n,1),X)
  n=dim(X)[1]; k=dim(X)[2]
  Y = as.matrix(unname(Hald[,5]))
  
  X2 = X[,-1]
  BigSigma22 <<- solve(t(X2)%*%(diag(n) - 1/n*ones(n,n))%*%X2)
  linearX <<- X
  
  data=matrix(Y,ncol = 1)
  dim=length(columns)+3 # dimension of the target likelihood function, this will change based on columns above (#columns + 3)
  N=2e4;X2=list();X2$N=2e4 # N: particle size of importance sampling; X2$N: sample size for adding a new mixture component
  Proposal=list();Proposal$M=10; # Number of mixing components in the initial proposal
  df=300  # degree of freedom of student's t distributions involved in all proposals
  # Proposal$W=rep(1/Proposal$M, Proposal$M) # Initial weight of each mixing component
  Proposal$W=matrix(1/Proposal$M,1,Proposal$M) # Initial weight of each mixing component
  Proposal$Mu=matrix(runif(Proposal$M*dim, min = -10, max = 10), nrow = Proposal$M, ncol = dim) # Initial mean of mixing components
  Proposal$Sigma=array(0, dim = c(dim, dim, Proposal$M)) # Initial covariance of mixing components
  
  Sigma_initial=1e2*diag(dim) # Initial covariance of added components
  for (i in 1:Proposal$M){
    Proposal$Sigma[,,i]=1e2*diag(dim); # Initial covariance of mixing components
  }
  
  gama=c(0.0001, 0.001, 0.01, seq(0.05,1,by=0.05)) # the annealing schedule
  gama=c(0.001, 0.01, 0.1, seq(0.5,1,by=.1)) # the annealing schedule
  print('Linear Model with Hald dataset')
  
  Proposal0=Proposal
  ESS_Arr=c()
  
  Cor_thr=0.8
  ESS_thr=0.6
  str0='-----Made by Bin Liu-----'
  w_thr=1e-4
  
  # tracking_length corresponds to how long it will take to break the add component while loop
  # I'd suggest something like 10 (or 20 if you'd like it to run longer)
  tracking_length=10
  tracking_thr=0.1 # something like 0.1 for faster runtimes, 0.001 for slower but more added components
  
  # determine the number of components you would like to add when adding components: 1 or 2
  # please note that adding 2 components is faster but may lead to poor and approximations by the algorithm
  add_component_number=1
}

for (j in 1:length(gama)){
  X$N=N
  X$Values=matrix(0,X$N,dim)
  X$Root=matrix(0,X$N,1)
  
  #---------one iteration of IS with EM------------#
  if (j==1){
    if (Target_function=='outerproduct'){
      X$Values=matrix(runif(X$N*dim, min = -10, max = 10), nrow = X$N, ncol = dim) # initial random samples
    }else if (Target_function=='Rastrigin'){
      X$Values=matrix(runif(X$N*dim, min = -5.12, max = 5.12), nrow = X$N, ncol = dim) # initial random samples
    }else if (Target_function=='helix'){
      # X$Values=matrix(runif(X$N*dim, min = -10, max = 10), nrow = X$N, ncol = dim) # initial random samples
      X$Values=cbind(matrix(runif(X$N*2, min = -50, max = 50), nrow = X$N, ncol = 2), matrix(runif(X$N*1, min = 0, max = 40), nrow = X$N, ncol = 1)) # initial random samples # initial random samples
    }else if (Target_function=='LinearModelHald'){
      X$Values=matrix(runif(X$N*dim, min = -10, max = 10), nrow = X$N, ncol = dim) # initial random samples
    }
    
    X$Resp=matrix(0,X$N,Proposal$M) # responsibility of each component with regard to each particle
    for (i in 1:Proposal$M){
      X$Resp[,i]=exp(log_t_pdf(X$Values,Proposal$Mu[i,],drop(Proposal$Sigma[,,i]),df)) # f in eqn 6
    }
    X$Proposal=X$Resp%*%t(Proposal$W) # eqn 6
    X$logProposal=log(X$Proposal)
    X$logProposal0=X$logProposal
    
    output = do.call(Target_function, args=list(X$Values,data))
    X$logPrior=output[[1]]; X$logLike=output[[2]] # logs of prior and likelihood, p in the paper
    X$logTarget=X$logPrior+X$logLike # p in the paper
    X$logAnnealTarget=X$logTarget*gama[j]+X$logProposal0*(1-gama[j]) # log of eqn 7
    
    X$logWeight=X$logAnnealTarget-X$logProposal # log of eqn 8
    logWeight_scaled=X$logWeight-max(X$logWeight)+10 # scaling the logs of the weights
    weight_temp=exp(logWeight_scaled)
    X$NormalizedWeight=weight_temp/sum(weight_temp) # eqn 9
    
    Proposal=t_mix_update_v2(X,Proposal,df) # initial EM algorithm to adjust proposal model, eqns 20-23
    Proposal$W=Proposal$W/sum(Proposal$W)
    X=t_mix_sample(Proposal,X,df) # sample from updated proposal model
    
    X$Resp=matrix(0,X$N,Proposal$M)
    for (i in 1:Proposal$M){
      X$Resp[,i]=exp(log_t_pdf(X$Values,Proposal$Mu[i,],drop(Proposal$Sigma[,,i]),df)) # f in eqn 6, updated proposal
    }
    X$Proposal=X$Resp%*%t(Proposal$W) # eqn 6
    X$logProposal=log(X$Proposal)
    
    rp_temp=matrix(0,X$N,Proposal0$M)
    for (i in 1:Proposal0$M){
      rp_temp[,i]=exp(log_t_pdf(X$Values,Proposal0$Mu[i,],drop(Proposal0$Sigma[,,i]),df)) # f in eqn 6, original proposal
    }
    X$logProposal0=log(rp_temp%*%t(Proposal0$W)) # log of q in paper
    
    output = do.call(Target_function, args=list(X$Values,data))
    X$logPrior=output[[1]]; X$logLike=output[[2]] # logs of prior and likelihood, p in the paper
    X$logTarget=X$logPrior+X$logLike # p in the paper
    X$logAnnealTarget=X$logTarget*gama[j]+X$logProposal0*(1-gama[j]) # log of eqn 7
    
    X$logWeight=X$logAnnealTarget-X$logProposal # log of eqn 8
    logWeight_scaled=X$logWeight-max(X$logWeight)+10  # scaling the logs of the weights
    weight_temp=exp(logWeight_scaled)
    X$NormalizedWeight=weight_temp/sum(weight_temp) # eqn 9
  } else {
    output=ISEM_astro(X,Proposal,Proposal0,df,Target_function,gama,j,data,w_thr) # one iteration of importance sampling and EM, section 3, step 3d
    X=output[[1]]; Proposal=output[[2]]
  }
  #------------------------------------------------#
  
  ESS[j]=1/sum(X$NormalizedWeight^2)/X$N
  string <- sprintf('ESS after EM: %.4f', ESS[j]); print(string)
  
  str1=as.character(j)
  str2='th iter--After EM'
  strf=paste0(str1, str2, str0, sep="")
  
  output=Merge(Proposal,X,Cor_thr,df,dim) # section 3, step 3b
  Proposal=output[[1]]; X=output[[2]]; Sign=output[[3]]
  if (Sign==T){ # if components were merged
    ESS[j]=1/sum(X$NormalizedWeight^2)/X$N
    
    str1=as.character(j)
    str2='th iter--After Merging--ESS='
    str3='%.4f'
    strf=paste0(str1, str2, str3, sep="")
    string <- sprintf(strf, ESS[j]); print(string)
    
    str2='th iter--After Merging'
    strf=paste0(str1, str2, str0, sep="")
  }
  
  #-------------add component---------------#
  j2=1
  ESS_add=ESS[j]
  while (ESS[j]<ESS_thr){ # while our ESS is below the threshold, add components (section 3.3)
    if(add_component_number==1){
      #-----initializing Proposal_temp stuff ----------------#
      Proposal_temp$M=1
      Proposal_temp$Mu=zeros(Proposal_temp$M,dim)
      Proposal_temp$Sigma=array(0, dim = c(dim,dim,Proposal_temp$M))
      rp_temp=matrix(0,X$N,Proposal_temp$M)
      
      ind_pt=which.max(X$NormalizedWeight) # index of maximum weight particle, section 3.3, step 1a
      Proposal_temp$Mu=X$Values[ind_pt,]
      Proposal_temp$Sigma[,,1]=Sigma_initial
      Proposal_temp$W=matrix(1,1,1)
      
      X2$Values=matrix(0,X2$N,dim)
      # output=ISEM_astro(X2,Proposal_temp,Proposal0,df,Target_function,gama,j,data,w_thr) # section 3.3, step 1b-c
      # X2=output[[1]]; Proposal_temp=output[[2]]
      X2=t_mix_sample(Proposal_temp,X2,df)
      #------#
      X2$Resp=zeros(X2$N,Proposal_temp$M)
      for(i in 1:Proposal_temp$M){
        X2$Resp[,i]=exp(log_t_pdf(X2$Values,Proposal_temp$Mu[i,],drop(Proposal_temp$Sigma[,,i]),df))
      }
      X2$Proposal=X2$Resp%*%t(Proposal_temp$W)
      X2$logProposal=log(X2$Proposal)
      
      rp_temp=zeros(X2$N,Proposal0$M)
      for(i in 1:Proposal0$M){
        rp_temp[,i]=exp(log_t_pdf(X2$Values,Proposal0$Mu[i,],drop(Proposal0$Sigma[,,i]),df))
      }
      X2$logProposal0=log(rp_temp%*%t(Proposal0$W))
      
      output = do.call(Target_function, args=list(X2$Values,data)) 
      X2$logPrior=output[[1]]; X2$logLike=output[[2]]
      X2$logTarget=X2$logPrior+X2$logLike
      X2$logAnnealTarget=X2$logTarget*gama[j]+X2$logProposal0*(1-gama[j])
      
      X2$logWeight=X2$logAnnealTarget-X2$logProposal
      logWeight_scaled=X2$logweight-max(X2$logWeight)+10
      weight_temp=exp(logWeight_scaled)
      X2$NormalizedWeight=weight_temp/sum(weight_temp)
      #--------#
      
      
      ESS_temp=1/sum(X2$NormalizedWeight^2)/X2$N
      Sigma_temp=Proposal$Sigma
      Sigma_temp=abind(Sigma_temp,Proposal_temp$Sigma,along=3)
      
      str1=as.character(j)
      str2='th iter--'
      str3='try to add a new component'
      strf=paste0(str1, str2, str3, str0, sep="")
      
      Proposal_try$M=Proposal$M+1
      Proposal_try$Mu=rbind(Proposal$Mu, Proposal_temp$Mu)
      Proposal_try$Sigma=Proposal$Sigma
      Proposal_try$Sigma=abind(Proposal_try$Sigma,Proposal_temp$Sigma,along=3)
      Proposal_try$W=cbind((1-X2$N/(X$N+X2$N))*Proposal$W, X2$N/(X$N+X2$N)*Proposal_temp$W) # section 3.3, step 2
      
      X_try$Values=rbind(X$Values,X2$Values)
      X_try$N=X$N+X2$N
      X_try$Resp=X$Resp
      
      X_try_Resp_temp=zeros(X2$N,Proposal$M)
      for (i in 1:Proposal$M){
        X_try_Resp_temp[,i]=exp(log_t_pdf(X2$Values,Proposal$Mu[i,],drop(Proposal$Sigma[,,i]),df)) # section 3.3, step 3 (this is q_trial)
      }
      X_try$Resp=rbind(X_try$Resp,X_try_Resp_temp)    
      
      rp_temp=exp(log_t_pdf(X$Values,Proposal_temp$Mu,drop(Proposal_temp$Sigma),df))
      X_try$Resp=cbind(X_try$Resp[1:(X$N+X2$N),], rbind(rp_temp,X2$Resp))
    }else if(add_component_number==2){
      #-----initializing Proposal_temp stuff ----------------#
      Proposal_temp$M=2 # we add 2 components instead of 1 like the paper says
      Proposal_temp$Mu=zeros(Proposal_temp$M,dim)
      Proposal_temp$Sigma=array(0, dim = c(dim,dim,Proposal_temp$M))
      rp_temp=matrix(0,X$N,Proposal_temp$M)
      ind_pt=which.max(X$NormalizedWeight) # index of maximum weight particle, section 3.3, step 1a
      Proposal_temp$Mu[1,]=X$Values[ind_pt,]
      Proposal_temp$Mu[2,]=-Proposal_temp$Mu[1,]
      Proposal_temp$Sigma[,,1]=Sigma_initial
      Proposal_temp$Sigma[,,2]=Sigma_initial
      Proposal_temp$W=matrix(0.5,1,2)
      
      X2$Values=matrix(0,X2$N,dim)
      output=ISEM_astro(X2,Proposal_temp,Proposal0,df,Target_function,gama,j,data,w_thr) # section 3.3, step 1b-c
      X2=output[[1]]; Proposal_temp=output[[2]]
      if(length(Proposal_temp$W)[1]<2){
        print(Proposal_temp$W)
        browser()
      }
      ESS_temp=1/sum(X2$NormalizedWeight^2)/X2$N
      Sigma_temp=Proposal$Sigma
      Sigma_temp=abind(Sigma_temp,Proposal_temp$Sigma,along=3)
      
      str1=as.character(j)
      str2='th iter--'
      str3='try to add a new component'
      strf=paste0(str1, str2, str3, str0, sep="")
      
      Proposal_try$M=Proposal$M+2
      Proposal_try$Mu=rbind(Proposal$Mu, Proposal_temp$Mu)
      Proposal_try$Sigma=Proposal$Sigma
      Proposal_try$Sigma=abind(Proposal_try$Sigma,Proposal_temp$Sigma,along=3)
      Proposal_try$W=cbind((1-X2$N/(X$N+X2$N))*Proposal$W, X2$N/(X$N+X2$N)*Proposal_temp$W) # section 3.3, step 2
      
      X_try$Values=rbind(X$Values,X2$Values)
      X_try$N=X$N+X2$N
      X_try$Resp=X$Resp
      
      X_try_Resp_temp=zeros(X2$N,Proposal$M)
      for (i in 1:Proposal$M){
        X_try_Resp_temp[,i]=exp(log_t_pdf(X2$Values,Proposal$Mu[i,],drop(Proposal$Sigma[,,i]),df)) # section 3.3, step 3 (this is q_trial)
      }
      X_try$Resp=rbind(X_try$Resp,X_try_Resp_temp)
      
      rp_temp[,1]=exp(log_t_pdf(X$Values,Proposal_temp$Mu[1,],drop(Proposal_temp$Sigma[,,1]),df))
      rp_temp[,2]=exp(log_t_pdf(X$Values,Proposal_temp$Mu[2,],drop(Proposal_temp$Sigma[,,2]),df))
      X_try$Resp=cbind(X_try$Resp[1:(X$N+X2$N),], rbind(rp_temp,X2$Resp))
    }
    
    index_survival=which(Proposal_try$W>w_thr) # delete components below weight threshold
    Proposal_try$M=length(index_survival)
    Proposal_try$Mu=Proposal_try$Mu[index_survival,]
    Proposal_try$Sigma=Proposal_try$Sigma[,,index_survival]
    Proposal_try$W=Proposal_try$W[,index_survival]
    Proposal_try$W=Proposal_try$W/sum(Proposal_try$W)
    
    X_try$Resp=X_try$Resp[,index_survival]
    X_try$Proposal=X_try$Resp%*%t(Proposal_try$W)
    X_try$logProposal=log(X_try$Proposal)
    X_try$logAnnealTarget=rbind(X$logAnnealTarget, X2$logAnnealTarget)
    X_try$logWeight=X_try$logAnnealTarget-X_try$logProposal
    
    logWeight_scaled=X_try$logWeight-max(X_try$logWeight)+10
    weight_temp=exp(logWeight_scaled)
    X_try$NormalizedWeight=weight_temp/sum(weight_temp)
    
    Proposal=Proposal_try
    X=X_try
    ESS[j]=1/sum(X_try$NormalizedWeight^2)/(X_try$N) # new ESS to compare to threshold (section 3.3, step 6)
    
    str2='th iter--'
    str3=as.character(j2)
    str4='th new component added:ESS='
    str5='%.4f'
    strf=paste0(str1,str2,str3,str4,str5,sep="")
    string <- sprintf(strf,ESS[j]); print(string)
    
    str2='th iter--'
    str3=as.character(j2)
    str4='th new component added'
    
    strf=paste0(str1,str2,str3,str4,str0,sep="")
    
    j2=j2+1
    ESS_add=cbind(ESS_add, ESS[j])
    
    # fail safe built in to make sure we aren't stuck in this loop forever
    # you can change the parameters tracking_length and tracking_thr at beginning of code
    if (j2%%10==0 & j2>=tracking_length){
      end=length(ESS_add)
      tmptmp=sum(ESS_add[(end-tracking_length+1):(end-as.integer(tracking_length/2))])
      if ((sum(ESS_add[(end-as.integer(tracking_length/2)+1):end])-tmptmp)/tmptmp<tracking_thr){
        break
      }
    }
  }
  
  while (T){ # does ISEM until we get a model with better ESS than the current one (not in paper)
    ESS_em=ESS[j]
    X$N=N
    X$Values=matrix(0,X$N,dim)
    X$Root=matrix(0,X$N,1)
    # Sample
    output=ISEM_astro(X,Proposal,Proposal0,df,Target_function,gama,j,data,w_thr)
    X=output[[1]];Proposal=output[[2]]
    
    ESS[j]=1/sum(X$NormalizedWeight^2)/(X$N)
    if (ESS[j]<ESS_em){
      break
    }
  }
  
  str2='th iter--'
  str3=as.character(j2)
  str4='th new component added'
  
  strf=paste0(str1,str2,str3,str4,str0,sep="")
  
  ESS_Arr=cbind(ESS_Arr, ESS[j])
  string <- sprintf(paste0('SAIS_final ', Target_function, ':j/Temp/ESS/M %.0f / %1.4f / %.4f / %.0f'), j,gama[j],ESS[j],Proposal$M); print(string)
}
toc()

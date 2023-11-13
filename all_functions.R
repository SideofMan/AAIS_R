# this is all helper functions in one document

#--------------lightspeed helper functions-----------------#

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

inv_triu <- function(U){
  # INV_TRIU     Invert upper triangular matrix.

  # Singularity test:
  # inv_triu([1 1; 0 0])

  x = solve_triu(U,diag(dim(U)[1]))
  return(x)
}

logdet <- function(A){
  A = as.matrix(A)
  U = chol(A)
  y = 2*sum(log(diag(U)))
  return(y)
}

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

solve_triu <- function(Tri_matrix,b){
  # SOLVE_TRIU      Left division by upper triangular matrix.
  # SOLVE_TRIU(T,b) is the same as T\b but requires T to be upper triangular
  # and runs faster.

  # note that this is Josh's bad version that isn't correct but should give the
  # sameish result

  return(solve(Tri_matrix)%*%b)
}

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

  # return(m)
  return(abs(m))
}

#------------AAIS helper functions----------------#

dsn <- function(x,location,scale,shape){
  #dsn psn qsn rsn
  #Skew-Normal Distribution
  #
  #DESCRIPTION
  #
  #Density function, distribution function, quantiles and random number
  #generation for the skew-normal (SN) distribution.
  #
  #USAGE
  #
  #dsn(x, location, scale, shape)
  #psn(q, location, scale, shape)
  #qsn(p, location, scale, shape, tol)
  #rsn(n, location, scale, shape)
  #
  #REQUIRED ARGUMENTS
  #
  #x	vector of quantiles. Missing values (NaN) are allowed.
  #q	vector of quantiles. Missing values (NaN) are allowed.
  #p	vector of probabilities. Missing values (NaN) are allowed.
  #
  #OPTIONAL ARGUMENTS
  #
  #location vector of location parameters (default is 0).
  #scale	  vector of (positive) scale parameters (default is 1).
  #shape	  vector of shape parameters. With 'psn' and 'qsn', it must
  #	  be of length 1 (default is 0).
  #n	  sample size (default is 1).
  #tol	  a scal value which regulates the accuracy of the result.
  #         (default is 1e-8)
  #
  #VALUE
  #
  #density (dsn), probability (psn), quantile (qsn) or  random sample (rsn)
  #from the skew-normal distribution with given location, scale and shape
  #parameters.
  #
  #BACKGROUND
  #
  #The family of skew-normal distributions is an extension of the normal
  #family, via the introdution of a shape parameter which regulates
  #skewness; when shape=0, the skew-normal distribution reduces to the
  #normal one.  The density of the SN distribution when location=0 and
  #scale=1 is 2*dnorm(x)*pnorm(shape*x). A multivariate version of the
  #distribution exists. See the references below for additional
  #information.
  #
  #DETAILS
  #
  #psn make use of function T_Owen
  #
  #REFERENCES
  #
  #Azzalini, A. (1985). A class of distributions which includes the normal
  #ones. Scand. J. Statist. 12, 171-178.
  #
  #Azzalini, A. and Dalla Valle, A. (1996). The multivariate skew-normal
  #distribution. Biometrika 83, 715-726.
  #
  #SEE ALSO
  #
  #T_Owen, dmsn, sn_mle
  #
  #EXAMPLES
  #
  #pdf = dsn(seq(-3,3,0.1),0,1,3)
  #cdf = psn(seq(-3,3,0.1),0,1,3)
  #qu = qsn(seq(0.1,0.9,0.1),0,1,-2)
  #rn = rsn(100, 5, 2, 5)

  nargin = length(match.call()) - 1

  if (nargin<4){
    shape = 0
  }
  if (nargin<3){
    scale = 1
  }
  if (nargin<2){
    location=0
  }
  if (nargin<1){
    stop('Argument x is missing');
  }
  if (is.nan(location)){
    location=0
  }
  if (is.nan(scale)){
    scale=1
  }
  if (is.nan(shape)){
    shape=0
  }
  d=2*dnorm( (x-location)/scale )*pnorm((shape*(x-location)/scale) )/scale

  return(d)
}

helix <- function(X,data){
  return(helix_correct(X,data))
  loglike=0*data

  # Target distribution being the flared Helix
  t=linspace(0,6*pi,1e4)
  xt=t*cos(t)
  yt=t*sin(t)
  zt=2*t
  mut=rbind(xt,yt,zt)
  S=diag(3)

  # n=dim(X)[1]; p=dim(X)[2]
  # pdf=zeros(n,1)
  # for(k in 1:n){
  #   pdf[k,1]=log(max(exp(-colSums((repmat(rbind(X[k,1], X[k,2], X[k,3]),1,dim(mut)[2])-mut)^2)/0.2)*(2*pi)^(-3/2)*3^(-1/2)))
  # }
  # logprior=pdf

  # faster way of doing it by 3x
  n=dim(X)[1]; p=dim(X)[2]
  pdf=zeros(n,1)
  # pdf=matrix(apply(X,MARGIN = 1,FUN = function(x) log(max(exp(-colSums((mut-x)^2)/0.2)*(2*pi)^(-3/2)*3^(-1/2)))))
  pdf=matrix(apply(X,MARGIN = 1,FUN = function(x) log(exp(-min(colSums((mut-x)^2))/0.2)*(2*pi)^(-3/2)*3^(-1/2))))
  logprior=pdf

  output=list(logprior,loglike)
  return(output)
}

helix_correct <- function(X,data){
  # Target distribution is the 3D flared helix from the paper
  loglike=0*data
  
  # n=dim(X)[1]; p=dim(X)[2]
  # pdf=zeros(n,1)
  # for(k in 1:n){
  #   if(abs(X[k,3])>30){
  #     pdf[k,1]=0
  #   }else{
  #     alpha=0+(6*pi)/60*(X[k,3]-(-30))
  #     x=(X[k,3]+35)*cos(alpha)
  #     y=(X[k,3]+35)*sin(alpha)
  #     pdf[k,1]=dmvnorm(c(X[k,1:2]), mean = c(x,y), sigma = diag(2), log = T)
  #   }
  # }
  # logprior=pdf
  
  sigma=1

  # faster way of doing it by a bit
  n=dim(X)[1]; p=dim(X)[2]
  pdf=zeros(n,1)
  # pdf=matrix(apply(X,MARGIN = 1,FUN = function(x) ifelse(x[3]<0 || x[3]>12*pi, 0, 1/(2*pi*sigma)*exp(-1/(2*sigma)*((x[1]-x[3]/2*cos(x[3]/2))^2 + (x[2]-x[3]/2*sin(x[3]/2))^2)))))
  # logprior = log(60/(12*pi)*pdf)
  pdf=matrix(apply(X,MARGIN = 1,FUN = function(x) ifelse(x[3]<0 || x[3]>60, log(0), dmvnorm(c(x[1], x[2]), mean = c((x[3]+5)*cos(x[3]*(6*pi/60)), (x[3]+5)*sin(x[3]*(6*pi/60))), sigma = diag(2), log = T))))
  logprior=pdf
  # pdf=matrix(apply(X,MARGIN = 1,FUN = function(x) ifelse(x[3]<(-30) || x[3]>30, log(0), dmvnorm(c(x[1], x[2]), mean = c((x[3]+35)*cos((x[3]-(-30))*(6*pi/60)), (x[3]+35)*sin((x[3]-(-30))*(6*pi/60))), sigma = diag(2), log = T))))
  # logprior = pdf
  
  
  output=list(logprior,loglike)
  return(output)
}

ISEM_astro <- function(X,Proposal,Proposal0,df,Target_function,gama,j,data,w_thr){  # One iteration of EM with IS
  X=t_mix_sample(Proposal,X,df)
  
  # this part is deleting components that have too low of a threshold
  index_survival=which(Proposal$W>w_thr)
  if(length(index_survival)<Proposal$M){
    Proposal$M=length(index_survival)
    Proposal$Mu=Proposal$Mu[index_survival,]
    Proposal$Sigma=Proposal$Sigma[,,index_survival]
    Proposal$W=Proposal$W[,index_survival]
    Proposal$W=Proposal$W/sum(Proposal$W)
    
    Root_temp=X$Root
    for(i in 1:length(index_survival)){
      temp=which(Root_temp==index_survival[i])
      X$Root[temp]=i*ones(length(temp),1)
    }
  }
  
  X$Resp=matrix(0,X$N,Proposal$M) # responsibility of each component with regard to each particle
  for(i in 1:Proposal$M){
    X$Resp[,i]=exp(log_t_pdf(X$Values,Proposal$Mu[i,],drop(Proposal$Sigma[,,i]),df))
  }
  X$Proposal=X$Resp%*%t(Proposal$W)
  X$logProposal=log(X$Proposal)
  
  if(j>1){
    rp_temp=matrix(0,X$N,Proposal0$M)
    for(i in 1:Proposal0$M){
      rp_temp[,i]=exp(log_t_pdf(X$Values,Proposal0$Mu[i,],drop(Proposal0$Sigma[,,i]),df))
    }
    X$logProposal0=log(rp_temp%*%t(Proposal0$W))
  } else{
    X$logProposal0=X$logProposal
  }
  
  output=do.call(Target_function, args=list(X$Values,data))
  X$logPrior=output[[1]]; X$logLike=output[[2]]
  X$logTarget=X$logPrior+X$logLike
  X$logAnnealTarget=X$logTarget*gama[j]+X$logProposal0*(1-gama[j])
  
  X$logWeight=X$logAnnealTarget-X$logProposal
  logWeight_scaled=X$logWeight-max(X$logWeight)+10
  weight_temp=exp(logWeight_scaled)
  X$NormalizedWeight=weight_temp/sum(weight_temp)
  
  Proposal=t_mix_update_v2(X,Proposal,df)
  Proposal$W=Proposal$W/sum(Proposal$W)
  X=t_mix_sample(Proposal,X,df)
  
  # this part is deleting components that have too low of a threshold
  index_survival=which(Proposal$W>w_thr)
  if(length(index_survival)<Proposal$M){
    Proposal$M=length(index_survival)
    Proposal$Mu=Proposal$Mu[index_survival,]
    Proposal$Sigma=Proposal$Sigma[,,index_survival]
    Proposal$W=Proposal$W[,index_survival]
    Proposal$W=Proposal$W/sum(Proposal$W)
    
    Root_temp=X$Root
    for(i in 1:length(index_survival)){
      temp=which(Root_temp==index_survival[i])
      X$Root[temp]=i*ones(length(temp),1)
    }
  }
  
  X$Resp=matrix(0,X$N,Proposal$M)
  for(i in 1:Proposal$M){
    X$Resp[,i]=exp(log_t_pdf(X$Values,Proposal$Mu[i,],drop(Proposal$Sigma[,,i]),df))
  }
  
  X$Proposal=X$Resp%*%t(Proposal$W)
  X$logProposal=log(X$Proposal)
  
  rp_temp=matrix(0,X$N,Proposal0$M)
  for(i in 1:Proposal0$M){
    rp_temp[,i]=exp(log_t_pdf(X$Values,Proposal0$Mu[i,],drop(Proposal0$Sigma[,,i]),df))
  }
  X$logProposal0=log(rp_temp%*%t(Proposal0$W))
  
  output=do.call(Target_function, args=list(X$Values,data))
  X$logPrior=output[[1]]; X$logLike=output[[2]]
  X$logTarget=X$logPrior+X$logLike
  X$logAnnealTarget=X$logTarget*gama[j]+X$logProposal0*(1-gama[j])
  
  X$logWeight=X$logAnnealTarget-X$logProposal
  logWeight_scaled=X$logWeight-max(X$logWeight)+10
  weight_temp=exp(logWeight_scaled)
  X$NormalizedWeight=weight_temp/sum(weight_temp)
  
  output=list(X,Proposal)
  return(output)
}

log_t_pdf <- function(X,Mu,Sigma,df){
  # Ref: Robust mixture modelling using the t distribution;ML estimation of the t distribution using EM and its extensions,ECM and ECME
  #  X and Mu are row vectors
  # The input to logdet must be a symmetric positive definite matrix.
  if(dim(Sigma)[1]==1){
    Sigma=matrix(Sigma)
  }
  n = dim(X)[1]; p = dim(X)[2]
  pdf=log(gamma((df+p)/2))-.5*logdet(Sigma)-(.5*p*log(pi*df)+log(gamma(df/2))+.5*(df+p)*log(1+sqdist(t(X),t(Mu),inv_posdef(Sigma))/df));
  
  if(any(is.nan(pdf))){
    print("There was an NaN inside log_t_pdf")
    browser()
  }
  
  pdf=Re(pdf);
  return(pdf)
}

Merge <- function(Proposal,X,Cor_thr,df,d){
  #--------Components Merging-------------#
  Mixture_temp=Proposal
  z=repmat(Mixture_temp$W,X$N,1)*X$Resp/repmat(X$Proposal,1,Proposal$M)

  z_mean=as.numeric(t(X$NormalizedWeight)%*%z)
  cor=matrix(0,Mixture_temp$M,Mixture_temp$M)
  for (r in 1:Mixture_temp$M){
    for (c in (r):Mixture_temp$M){
      if (r == c){
        cor[r,c] = 0
      } else {
        cor[r,c]=sum((z[,r]-z_mean[r])*(z[,c]-z_mean[c])*X$NormalizedWeight)/(sqrt(sum((z[,r]-z_mean[r])^2*X$NormalizedWeight))*sqrt(sum((z[,c]-z_mean[c])^2*X$NormalizedWeight)))
      }
    }
  }

  m_ind=matrix(1,1,Mixture_temp$M)
  if (length(m_ind) > 1){
    for (r in 1:Mixture_temp$M){
      if (m_ind[r] == 1){
        for (c in r:Mixture_temp$M){
          if (m_ind[c] == 1){
            if (cor[r,c]>Cor_thr){
              m_ind[c]=0
              ind_pts_c=which(X$Root==c)
              Mixture_temp$Mu[r,]=(Mixture_temp$Mu[r,]*as.numeric(Mixture_temp$W[1,r])+Mixture_temp$Mu[c,]*as.numeric(Mixture_temp$W[1,c]))/as.numeric(Mixture_temp$W[1,r]+Mixture_temp$W[1,c])
              Mixture_temp$Sigma[,,r]=(drop(Mixture_temp$Sigma[,,r])*as.numeric(Mixture_temp$W[1,r])+drop(Mixture_temp$Sigma[,,c])*as.numeric(Mixture_temp$W[1,c]))/as.numeric(Mixture_temp$W[1,r]+Mixture_temp$W[1,c])

              if(!isposdef(drop(Mixture_temp$Sigma[,,r]))){
                if(d==7){
                  Mixture_temp$Sigma[,,r]=1e6*diag(c(1e-6, 1e-6, 1e-10, 1e-6, 1e-6, 1e-6, 1e-6))
                } else if(d==12){
                  Mixture_temp$Sigma[,,r]=1e6*diag(c(1e-6, 1e-6, 1e-10, 1e-6, 1e-6, 1e-6, 1e-6, 1e-10, 1e-6, 1e-6, 1e-6, 1e-6))
                } else if(d==17){
                  Mixture_temp$Sigma[,,r]=1e6*diag(c(1e-6, 1e-6, 1e-10, 1e-6, 1e-6, 1e-6, 1e-6, 1e-10, 1e-6, 1e-6, 1e-6, 1e-6, 1e-10, 1e-6, 1e-6, 1e-6, 1e-6))
                }
              }

              Mixture_temp$W[1,r]=sum(Mixture_temp$W[,c(r,c)])
              X$Resp[,r]=exp(log_t_pdf(X$Values,Mixture_temp$Mu[r,],drop(Mixture_temp$Sigma[,,r]),df))
              X$Root[ind_pts_c,1]=r*matrix(1,length(ind_pts_c),1)
            }
          }
        }
      }
    }
  }

  Sign=F
  index_remain=which(m_ind==1)
  M_temp=length(index_remain)
  if(M_temp<Proposal$M){
    Proposal$M=length(index_remain)
    Proposal$Mu=Mixture_temp$Mu[index_remain,]
    Proposal$Sigma=Mixture_temp$Sigma[,,index_remain]
    Proposal$W=Mixture_temp$W[,index_remain]
    X$Resp=X$Resp[,index_remain]
    X$Proposal=X$Resp%*%t(Proposal$W)
    X$logProposal=log(X$Proposal)
    X$logWeight=X$logAnnealTarget-X$logProposal
    Root_temp=X$Root
    for(i in 1:M_temp){
      ind_temp=which(Root_temp==index_remain[i])
      X$Root[ind_temp,1]=i*matrix(1,length(ind_temp),1)
    }
    logWeight_scaled=X$logWeight-max(X$logWeight)+10 # scaling
    weight_temp=exp(logWeight_scaled)
    X$NormalizedWeight=weight_temp/sum(weight_temp)
    Sign=T
  }
  output=list(Proposal,X,Sign)
  return(output)
}

t_mix_sample <- function(Proposal,X,df){
  M=Proposal$M
  W_m=Proposal$W
  Mu=Proposal$Mu
  # Sigma=drop(Proposal$Sigma)
  Sigma=Proposal$Sigma
  d=dim(Mu)[2]
  pts=matrix(0,X$N,d)
  root=matrix(0,X$N,1)
  for(i in 1:X$N){
    output=t_mixture_sampling(M,W_m,Mu,Sigma,df)
    pts[i,]=output[[1]]; root[i,1]=output[[2]]
  }
  X$Values=pts
  X$Root=root
  return(X)
}

t_mix_update_v2 <- function(X, Proposal, df){
  M=Proposal$M
  W_m=Proposal$W
  Mu=Proposal$Mu
  Sigma=Proposal$Sigma
  q=X$NormalizedWeight
  rp=X$Resp
  proposal=X$Proposal
  XV=X$Values
  N = dim(XV)[1]; d = dim(XV)[2]

  W_m_minus=matrix(0,1,M)
  Mu_minus=matrix(0,M,d)
  Sigma_minus=array(0, dim=c(d,d,M))

  t_mix_update_v2_varlist = c("W_m", "rp", "proposal", "q", "df", "d", "sqdist",
                              "XV", "Mu", "inv_posdef", "Sigma", "N", "cholproj",
                              "inv_triu", "logdet", "solve_triu", "sqdist",
                              "drop", "repmat")

  multiResultClass <- function(rou=NULL,W_m_minus=NULL,u=NULL,Mu_minus=NULL,Sigma_minus=NULL){
    # https://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop

    result = list(
      rou=rou,
      W_m_minus=W_m_minus,
      u=u,
      Mu_minus=Mu_minus,
      Sigma_minus=Sigma_minus
    )

    class(result) <- append(class(result), "multiResultClass")
    return(result)
  }

  cl <- makeCluster(4)
  registerDoParallel(cl)
  clusterExport(cl, varlist = t_mix_update_v2_varlist, envir = environment())

  oper <- foreach(ii = 1:M) %dopar% {
    # result <- multiResultClass()
    result = list()
    result$rou = as.numeric(W_m[,ii])*rp[,ii,drop=F]/proposal
    result$W_m_minus = as.numeric(t(result$rou)%*%q)
    result$u = (df+d)/(df + sqdist(t(XV), t(Mu[ii,,drop=F]), inv_posdef(drop(Sigma[,,ii]))))
    weight_temp = (result$rou*result$u*q)/as.numeric(t(result$rou*result$u)%*%q)
    (result$Mu_minus = t(weight_temp)%*%XV)
    matrix_temp=XV-repmat(result$Mu_minus,N,1)
    result$Sigma_minus=t(matrix_temp)%*%(matrix_temp*repmat(weight_temp,1,d))
    result$Sigma_minus = 0.5*(result$Sigma_minus + t(result$Sigma_minus))
    return(result)
  }
  stopCluster(cl)

  W_m_minus=matrix(oper[[1]]$W_m_minus)
  rou=oper[[1]]$rou
  u=oper[[1]]$u
  Mu_minus=oper[[1]]$Mu_minus
  Sigma_minus=array(oper[[1]]$Sigma_minus, dim = c(dim(Sigma_minus)[1], dim(Sigma_minus)[2], 1))

  if(M>1){
    for (i in 2:M){
      W_m_minus=cbind(W_m_minus,oper[[i]]$W_m_minus)
      rou=cbind(rou,oper[[i]]$rou)
      u=cbind(u,oper[[i]]$u)
      Mu_minus=rbind(Mu_minus,oper[[i]]$Mu_minus)
      Sigma_minus=abind(Sigma_minus,oper[[i]]$Sigma_minus,along=3)
    }
  }


  MixtureUpdated=list()
  MixtureUpdated$M=M
  MixtureUpdated$W=W_m_minus
  MixtureUpdated$Mu=Mu_minus
  MixtureUpdated$Sigma=Sigma_minus

  return(MixtureUpdated)
}

t_mix_update_v2_no_parallel <- function(X, Proposal, df){
  M=Proposal$M
  W_m=Proposal$W
  Mu=Proposal$Mu
  Sigma=Proposal$Sigma
  q=X$NormalizedWeight
  rp=X$Resp
  proposal=X$Proposal
  XV=X$Values

  N = dim(XV)[1]; d = dim(XV)[2]

  W_m_minus=matrix(0,1,M)
  Mu_minus=matrix(0,M,d)
  Sigma_minus=array(0, dim=c(d,d,M))
  rou=matrix(0,N,M)
  u=matrix(0,N,M)

  for(ii in 1:M){
    rou[,ii] = as.numeric(W_m[,ii])*rp[,ii]/proposal
    W_m_minus[,ii] = sum(t(rou[,ii])%*%q)
    u[,ii] = (df+d)/(df + sqdist(t(XV), t(Mu[ii,]), inv_posdef(drop(Sigma[,,ii]))))
    weight_temp = (rou[,ii]*u[,ii]*q)/(as.numeric(t(rou[,ii]*u[,ii])%*%q))
    Mu_minus[ii,] = t(weight_temp)%*%XV
    matrix_temp=XV-repmat(Mu_minus[ii,],N,1)
    Sigma_minus[,,ii]=t(matrix_temp)%*%(matrix_temp*repmat(weight_temp,1,d))
    Sigma_minus[,,ii] = 0.5*(drop(Sigma_minus[,,ii]) + t(drop(Sigma_minus[,,ii])))
  }

  MixtureUpdated=list()
  MixtureUpdated$M=M
  MixtureUpdated$W=W_m_minus
  MixtureUpdated$Mu=Mu_minus
  MixtureUpdated$Sigma=Sigma_minus

  return(MixtureUpdated)
}

t_mixture_pdf <- function(X, M, W, Mu, Sigma, df){
  #Ref: Robust mixture modelling using the t distribution; ML estimation of the t distribution using EM and its extensions,ECM and ECME
  n = dim(X)[1]

  cl <- makeCluster(4)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("t_pdf", "X", "Mu", "Sigma", "df", "drop"), envir = environment())

  Cpt_pdf = foreach(i = 1:M, .combine = "cbind") %dopar% {
    t_pdf(X, Mu[i,,drop=F], drop(Sigma[,,i]), df)
  }

  stopCluster(cl)

  return(as.numeric(Cpt_pdf%*%matrix(W, ncol = 1)))
}

t_mixture_pdf_no_parallel <- function(X, M, W, Mu, Sigma, df){
  #Ref: Robust mixture modelling using the t distribution; ML estimation of the t distribution using EM and its extensions,ECM and ECME
  n = dim(X)[1]

  Cpt_pdf = matrix(0, nrow = n, ncol = M)
  for (i in 1:M){
    Cpt_pdf[,i] = t_pdf(X, Mu[i,], drop(Sigma[,,i]), df)
  }
  return(as.numeric(Cpt_pdf%*%matrix(W, ncol = 1)))
}

t_mixture_sampling <- function(M,W_m,Mu,Sigma,df){
  `[`<-base::`[`
  r=W_m[1]
  d=dim(Mu)[2]
  rand_num=runif(1)
  for(ii in 1:M){
    if(rand_num <= r){
      tao=rgamma(1, shape=df/2, scale=1/(df/2))
      Sigma_temp=Sigma[,,ii]/tao
      if(is.vector(Sigma_temp)){
        Sigma_temp=matrix(Sigma_temp) # if 1D, make it a matrix
      }
      if(isposdef(Sigma_temp)){
        X=t(randnorm(1,t(Mu[ii,,drop=F]),NULL,Sigma_temp))
      } else{
        eig_values=eig(Sigma_temp)
        delta=runif(1,mean(eig_values),max(eig_values))
        Sigma_temp=Sigma_temp+1e-10*delta*diag(d)
        X=t(randnorm(1,t(Mu[ii,,drop=F]),NULL,Sigma_temp))
      }
      cmp_ind=ii
      break
    } else{
      r=r+W_m[ii+1]
    }
  }
  output=list(X,cmp_ind)
  return(output)
}

t_pdf <- function(X, Mu, Sigma, df){
  #Ref: Robust mixture modelling using the t distribution;ML estimation of the t distribution using EM and its extensions,ECM and ECME
  # X and Mu are row vectors
  dimX = dim(X)
  n = dimX[1]
  p = dimX[2]
  pdf = matrix(0, nrow = n, ncol = 1)
  invSig = solve(Sigma)
  for (k in 1:n){
    pdf[k,1] = gamma((df+p)/2)*(det(Sigma))^(-.5)/((pi*df)^(p/2)*gamma(df/2)*(1+((X[k,]-Mu)%*%invSig%*%t(X[k,]-Mu))/df)^((df+p)/2))
  }
  return(pdf)
}

#--------Target functions-----------#

Target_pdf_out_prod_density <- function(X){
  # Target density being a outer product of 7 univariate densities
  # Ref: Posterior-guided importance sampling for calculating marginal likelihoods
  dimX = dim(X)
  n = dimX[1]
  p = dimX[2]
  pdf = matrix(0, nrow = n, ncol = 1)
  for (k in 1:n){
    pdf[k,1]=1;
    pdf[k,1]=pdf[k,1]*(3/5*dgamma(X[k,1]+10,2,scale = 3)+2/5*dgamma(10-X[k,1],2,scale = 5));
    pdf[k,1]=pdf[k,1]*(3/4*dsn(X[k,2],3,1,5)+1/4*dsn(X[k,2],-3,3,-6))
    pdf[k,1]=pdf[k,1]*dt(X[k,3],4)
    pdf[k,1]=pdf[k,1]*(1/2*dbeta(X[k,4]+3,3,3)+1/2*dnorm(X[k,4],0,sqrt(1)))
    pdf[k,1]=pdf[k,1]*(1/2*dexp(X[k,5],1)+1/2*dexp(-X[k,5],1))
    pdf[k,1]=pdf[k,1]*dsn(X[k,6],0,8,-3)
    pdf[k,1]=pdf[k,1]*(1/8*dnorm(X[k,7],-10,sqrt(.1))+1/4*dnorm(X[k,7],0,sqrt(.15))+5/8*dnorm(X[k,7],7,sqrt(.2)))
  }
  return(pdf)
}

my_gaussian_mixture <- function(X, data){
  return(my_gaussian_mixture_fast(X,data)) # delete later
  
  loglike=0*data
  # Target density being a mixture of 2D gaussian distributions scaled to go
  # down to zero as you get further away from the origin
  
  n=dim(X)[1]; p=dim(X)[2]
  gridlength=4
  testx=matrix(linspace(-5.12,5.12,gridlength), ncol=1)
  testy=matrix(linspace(-5.12,5.12,gridlength), ncol=1)
  XY=meshgrid(testx,testy); testX=XY$X; testY=XY$Y
  Mewlikethecat=cbind(matrix(c(testX),ncol=1), matrix(c(testY),ncol=1))
  mySiggy=1e0/1*diag(p)
  pdf=matrix(0,n,1)
  
  #--------Cutoff version--------------#
  for(k in 1:n){
    if(norm(X[k,],type = "I")>5.12){ # this should be type = "I" but R is wrong
      pdf[k,1]=0
    } else{
      for (i in 1:dim(Mewlikethecat)[1]){
        pdf[k,1]=pdf[k,1]+dmvnorm(X[k,],Mewlikethecat[i,],mySiggy)
      }
    }
    pdf[k,1]=log(pdf[k,1])
  }
  
  # #-------Not cutoff version-----------#
  # for(k in 1:n){
  #   for (i in 1:dim(Mewlikethecat)[1]){
  #     pdf[k,1]=pdf[k,1]+dmvnorm(X[k,],Mewlikethecat[i,],mySiggy)
  #   }
  #   pdf[k,1]=log(pdf[k,1])
  # }
  
  logprior=pdf
  output=list(logprior,loglike)
  return(output)
}

my_gaussian_mixture_1d <- function(X, data){
  loglike=0*data
  X=matrix(X)
  # Target density being a mixture of 2D gaussian distributions scaled to go
  # down to zero as you get further away from the origin
  
  n=dim(X)[1]; p=dim(X)[2]
  gridlength=3
  testx=matrix(linspace(-5.12,5.12,gridlength), ncol=1)
  Mewlikethecat=testx
  mySiggy=1e0/1*diag(p)
  pdf=matrix(0,n,1)
  
  # #--------Cutoff version--------------#
  # for(k in 1:n){
  #   if(norm(X[k,],type = "I")>5.12){
  #     pdf[k,1]=0
  #   } else{
  #     for (i in 1:dim(Mewlikethecat)[1]){
  #       pdf[k,1]=pdf[k,1]+dmvnorm(X[k,],Mewlikethecat[i,],mySiggy)
  #     }
  #   }
  #   pdf[k,1]=log(pdf[k,1])
  # }
  
  #-------Not cutoff version-----------#
  for(k in 1:n){
    for (i in 1:dim(Mewlikethecat)[1]){
      pdf[k,1]=pdf[k,1]+dmvnorm(X[k,],Mewlikethecat[i,],mySiggy)
    }
    pdf[k,1]=log(pdf[k,1])
  }
  logprior=pdf
  
  output=list(logprior,loglike)
  return(output)
}

my_gaussian_mixture_fast <- function(X, data){
  loglike=0*data
  # Target density being a mixture of 2D gaussian distributions scaled to go
  # down to zero as you get further away from the origin
  
  n=dim(X)[1]; p=dim(X)[2]
  gridlength=4
  testx=matrix(linspace(-5.12,5.12,gridlength), ncol=1)
  testy=matrix(linspace(-5.12,5.12,gridlength), ncol=1)
  XY=meshgrid(testx,testy); testX=XY$X; testY=XY$Y
  Mewlikethecat=cbind(matrix(c(testX),ncol=1), matrix(c(testY),ncol=1))
  mySiggy=1e0/1*diag(p)
  Siggylist=replicate(dim(Mewlikethecat)[1],mySiggy,simplify = F)
  pdf=zeros(n,1)
  
  #--------Cutoff version--------------#
  for(k in 1:n){
    if(norm(X[k,],type = "1")>5.12){ # this should be type = "I" but R is wrong
      pdf[k,1]=0
    } else{
      pdf[k,1]=sum(mmvnpdfC(t(X[k,]),t(Mewlikethecat),Siggylist,Log=F))
    }
    pdf[k,1]=log(pdf[k,1])
  }
  logprior=pdf
  
  # #-------Not cutoff version-----------#
  # values=mmvnpdfC(t(X),t(Mewlikethecat),Siggylist,Log = F)
  # pdf=matrix(colSums(values), ncol=1)
  # logprior=log(pdf)
  
  output=list(logprior,loglike)
  return(output)
}

outerproduct <- function(X, data){
  `[`<-base::`[`
  loglike=0*data
  # Target density being an outer product of 1D densities
  
  dimX = dim(X)
  n = dimX[1]; p = dimX[2]
  pdf = matrix(0, nrow = n, ncol = 1)
  for (k in 1:n){
    pdf[k,1]=1;
    pdf[k,1]=pdf[k,1]*(3/5*dgamma(X[k,1]+10,2,scale = 3)+2/5*dgamma(10-X[k,1],2,scale = 5));
    pdf[k,1]=pdf[k,1]*(3/4*dsn(X[k,2],3,1,5)+1/4*dsn(X[k,2],-3,3,-6))
    pdf[k,1]=pdf[k,1]*dt(X[k,3],4)
    pdf[k,1]=pdf[k,1]*(1/2*dbeta(X[k,4]+3,3,3)+1/2*dnorm(X[k,4],0,sqrt(1)))
    pdf[k,1]=pdf[k,1]*(1/2*dexp(X[k,5],1)+1/2*dexp(-X[k,5],1))
    pdf[k,1]=pdf[k,1]*dsn(X[k,6],0,8,-3)
    pdf[k,1]=pdf[k,1]*(1/8*dnorm(X[k,7],-10,sqrt(.1))+1/4*dnorm(X[k,7],0,sqrt(.15))+5/8*dnorm(X[k,7],7,sqrt(.2)))
  }
  
  logprior=log(pdf)
  output=list(logprior,loglike)
  return(output)
}

Rastrigin <- function(X,data){
  loglike=0*data
  
  n=dim(X)[1]; p=dim(X)[2];
  pdf=10*p*ones(n,1)
  
  `[` <- function(...) base::`[`(...,drop=FALSE)
  
  for(k in 1:n){
    if(norm(X[k,],"O")>5.12){
      pdf[k]=0
    }else{
      for(i in 1:p){
        pdf[k]=pdf[k]+(X[k,i]^2-10*cos(2*pi*X[k,i]))
      }
    }
    pdf[k]=log(pdf[k])
  }
  logprior=pdf
  
  output=list(logprior,loglike)
  return(output)
}

LinearModelML_exact <- function(X,Y){
  S2=var(Y)*(n-1)
  Q=1/S2*t(Y)%*%(diag(n) - X%*%solve(t(X)%*%X)%*%t(X))%*%Y; Q=as.numeric(Q)
  
  output = gamma((n-1)/2)/(k*sqrt(n)*(pi*S2*Q)^((n-1)/2))*(k/(n+1))^((k-1)/2)*hypergeo(k/2,(n-1)/2,(k+2)/2,k*(1-1/Q)/(n+1))
  
  return(as.numeric(Re(output)))
}

LinearModelHald <- function(thetas, data){
  loglike=log(Linearmodel_likelihood(thetas,data))
  logprior=log(Linearmodel_prior(thetas))
  
  output=list(logprior,loglike)
  return(output)
}

LinearModel_likelihood <- function(thetas,data){
  X = matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 1, 11, 
               11, 7, 11, 3, 1, 2, 21, 1, 11, 10, 26, 29, 56, 31, 52, 55, 71, 
               31, 54, 47, 40, 66, 68, 6, 15, 8, 8, 6, 9, 17, 22, 18, 4, 23, 
               9, 8, 60, 52, 20, 47, 33, 22, 6, 44, 22, 26, 34, 12, 12),ncol=5)
  # X2 = X[,-1]
  BigSigma22 = matrix(c(0.0927104018856017, 0.0856862094001623, 0.0926373565871297, 
                        0.084454955318979, 0.0856862094001621, 0.0875602572156389, 0.0878666396937198, 
                        0.0855980995271992, 0.0926373565871298, 0.0878666396937201, 0.0952014097401136, 
                        0.0863919188425777, 0.0844549553189787, 0.0855980995271992, 0.0863919188425775, 
                        0.0840311911923554),
                      nrow = 4)
  n = 13; k = 5;
  
  Y = matrix(data,ncol=1)
  
  # X = as.matrix(unname(Hald[,1:4]))
  # n=dim(X)[1]
  # X=cbind(ones(n,1),X)
  # n=dim(X)[1]; k=dim(X)[2]
  
  # Beta=thetas[,1:5]
  # sigmasq=matrix(thetas[,6])
  # g=matrix(thetas[,7])
  
  Beta=thetas[,1:5]
  eta=matrix(thetas[,6])
  xi=matrix(thetas[,7])
  
  # transform the data
  sigmasq=exp(eta)
  g=exp(xi)-1+(n+1)/k
  
  return(1/(2*pi*sigmasq)^(n/2)*exp(-1/(2*sigmasq)*norm(Y-X%*%Beta, type = "2")^2))
}

LinearModel_prior <- function(thetas){
  X = matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 1, 11, 
               11, 7, 11, 3, 1, 2, 21, 1, 11, 10, 26, 29, 56, 31, 52, 55, 71, 
               31, 54, 47, 40, 66, 68, 6, 15, 8, 8, 6, 9, 17, 22, 18, 4, 23, 
               9, 8, 60, 52, 20, 47, 33, 22, 6, 44, 22, 26, 34, 12, 12),ncol=5)
  # X2 = X[,-1]
  BigSigma22 = matrix(c(0.0927104018856017, 0.0856862094001623, 0.0926373565871297, 
                        0.084454955318979, 0.0856862094001621, 0.0875602572156389, 0.0878666396937198, 
                        0.0855980995271992, 0.0926373565871298, 0.0878666396937201, 0.0952014097401136, 
                        0.0863919188425777, 0.0844549553189787, 0.0855980995271992, 0.0863919188425775, 
                        0.0840311911923554),
                      nrow = 4)
  n = 13; k = 5;
  # X = as.matrix(unname(Hald[,1:4]))
  # n=dim(X)[1]
  # X=cbind(ones(n,1),X)
  # n=dim(X)[1]; k=dim(X)[2]
  
  # Beta=thetas[,1:5]
  # sigmasq=matrix(thetas[,6])
  # g=matrix(thetas[,7])
  
  Beta=thetas[,1:5]
  eta=matrix(thetas[,6])
  xi=matrix(thetas[,7])
  
  # transform the data
  sigmasq=exp(eta)
  g=exp(xi)-1+(n+1)/k
  
  if (g > (1+n)/k - 1){
    indicatorg = 1
  } else {
    indicatorg = 0
  }
  
  return(1/sigmasq*dmvnorm(Beta2,sigma = g*sigmasq*BigSigma22)*(1+g)^(-3/2)*indicatorg*exp(xi))
}
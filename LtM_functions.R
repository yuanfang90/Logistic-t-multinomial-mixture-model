generate_data <- function(num_grp,n,true_mu,true_Sig,M,nu){
  Y <- list()
  W <- list()
  U<-list()
  
  z<-t(rmultinom(n=n,size=1,prob=pi_g))
  ng<-colSums(z)
  true_lab<-c(rep(1,ng[1]),rep(2,ng[2]))
  
  
  for (g in 1:num_grp) {
    # Y[[g]] <- rmvt(num_observation[g],mu = true_mu[[g]],S = true_Sig[[g]]/nu[[g]], df = nu[[g]])
    U[[g]]<-as.matrix(rgamma(ng[g],shape=nu[g]/2,rate=nu[g]/2))
    Y[[g]]<-matrix(nrow=ng[g],ncol=K)
    W[[g]]<-matrix(nrow=ng[g],ncol=K+1)
    for(i in 1:ng[g]){
      Y[[g]][i,] <- rmvnorm(1, true_mu[[g]], true_Sig[[g]]/U[[g]][i,])
      theta <- c(exp(Y[[g]][i,]),1)/sum(exp(Y[[g]][i,])+1)
      W[[g]][i,] <- rmultinom(1,sample(c(M/2):M,1),prob=theta)
    }
  }
  
  Y <- Reduce(rbind,Y)
  W <- Reduce(rbind,W)
  U<-Reduce(rbind,U)
  return(list(W=W,Y=Y,label=true_lab,U=U))
}


### update V and m with newton method
## first and second derivative of elbo wrt v for the kth dim
vk_fun <- function(x,iSig_kk,m_k,M,xi,K,u){
  return(1/x-x*iSig_kk*K*u-x*M*exp(m_k+x^2/2)/xi) 
}
vk_fun.prime <- function(x,iSig_kk,m_k,M,xi,K,u){
  return(-1/x^2-iSig_kk*K*u-(x^2+1)*M*exp(m_k+x^2/2)/xi)
}

## first and second derivative of elbo wrt m as a vector
mk_fun <- function(x,W_i,Sig,mu,vsquare,M,xi,iSig,u){
  return(W_i-iSig%*%(x-mu)*u-M*exp(x+vsquare/2)/xi)
}
mk_fun.jacobian <- function(x,Sig,vsquare,M,xi,iSig,u){
  return(-iSig*u-M*diag(exp(x+vsquare/2))/xi)
}

## newton's method with one update
newton_up_m <- function(xold,W_i,Sig,mu,vsquare,M,xi,K,iSig,u){
  xnew <- c(xold-solve(mk_fun.jacobian(xold,Sig,vsquare,M,xi,iSig,u),tol=.Machine$double.neg.eps)%*%mk_fun(xold,W_i,Sig,mu,vsquare,M,xi,iSig,u))
  xnew[(K+1)] <- 0
  return(xnew)
}

newton_up_v <- function(xold,iSig_kk,m_k,M,xi,K,u){
  xnew <- xold-vk_fun(xold,iSig_kk,m_k,M,xi,K,u)/vk_fun.prime(xold,iSig_kk,m_k,M,xi,K,u)
  return(xnew)
}

## update variational prameter
varpar_fun <- function(W,m,V,mu,Sig,K,z,iSig,u){
  n=nrow(W)
  xi_hat <- numeric(n)
  m_hat <- matrix(NA,nrow=n,ncol=K+1)
  V_hat <- vector("list",n)
  Vmat<-matrix(NA,nrow=n,ncol=K+1)
  # forp <- colSums(z,na.rm=TRUE)/sum(z,na.rm=TRUE)
  for(i in 1:n){
    # g=ifelse(any(is.na(z[i,])),forp,which.max(z[i,]))
    g=which.max(z[i,])
    ## update xi
    xi_hat[i] <- sum(exp(m[i,]+diag(V[[i]])/2))
    ## update m
    m_hat[i,] <- newton_up_m(xold=m[i,],W_i=W[i,],Sig=Sig[[g]],mu=mu[g,],vsquare=diag(V[[i]]),M=sum(W[i,]),xi=xi_hat[i],K=K,iSig=iSig[[g]],u=u[i,g])
    ## update V
    for_V <- rep(0,(K+1))
    for(k in 1:K){
      for_V[k] <- newton_up_v(xold=sqrt(diag(V[[i]])[k]),iSig_kk=iSig[[g]][k,k],m_k=m[i,k],M=sum(W[i,]),xi=xi_hat[i],K=K,u=u[i,g])
    }
    V_hat[[i]] <- diag((for_V)^2)
    Vmat[i,]<-diag(V_hat[[i]])
  }
  
  return(list(xi=xi_hat,m=m_hat,V=V_hat,Vmat=Vmat))
}

#### M-step
mst_fun <- function(m,V,G,K,z,u){
  mu <- matrix(NA,nrow=G,ncol=(K+1))
  Sig <- list()
  for(g in 1:G){
    temp <- try(cov.wt(m, wt=(z[,g]*u[,g]), center=TRUE, method="ML"), silent=TRUE)
    if(!class(temp) == "try-error"){
      mu[g,] <- temp$center
      forsig <- Reduce(`+`,Map(`*`,V,(z[,g]*u[,g]/sum(z[,g]*u[,g]))))
      Sig[[g]] <- temp$cov+(K+1)*forsig
    }
  }
  return(list(mu=mu,Sig=Sig))
}


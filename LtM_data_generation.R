rm(list=ls())

require(mvtnorm)
library(MASS)
library(mclust)

G <- 2
true_nu <- c(10, 200)
K <- 3
t_sigma_1 <- matrix(c(0.5,-0.22,0,-0.22,0.3,-0.15,0,-0.15,0.3),K,K)
t_sigma_2 <- matrix(c(0.7,-0.2,0.5,-0.2,1,0,0.5,0,0.7),K,K)
true_sig<-list(t_sigma_1,t_sigma_2)

t_mu_1 <- c(2,3,4)
t_mu_2 <- c(2,3,0)
true_mu <- list(t_mu_1, t_mu_2)

pi_g<-c(0.5,0.5)
n <- 1000


Data <- generate_data(num_grp=G,n=n,true_mu=true_mu,true_Sig=true_sig,M=5000,nu=true_nu)

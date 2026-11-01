
rm(list=ls())
library(snowfall)
library(rlecuyer)

parallel_func=function(arg){
library(copula)
library(nleqslv)
library(mvtnorm)
library(QRM)
  ### Kernel function Kbar for standard normal survival function
  Kbar <- function(z) {1 - pnorm(z)}
  
  ### Functions for estimation
  
  
 
# (ii) Estimate VaR using the given Kbar and h
# X: n by d; 
  estimate_VaR <- function(c_vec,X, alpha, h, wgt) {
n=dim(X)[1];d=dim(X)[2]
c_vec=matrix(c_vec,nrow=d,ncol=1)
c_vec[is.na(c_vec)]=0
wgt=matrix(wgt,nrow=n,ncol=1)
values=X%*%(1-c_vec) 
ini_value=quantile(values, probs=alpha)
#print(c_vec)
    target <- 1 - alpha
    f <- function(x) {t(wgt)%*%matrix((Kbar((x - values)/h)-target),nrow=n)}
    lower_bound <- min(values) - 10*h
    upper_bound <- max(values) + 10*h
    sol <- try(uniroot(f, interval=c(lower_bound, upper_bound)),silent=T)
    res=ini_value
    if(length(sol)>1){
    res=sol$root}
    return(res)
    
    }
  
# (iii*) Solve for (c, lambda_I) for the individual reinsurance (ES) with jacobian (fast)
  solve_c_lambda_I <- function(X, beta_vec, tilde_beta_vec, mu0, alpha, h, 
                               c_init, lambda_init, wgt) {
    n<- dim(X)[1];d=dim(X)[2]
    beta_vec=matrix(beta_vec,nrow=d,ncol=1)
    tilde_beta_vec=matrix(tilde_beta_vec,nrow=d,ncol=1)
wgt=matrix(wgt,nrow=n,ncol=1)
    mu_hat <- t(X)%*%wgt/sum(wgt)
   
    obj_fun <- function(par) {
      c_vec <- matrix(exp(par[1:d])/(1+exp(par[1:d])),nrow=d,ncol=1)
      lambda_I <- par[d+1]
      gamma_I <- estimate_VaR(c_vec,X, alpha, h, wgt)
     
      weights <- matrix(t(wgt)*matrix(Kbar((gamma_I - X%*%(1-c_vec))/h),nrow=1),nrow=n)
     
       # L_{I,es,1}: vector of length d
      L_I_es_1 <- -t(X)%*%weights/sum(wgt)/(1-alpha)+mu_hat+lambda_I*diag(as.vector(tilde_beta_vec))%*%mu_hat


      # L_{I,es,2}: scalar
      L_I_es_2 <- t(c_vec)%*%diag(as.vector(tilde_beta_vec))%*%mu_hat-mu0-t(beta_vec)%*%mu_hat

      return(c(L_I_es_1, L_I_es_2))
    }
    
    start_val <- c(log(c_init/(1-c_init)), lambda_init)
    
    c_vec_hat=c_init
    
    sol <- try(nleqslv(start_val, obj_fun),silent=T)
   
    if(length(sol)>1){
    c_vec_hat <- exp(sol$x[1:d])/(1+exp(sol$x[1:d]))
    }
    
return(c_vec_hat)    
#return(list(sol=sol, c_hat=c_hat, lambda_I_hat=lambda_I_hat, convergence=sol$convergence))
  }
  
estimate_theta=function(X,beta,tilde_beta,beta_vec, tilde_beta_vec, mu0, alpha, wgt){
n=dim(X)[1];d=dim(X)[2]
beta_vec=matrix(beta_vec,nrow=d,ncol=1)
tilde_beta_vec=matrix(tilde_beta_vec,nrow=d,ncol=1)
wgt=matrix(wgt,nrow=n,ncol=1)
m=as.integer(n*0.5);
#s_L=(1:m)*2-1;s_U=(1:m)*2;
s_L=1:m;s_U=(1+m):n
X_L=X[s_L,];X_U=X[s_U,]
wgt_L=wgt[s_L,]
wgt_U=wgt[s_U,]
mu_hat <- t(X)%*%wgt/sum(wgt)



##the same premium
c_hat=(sum(beta_vec*mu_hat)+mu0)/(tilde_beta*sum(mu_hat))


h_n<- sd(rowSums(X))*(1-c_hat)*n^(-0.33)
h_m=sd(rowSums(X))*(1-c_hat)*m^(-0.33)

gamma_hat_J=estimate_VaR(matrix(c_hat,nrow=d,ncol=1),X,alpha,h_n,wgt)
gamma_L_hat_J=estimate_VaR(matrix(c_hat,nrow=d,ncol=1),X_L,alpha,h_m,wgt_L)
X_J=matrix(apply((1-c_hat)*X,1,sum),ncol=1)
X_L_J=matrix(apply((1-c_hat)*X_L,1,sum),ncol=1)

part1=1/(1-alpha)*sum(wgt*X_J*Kbar((gamma_hat_J-X_J)/h_n))/sum(wgt)
part1_L=1/(1-alpha)*sum(wgt_L*X_L_J*Kbar((gamma_L_hat_J-X_L_J)/h_m))/sum(wgt_L)
part2=-(1-c_hat)*sum(mu_hat)+mu0


c_init_vec <- rep(c_hat, d)
lambda_init <- 1

c_vec_hat=solve_c_lambda_I(X, beta_vec, tilde_beta_vec, mu0, alpha, h_n, c_init_vec, lambda_init, wgt)


gamma_hat_I=estimate_VaR(c_vec_hat,X,alpha,h_n,wgt)
gamma_U_hat_I=estimate_VaR(c_vec_hat,X_U,alpha,h_m,wgt_U)



X_I=X%*%matrix(1-c_vec_hat,ncol=1)
X_U_I=X_U%*%matrix(1-c_vec_hat,ncol=1)



part3=1/(1-alpha)*sum(wgt*X_I*Kbar((gamma_hat_I-X_I)/h_n))/sum(wgt)
part3_U=1/(1-alpha)*sum(wgt_U*X_U_I*Kbar((gamma_U_hat_I-X_U_I)/h_m))/sum(wgt_U)
part4=-sum((1-c_vec_hat)*mu_hat)+mu0


rho_I=part3+part4
rho_U_I=part3_U+part4
rho_J=part1+part2
rho_L_J=part1_L+part2
theta_hat1=log(rho_J/rho_I)
theta_hat2=log(rho_L_J/rho_U_I)
theta_hat=1/2*(theta_hat1+theta_hat2)

res=(c(theta_hat1,theta_hat2,theta_hat,rho_J,rho_L_J,rho_I,rho_U_I,c_hat,c_vec_hat))

return(res)
} 


 
#####
set.seed(102025)

d=3;alpha=0.95;df.nu=3;



####
simu=function(n,N.sim,Boot,Case,Sigma,Scale,p.shape){
  if(Sigma==1){
    p.sigma <- array(c(1,0.25,0.25,0.25,1,0.25,0.25,0.25,1),dim=c(3,3))
  }
  
  if(Scale==1){p.scale <- c(1000, 1000, 1000);}
  if(Scale==2){p.scale <- c(800, 1000, 1200);}
  
  if(Case==1){
    beta_scalar <- 0.05*0; 
    beta_vec <- c(0.05, 0.05, 0.05);
    tilde_beta_scalar <- 0.1;
    tilde_beta_vec <- c(0.1, 0.1, 0.1);
  }
  if(Case==2){
    beta_scalar <- 0; 
    beta_vec <-c(0.05,0.05,0.05);
    tilde_beta_scalar <- 0.1; tilde_beta_vec <- c(0.1+0.01*p.scale[2]/p.scale[1], 0.1, 0.1+0.01*p.scale[2]/p.scale[3]);
    }
  if(Case==3){beta_scalar <- 0; 
  beta_vec <-c(0.05,0.05,0.05);
  tilde_beta_scalar <- 0.1; tilde_beta_vec <- c(0.1-0.01*p.scale[2]/p.scale[1], 0.1, 0.1-0.01*p.scale[2]/p.scale[3]);
  }
  
  if(Case==4){
    beta_scalar <- 0; 
    beta_vec <-c(0.05,0.05,0.05);
    tilde_beta_scalar <- 0.1; tilde_beta_vec <- c(0.1+0.01*p.scale[2]/p.scale[1], 0.1, 0.1-0.01*p.scale[2]/p.scale[3]);
  }
  
  
  mu_true=p.scale/(p.shape-1)
mu0=0.45*tilde_beta_scalar*sum(mu_true)-sum(beta_vec*mu_true)
  

  res=matrix(NA, nrow=N.sim,ncol=(8+d)*(1+Boot))
for(j in 1:N.sim){
        #U <- rmvnorm(n, mean=rep(0,d), sigma=Sigma)
        # Transform to uniform(0,1) using pnorm
        #U_unif <- pnorm(U)
  
  U=rcopula.t(n, df=df.nu,p.sigma)
  U_unif=U
X <- matrix(NA, nrow=n,ncol=d)
        for (s in 1:d) {X[,s] <- p.scale[s]*((1-U_unif[,s])^(-1/p.shape)-1)}
        wgt=rep(1,n) # not bootstrap yet!
res[j,1:(8+d)]=estimate_theta(X,beta_scalar,tilde_beta_scalar,beta_vec, tilde_beta_vec, mu0, alpha,wgt)
if(Boot>0){
for(b in 1:Boot){
W=rexp(n)
res[j,(1+(8+d)*b):((8+d)*(1+b))]=estimate_theta(X,beta_scalar,tilde_beta_scalar,beta_vec, tilde_beta_vec, mu0, alpha,W)
}
}
#print(c(j,res[j,1:(8+d)]))
}
return(res)
}
####
n=100000;N.sim=10000;Boot=0;
sim.true=matrix(NA,nrow=1,ncol=8+d);
if(arg==1){
Case=1;Sigma=1;Scale=1;p.shape=3.5;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'));

par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()

pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()




colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}

if(arg==2){
Case=2;Sigma=1;Scale=1;p.shape=3.5;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}

if(arg==3){
Case=3;Sigma=1;Scale=1;p.shape=3.5;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}
if(arg==4){
  Case=4;Sigma=1;Scale=1;p.shape=3.5;
  res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
  sim.true[1,]=apply(res,2,mean)
  pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
  par(mfrow=c(1,3))
  plot(density(res[,1]),xlab="",main="hat_theta1")
  plot(density(res[,2]),xlab="",main="hat_theta2")
  plot(density(res[,3]),xlab="",main="hat_theta3")
  dev.off()
  pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
  print(sim.true)
  write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}
if(arg==5){
Case=1;Sigma=1;Scale=2;p.shape=3.5;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}
if(arg==6){
Case=2;Sigma=1;Scale=2;p.shape=3.5;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}

if(arg==7){
Case=3;Sigma=1;Scale=2;p.shape=3.5;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}
if(arg==8){
  Case=4;Sigma=1;Scale=2;p.shape=3.5;
  res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
  sim.true[1,]=apply(res,2,mean)
  pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
  par(mfrow=c(1,3))
  plot(density(res[,1]),xlab="",main="hat_theta1")
  plot(density(res[,2]),xlab="",main="hat_theta2")
  plot(density(res[,3]),xlab="",main="hat_theta3")
  dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

  
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
  print(sim.true)
  write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}
if(arg==9){
Case=1;Sigma=1;Scale=1;p.shape=6;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}
if(arg==10){
Case=2;Sigma=1;Scale=1;p.shape=6;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}
if(arg==11){
Case=3;Sigma=1;Scale=1;p.shape=6;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}
if(arg==12){
  Case=4;Sigma=1;Scale=1;p.shape=6;
  res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
  sim.true[1,]=apply(res,2,mean)
  pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
  par(mfrow=c(1,3))
  plot(density(res[,1]),xlab="",main="hat_theta1")
  plot(density(res[,2]),xlab="",main="hat_theta2")
  plot(density(res[,3]),xlab="",main="hat_theta3")
  dev.off()
  
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
  print(sim.true)
  write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}

if(arg==13){
Case=1;Sigma=1;Scale=2;p.shape=6;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}

if(arg==14){
Case=2;Sigma=1;Scale=2;p.shape=6;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()

pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}
if(arg==15){
Case=3;Sigma=1;Scale=2;p.shape=6;
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()


colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}

if(arg==16){
  Case=4;Sigma=1;Scale=2;p.shape=6;
  res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
  sim.true[1,]=apply(res,2,mean)
  pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf'))
  par(mfrow=c(1,3))
  plot(density(res[,1]),xlab="",main="hat_theta1")
  plot(density(res[,2]),xlab="",main="hat_theta2")
  plot(density(res[,3]),xlab="",main="hat_theta3")
  dev.off()
  
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ1.pdf'));
qqnorm(res[,1]);qqline(res[,1])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ2.pdf'));
qqnorm(res[,2]);qqline(res[,2])
dev.off()
pdf(file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'_QQ3.pdf'));
qqnorm(res[,3]);qqline(res[,3])
dev.off()


colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
  print(sim.true)
  write.csv(sim.true,file=paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv'))
}
  
 }
#job_id=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#sfInit(parallel = TRUE, cpus =1, type = 'SOCK')
#sfClusterSetupRNGstream(seed = rep(job_id,4))
#res <- sfLapply(job_id, parallel_func)
#sfStop()
for(job_id in 1:16){
  set.seed(rep(job_id,4));parallel_func(job_id)
}


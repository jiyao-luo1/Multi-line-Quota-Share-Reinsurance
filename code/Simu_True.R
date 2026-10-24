# Portable public copy; estimation and historical random generation are unchanged.
# Run from the repository root with Rscript code/Simu_True.R.
script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file with Rscript.")
repository_root <- dirname(dirname(normalizePath(sub("^--file=", "", script_argument))))
env_integer <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || !is.finite(value) || value != floor(value) || value < minimum)
    stop(name, " must be an integer >= ", minimum, ".")
  as.integer(value)
}
env_number <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || !is.finite(value)) stop(name, " must be finite and numeric.")
  value
}
output_directory <- function(default) {
  value <- Sys.getenv("REINSURANCE_OUT", default)
  if (!grepl("^(/|[A-Za-z]:[/\\])", value)) value <- file.path(repository_root, value)
  dir.create(value, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(value)) stop("Cannot create output directory: ", value)
  value
}
for (package in c("QRM", "nleqslv")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing R package: ", package)
}
library(QRM)
library(nleqslv)

OUT <- output_directory("outputs/true")



parallel_func=function(arg){
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
    
    #lambda_I_hat <- sol$x[d+1]
    #print(is.na(c_vec_hat))
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


#c_hat=(beta*sum(mu_hat)+mu0)/(tilde_beta*sum(mu_hat))
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
  if(Sigma==2){p.sigma <- array(c(1,0,-0.25,0,1,-0.25,-0.25,-0.25,1),dim=c(3,3))
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
    beta_scalar <- 0.05*0; 
    beta_vec <-c(0.05,0.05,0.05)+0*c(0.05+0.02*p.scale[2]/p.scale[1], 0.05-0.01, 0.05+0.02*p.scale[2]/p.scale[3]);
    tilde_beta_scalar <- 0.1; tilde_beta_vec <- c(0.1+0.01*p.scale[2]/p.scale[1], 0.1, 0.1+0.01*p.scale[2]/p.scale[3]);
    }
  if(Case==3){beta_scalar <- 0.05*0; 
  beta_vec <-c(0.05,0.05,0.05)+0*c(0.05-0.02*p.scale[2]/p.scale[1], 0.05+0.01, 0.05-0.02*p.scale[2]/p.scale[3]);
  tilde_beta_scalar <- 0.1; tilde_beta_vec <- c(0.1-0.01*p.scale[2]/p.scale[1], 0.1, 0.1-0.01*p.scale[2]/p.scale[3]);
  }
  
  if(Case==4){
    beta_scalar <- 0.05*0; 
    beta_vec <-c(0.05,0.05,0.05)+0*c(0.05+0.05*p.scale[2]/p.scale[1], 0.05-0.01, 0.05+0.05*p.scale[2]/p.scale[3]);
    tilde_beta_scalar <- 0.1; tilde_beta_vec <- c(0.1+0.01*p.scale[2]/p.scale[1], 0.1, 0.1-0.01*p.scale[2]/p.scale[3]);
  }
  
  
  
  mu0=(0.2*tilde_beta_scalar-beta_scalar)/(p.shape-1)*sum(p.scale)

  res=matrix(NA, nrow=N.sim,ncol=(8+d)*(1+Boot))
for(j in 1:N.sim){
        #U <- rmvnorm(n, mean=rep(0,d), sigma=Sigma)
        # Transform to uniform(0,1) using pnorm
        #U_unif <- pnorm(U)
  
  U=rcopula.t(n, df=df.nu,p.sigma)
  U_unif=pt(U,df=df.nu)
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
n <- env_integer("REINSURANCE_TRUE_N", 100000L, 4L)
N.sim <- env_integer("REINSURANCE_TRUE_REPS", 10000L, 2L)
Boot <- 0L
if (n %% 2L != 0L) stop("REINSURANCE_TRUE_N must be even for equal halves.")
sim.true=matrix(NA,nrow=1,ncol=8+d);
if(arg==1){
Case <- env_integer("REINSURANCE_CASE", 1L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 1L)
p.shape <- env_number("REINSURANCE_SHAPE", 3.5)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')));

par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()

colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}

if(arg==2){
Case <- env_integer("REINSURANCE_CASE", 2L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 1L)
p.shape <- env_number("REINSURANCE_SHAPE", 3.5)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}

if(arg==3){
Case <- env_integer("REINSURANCE_CASE", 3L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 1L)
p.shape <- env_number("REINSURANCE_SHAPE", 3.5)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}
if(arg==4){
  Case <- env_integer("REINSURANCE_CASE", 4L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 1L)
p.shape <- env_number("REINSURANCE_SHAPE", 3.5)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
  res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
  sim.true[1,]=apply(res,2,mean)
  pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
  par(mfrow=c(1,3))
  plot(density(res[,1]),xlab="",main="hat_theta1")
  plot(density(res[,2]),xlab="",main="hat_theta2")
  plot(density(res[,3]),xlab="",main="hat_theta3")
  dev.off()
  colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
  print(sim.true)
  write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}
if(arg==5){
Case <- env_integer("REINSURANCE_CASE", 1L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 2L)
p.shape <- env_number("REINSURANCE_SHAPE", 3.5)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}
if(arg==6){
Case <- env_integer("REINSURANCE_CASE", 2L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 2L)
p.shape <- env_number("REINSURANCE_SHAPE", 3.5)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}

if(arg==7){
Case <- env_integer("REINSURANCE_CASE", 3L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 2L)
p.shape <- env_number("REINSURANCE_SHAPE", 3.5)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}
if(arg==8){
  Case <- env_integer("REINSURANCE_CASE", 4L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 2L)
p.shape <- env_number("REINSURANCE_SHAPE", 3.5)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
  res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
  sim.true[1,]=apply(res,2,mean)
  pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
  par(mfrow=c(1,3))
  plot(density(res[,1]),xlab="",main="hat_theta1")
  plot(density(res[,2]),xlab="",main="hat_theta2")
  plot(density(res[,3]),xlab="",main="hat_theta3")
  dev.off()
  colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
  print(sim.true)
  write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}
if(arg==9){
Case <- env_integer("REINSURANCE_CASE", 1L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 1L)
p.shape <- env_number("REINSURANCE_SHAPE", 6)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}
if(arg==10){
Case <- env_integer("REINSURANCE_CASE", 2L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 1L)
p.shape <- env_number("REINSURANCE_SHAPE", 6)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}
if(arg==11){
Case <- env_integer("REINSURANCE_CASE", 3L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 1L)
p.shape <- env_number("REINSURANCE_SHAPE", 6)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}
if(arg==12){
  Case <- env_integer("REINSURANCE_CASE", 4L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 1L)
p.shape <- env_number("REINSURANCE_SHAPE", 6)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
  res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
  sim.true[1,]=apply(res,2,mean)
  pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
  par(mfrow=c(1,3))
  plot(density(res[,1]),xlab="",main="hat_theta1")
  plot(density(res[,2]),xlab="",main="hat_theta2")
  plot(density(res[,3]),xlab="",main="hat_theta3")
  dev.off()
  colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
  print(sim.true)
  write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}

if(arg==13){
Case <- env_integer("REINSURANCE_CASE", 1L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 2L)
p.shape <- env_number("REINSURANCE_SHAPE", 6)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}

if(arg==14){
Case <- env_integer("REINSURANCE_CASE", 2L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 2L)
p.shape <- env_number("REINSURANCE_SHAPE", 6)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}
if(arg==15){
Case <- env_integer("REINSURANCE_CASE", 3L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 2L)
p.shape <- env_number("REINSURANCE_SHAPE", 6)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
sim.true[1,]=apply(res,2,mean)
pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
par(mfrow=c(1,3))
plot(density(res[,1]),xlab="",main="hat_theta1")
plot(density(res[,2]),xlab="",main="hat_theta2")
plot(density(res[,3]),xlab="",main="hat_theta3")
dev.off()
colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
print(sim.true)
write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}

if(arg==16){
  Case <- env_integer("REINSURANCE_CASE", 4L)
Sigma <- env_integer("REINSURANCE_SIGMA", 1L)
Scale <- env_integer("REINSURANCE_SCALE", 2L)
p.shape <- env_number("REINSURANCE_SHAPE", 6)
if (!Case %in% 1:4 || !Sigma %in% 1:2 || !Scale %in% 1:2 || p.shape <= 1)
  stop("Require CASE in 1:4, SIGMA and SCALE in 1:2, and SHAPE > 1.")
  res=simu(n,N.sim,Boot,Case,Sigma,Scale,p.shape)
  sim.true[1,]=apply(res,2,mean)
  pdf(file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.pdf')))
  par(mfrow=c(1,3))
  plot(density(res[,1]),xlab="",main="hat_theta1")
  plot(density(res[,2]),xlab="",main="hat_theta2")
  plot(density(res[,3]),xlab="",main="hat_theta3")
  dev.off()
  colnames(sim.true)=c("theta1","theta2","theta","rho_J","rho_L_J","rho_I","rho_U_I","c_0","c_es1","c_es2","c_es3")
  print(sim.true)
  write.csv(sim.true,file=file.path(OUT,paste0('True_Case',Case,'_Sigma',Sigma,'_Scale',Scale,'_Shape',p.shape,'.csv')))
}
  
 }
# By default calculate all sixteen historical settings. A job ID selects one.
# Any scenario override selects branch 1 unless REINSURANCE_JOB_ID is supplied.
scenario_variables <- c("REINSURANCE_CASE", "REINSURANCE_SIGMA",
                        "REINSURANCE_SCALE", "REINSURANCE_SHAPE")
selected_job <- Sys.getenv("REINSURANCE_JOB_ID", "")
if (nzchar(selected_job)) {
  jobs <- env_integer("REINSURANCE_JOB_ID", 1L)
} else if (any(nzchar(Sys.getenv(scenario_variables, "")))) {
  jobs <- 1L
} else {
  jobs <- seq_len(16L)
}
if (any(!jobs %in% 1:16)) stop("REINSURANCE_JOB_ID must be in 1:16 for Simu_True.R.")
for (job_id in jobs) parallel_func(job_id)

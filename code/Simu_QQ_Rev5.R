# New Q-Q-only version of the supplied simulation script.
# Estimation formulas below are preserved from Code/Simu_True.R.
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
attr(c_vec_hat, "solver_diagnostics") <- if (inherits(sol, "try-error")) {
  c(termcd=NA_real_, max_residual=Inf)
} else {
  c(termcd=sol$termcd, max_residual=max(abs(sol$fvec)))
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

attr(res, "solver_diagnostics") <- attr(c_vec_hat, "solver_diagnostics")
return(res)
} 


 
# Q-Q diagnostics only. The original estimator, double-pt generator, and mu0
# are deliberately retained for traceability to the historical implementation.
# See PORTABILITY_NOTES.md before interpreting this generator as the stated Pareto model.
# Run from any directory, for example:
# Rscript code/Simu_QQ_Rev5.R --n=3000 --reps=2000 --workers=4
# Reuse saved draws to redraw: Rscript code/Simu_QQ_Rev5.R --plot-only
# A larger run is configurable, e.g. --n=100000 --reps=10000 --out=QQ_Rev5_large.

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value=TRUE)
  if (!length(hit)) default else sub(paste0("^--",key,"="),"",tail(hit,1))
}
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1]))
project_root <- dirname(dirname(script_path))
cfg <- list(n=as.integer(get_arg("n",Sys.getenv("REINSURANCE_N","3000"))),
            reps=as.integer(get_arg("reps",Sys.getenv("REINSURANCE_REPS","2000"))),
            seed=as.integer(get_arg("seed",Sys.getenv("REINSURANCE_SEED","20261021"))),
            cases=as.integer(strsplit(get_arg("cases",Sys.getenv("REINSURANCE_CASES",Sys.getenv("REINSURANCE_CASE","1,2,3,4"))),",")[[1]]),
            shapes=as.numeric(strsplit(get_arg("shapes",Sys.getenv("REINSURANCE_SHAPES",Sys.getenv("REINSURANCE_SHAPE","3.5,6"))),",")[[1]]),
            scale=as.integer(get_arg("scale",Sys.getenv("REINSURANCE_SCALE","1"))),
            mode="original_implementation")
workers <- as.integer(get_arg("workers",Sys.getenv("REINSURANCE_WORKERS",Sys.getenv("REINSURANCE_NCORE","4"))))
stopifnot(cfg$n>=100, cfg$n%%2==0, cfg$reps>=20,
          workers>=1,workers<=8,all(cfg$cases%in%1:4),cfg$scale%in%1:2)
out_arg <- get_arg("out",Sys.getenv("REINSURANCE_OUT","outputs/qq"))
out <- if (grepl("^(/|[A-Za-z]:[/\\])",out_arg)) out_arg else file.path(project_root,out_arg)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
raw_dir <- file.path(out,"replicates")
dir.create(raw_dir,recursive=TRUE,showWarnings=FALSE)
config_path <- file.path(out,"configuration.rds")
if (file.exists(config_path)) {
  if (!identical(readRDS(config_path),cfg)) stop("Output directory has a different configuration; choose a new --out.")
} else saveRDS(cfg,config_path)
plot_only <- "--plot-only"%in%args

for (pkg in c("QRM","nleqslv"))
  if (!requireNamespace(pkg,quietly=TRUE)) stop("Missing package: ",pkg)
library(QRM)
library(nleqslv)

generate_replicate <- function(seed,case,shape) {
  set.seed(seed)
  scales <- if(cfg$scale==1)c(1000,1000,1000) else c(800,1000,1200)
  R <- matrix(.25,3,3);diag(R)<-1
  # Exact original generator. QRM::rcopula.t already returns uniforms;
  # the second pt is retained, NOT silently corrected in a plotting-only change.
  U <- QRM::rcopula.t(cfg$n,df=3,R)
  U_unif <- pt(U,df=3)
  X <- sweep((1-U_unif)^(-1/shape)-1,2,scales,"*")
  beta_vec <- rep(.05,3)
  tilde_beta_vec <- switch(case,
    rep(.1,3),
    c(.1+.01*scales[2]/scales[1],.1,.1+.01*scales[2]/scales[3]),
    c(.1-.01*scales[2]/scales[1],.1,.1-.01*scales[2]/scales[3]),
    c(.1+.01*scales[2]/scales[1],.1,.1-.01*scales[2]/scales[3]))
  mu0 <- .2*.1/(shape-1)*sum(scales) # beta_scalar=0, as in the original.
  warns <- character()
  result <- withCallingHandlers(
    estimate_theta(X,0,.1,beta_vec,tilde_beta_vec,mu0,.95,rep(1,cfg$n)),
    warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
  list(estimates=as.numeric(result),
       solver=attr(result,"solver_diagnostics"),warnings=unique(warns))
}

scenario_id <- function(case,shape)paste0("Case",case,"_Scale",cfg$scale,"_Shape",shape)
all_draws <- list()
audit <- list()
for(shape in cfg$shapes)for(case in cfg$cases) {
  id <- scenario_id(case,shape)
  checkpoint <- file.path(raw_dir,paste0(id,".rds"))
  draws <- if(file.exists(checkpoint))readRDS(checkpoint) else list()
  if(length(draws)<cfg$reps && plot_only)stop("Incomplete saved replicates: ",id)
  if(length(draws)<cfg$reps) {
    message(id,": ",length(draws),"/",cfg$reps)
    starts <- seq.int(length(draws)+1,cfg$reps,by=100L)
    for(start in starts) {
      ids <- seq.int(start,min(start+99L,cfg$reps))
      # Independent reproducible seeds; the same streams across scenarios
      # implement common random numbers, as did the original per-case seed reset.
      work <- function(i)generate_replicate(cfg$seed+i,case,shape)
      chunk <- if(workers==1)lapply(ids,work) else
        if (.Platform$OS.type == "windows") {
          message("Windows: using one worker for portable execution.")
          lapply(ids,work)
        } else parallel::mclapply(ids,work,mc.cores=workers,mc.set.seed=FALSE)
      if(any(vapply(chunk,inherits,logical(1),"try-error")))stop("A worker failed; saved draws remain resumable.")
      draws <- c(draws,chunk)
      saveRDS(draws,checkpoint)
      message(id,": ",length(draws),"/",cfg$reps)
    }
  }
  mat <- do.call(rbind,lapply(draws,`[[`,"estimates"))
  colnames(mat)<-c("theta1","theta2","theta3","rho_J","rho_L_J","rho_I","rho_U_I","c_joint","c1","c2","c3")
  if(any(!is.finite(mat)))stop("Non-finite estimates: no observations will be silently removed: ",id)
  solver <- do.call(rbind,lapply(draws,`[[`,"solver"))
  all_draws[[id]]<-mat
  write.csv(mat,file.path(raw_dir,paste0(id,".csv")),row.names=FALSE)
  audit[[id]]<-data.frame(scenario=id,n=cfg$n,replications=nrow(mat),
    max_residual=max(solver[,"max_residual"]),
    residual_gt_1e_4=sum(solver[,"max_residual"]>1e-4),
    solver_not_code_1=sum(is.na(solver[,"termcd"])|solver[,"termcd"]!=1),
    warning_replications=sum(lengths(lapply(draws,`[[`,"warnings"))>0))
}
write.csv(do.call(rbind,audit),file.path(out,"solver_audit.csv"),row.names=FALSE)
saveRDS(list(configuration=cfg,draws=all_draws,session=sessionInfo()),file.path(out,"qq_results.rds"))

qq_panel <- function(values,estimator,case,shape,zoom=FALSE) {
  # No winsorization, tail trimming, selective exclusions or normality p-values.
  z <- sort((values-mean(values))/sd(values))
  p <- ppoints(length(z));x<-qnorm(p)
  keep <- if(zoom)p>=.90 else rep(TRUE,length(p))
  cols <- ifelse(p>=.99,"#D55E00",ifelse(p>=.95,"#0072B2","#555555"))
  xlim <- if(zoom)c(qnorm(.90),max(x)) else range(x)
  ylim <- range(c(z[keep],x[keep]))
  plot(x[keep],z[keep],type="n",xlim=xlim,ylim=ylim,
       xlab="Standard normal quantile",ylab="Standardized estimate quantile",
       main="",
       cex.main=1,cex.lab=.86,cex.axis=.82,bty="l")
  abline(a=0,b=1,col="#222222",lwd=1.1)
  abline(v=qnorm(c(.95,.99)),col="#b0b0b0",lty=3,lwd=.65)
  points(x[keep],z[keep],pch=16,cex=.35,col=cols[keep])
}

# Write each Q-Q plot to its own PDF; no title or legend is drawn.
draw_single_plot <- function(label, case, shape, estimator, zoom=FALSE) {
  filename <- paste0("QQ_",label,"_Case",case,"_Est",estimator,".pdf")
  pdf(file.path(out,filename),width=4,height=2.7,
      family="Times",useDingbats=FALSE)
  par(mar=c(3.5,3.8,1.2,.6),oma=rep(0,4),family="serif",
      mgp=c(2.2,.65,0),tcl=-.25)
  qq_panel(all_draws[[scenario_id(case,shape)]][,estimator],
           estimator,case,shape,zoom)
  dev.off()
}

groups <- list(Main=list(cases=1:2,shape=3.5),
               Supp_A=list(cases=3:4,shape=3.5),
               Supp_B=list(cases=1:2,shape=6),
               Supp_C=list(cases=3:4,shape=6))
available <- function(g)all(g$cases%in%cfg$cases)&&g$shape%in%cfg$shapes
for(label in names(groups))if(available(groups[[label]])) {
  g <- groups[[label]]
  for(case in g$cases)for(estimator in 1:3)
    draw_single_plot(label,case,g$shape,estimator)
}
if(available(groups$Main))for(case in 1:2)for(estimator in 1:3)
  draw_single_plot("Main_UpperTail",case,3.5,estimator,zoom=TRUE)

summary_rows <- list()
for(id in names(all_draws))for(j in 1:3) {
  y<-all_draws[[id]][,j];z<-(y-mean(y))/sd(y)
  summary_rows[[paste0(id,"_",j)]]<-data.frame(scenario=id,estimator=j,
    mean=mean(y),sd=sd(y),skewness=mean(z^3),
    q95=unname(quantile(z,.95)),normal_q95=qnorm(.95),
    q99=unname(quantile(z,.99)),normal_q99=qnorm(.99))
}
write.csv(do.call(rbind,summary_rows),file.path(out,"qq_summary.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
message("Completed Q-Q figures and saved all replicates in ",out)

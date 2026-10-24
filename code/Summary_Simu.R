# Portable public copy; estimation and historical random generation are unchanged.
# Run from the repository root with Rscript code/Summary_Simu.R.
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

myn <- env_integer("REINSURANCE_N", 3000L, 4L)
mycase <- env_integer("REINSURANCE_CASE", 4L)
mysigma <- env_integer("REINSURANCE_SIGMA", 1L)
myscale <- env_integer("REINSURANCE_SCALE", 2L)
myshape <- env_number("REINSURANCE_SHAPE", 6)
myboot <- env_integer("REINSURANCE_BOOT", 2000L, 40L)
if (myboot %% 40L != 0L) stop("REINSURANCE_BOOT must be a multiple of 40 for the original interval order statistics.")
if (!mycase %in% 1:4 || !mysigma %in% 1:2 || !myscale %in% 1:2 || myshape <= 1)
  stop("Invalid simulation scenario.")
resolve_directory <- function(name, default) {
  value <- Sys.getenv(name, default)
  if (!grepl("^(/|[A-Za-z]:[/\\])", value)) value <- file.path(repository_root, value)
  value
}
simulation_directory <- resolve_directory("REINSURANCE_SIM_DIR", "outputs/simulation")
truth_directory <- resolve_directory("REINSURANCE_TRUE_DIR", "outputs/true")
OUT <- output_directory("outputs/summary")

Myname=paste0('Simu_n',myn,'_Case',mycase,'_Sigma',mysigma,'_Scale',myscale,'_Shape',myshape);
files <- list.files(simulation_directory, pattern = paste0("^", Myname, "_.*\\.Rdata$"), full.names = TRUE)
if (!length(files)) stop("No matching simulation outputs in ", simulation_directory, "; run Simu.R for this scenario first.")

#Myname.TV=paste0('True_Case',mycase,'_Sigma',mysigma,'_Scale',myscale,'_Shape',myshape,'.csv');

Myname.TV <- file.path(truth_directory,paste0('True_Case',mycase,'_Sigma',mysigma,'_Scale',myscale,'_Shape',myshape,'.csv'))
if (!file.exists(Myname.TV)) stop("Missing reference CSV: ", Myname.TV, "; run Simu_True.R for this scenario first.")
TrueValue=read.csv(Myname.TV)[,-1]
print(TrueValue)
Myres=NULL;

for(myi in 1:length(files)){
info=file.info(files[myi]);
if(info$size!=0){
loaded <- new.env(parent = emptyenv())
load(files[myi], envir = loaded)
if (!exists("res", envir = loaded, inherits = FALSE)) stop("Missing object res in ", files[myi])
res <- loaded$res
if (!is.matrix(res) || !is.numeric(res) || ncol(res) != 11L * (myboot + 1L))
  stop("Unexpected simulation dimensions in ", files[myi], "; REINSURANCE_BOOT must match the generated outputs.")
Myres=rbind(res,Myres);
}
}
if (is.null(Myres) || !nrow(Myres)) stop("No nonempty simulation outputs were loaded.")
Myres=na.omit(Myres)
if (!nrow(Myres)) stop("No complete simulation rows remain.")
if (any(!is.finite(Myres))) stop("Simulation outputs contain non-finite values.")
mydim=dim(Myres)[1]
limit <- env_integer("REINSURANCE_SUMMARY_LIMIT", 1000L)
if(mydim>limit){
Myres=Myres[seq_len(limit),,drop=FALSE]
}
mydim=dim(Myres)[1]
cover=matrix(NA,nrow=mydim,ncol=6)

for(myi in 1:mydim){
tmp1=rep(NA,myboot)
tmp2=tmp1
tmp3=tmp1
for(myb in 1:myboot){
tmp1[myb]=Myres[myi,1+11*myb]-Myres[myi,1]
tmp2[myb]=Myres[myi,2+11*myb]-Myres[myi,2]
tmp3[myb]=Myres[myi,3+11*myb]-Myres[myi,3]
}
tmp1.s=sort(tmp1)
tmp2.s=sort(tmp2)
tmp3.s=sort(tmp3)
tmp0.1=Myres[myi,1]-TrueValue[1]
tmp0.2=Myres[myi,2]-TrueValue[1]
tmp0.3=Myres[myi,3]-TrueValue[1]
cover[myi,1]=(tmp0.1>tmp1.s[myboot*0.05])*(tmp0.1<tmp1.s[myboot*0.95])
cover[myi,2]=(tmp0.2>tmp2.s[myboot*0.05])*(tmp0.2<tmp2.s[myboot*0.95])
cover[myi,3]=(tmp0.3>tmp3.s[myboot*0.05])*(tmp0.3<tmp3.s[myboot*0.95])
cover[myi,4]=(tmp0.1>tmp1.s[myboot*0.025])*(tmp0.1<tmp1.s[myboot*0.975])
cover[myi,5]=(tmp0.2>tmp2.s[myboot*0.025])*(tmp0.2<tmp2.s[myboot*0.975])
cover[myi,6]=(tmp0.3>tmp3.s[myboot*0.025])*(tmp0.3<tmp3.s[myboot*0.975])
}
print(mydim)
print(apply(cover,2,mean))

coverage_summary <- data.frame(
  estimator = rep(c("theta1", "theta2", "theta3"), 2L),
  nominal_level = rep(c(0.90, 0.95), each = 3L),
  coverage = colMeans(cover), replications = mydim, bootstrap_replicates = myboot)
write.csv(coverage_summary, file.path(OUT, paste0(Myname, "_coverage.csv")), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, paste0(Myname, "_sessionInfo.txt")))

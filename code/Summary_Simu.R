
myn=3000;mycase=4;mysigma=1;myscale=2;myshape=6;
setwd("/Users/lpeng/2026/Quota_ASTIN_Rev2/Simu_Data")
Myname=paste0('Simu_n',myn,'_Case',mycase,'_Sigma',mysigma,'_Scale',myscale,'_Shape',myshape);
files=list.files(pattern =paste0(Myname,"_.*\\.Rdata"));

#Myname.TV=paste0('True_Case',mycase,'_Sigma',mysigma,'_Scale',myscale,'_Shape',myshape,'.csv');

Myname.TV=paste0('~/2026/Quota_ASTIN_Rev2/Simu_True/True_Case',mycase,'_Sigma',mysigma,'_Scale',myscale,'_Shape',myshape,'.csv');
TrueValue=read.csv(Myname.TV)[,-1]
print(TrueValue)
Myres=NULL;

for(myi in 1:length(files)){
info=file.info(files[myi]);
if(info$size!=0){
load(files[myi]);
Myres=rbind(res,Myres);
}
}
Myres=na.omit(Myres)
mydim=dim(Myres)[1]
if(mydim>1000){
Myres=Myres[1:1000,]
}
mydim=dim(Myres)[1]
cover=matrix(NA,nrow=mydim,ncol=6)
myboot=2000
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
print(apply(cover,2,mean)-c(0.9,0.9,0.9,0.95,0.95,0.95))

setwd("/sc/arion/scratch/wangm08/liver_gtex/outdir/rnavalid_chr/")
i=1
resRaw=resRaw+get(paste("res_chr",i,sep=""))$resRaw
resFilter1=resFilter1+get(paste("res_chr",i,sep=""))$resFilter1
resFilter2=resFilter2+get(paste("res_chr",i,sep=""))$resFilter2

for(i in 2:20){
 file=paste("chr",i,".rda",sep="")
 load(file)
 resRaw=resRaw+get(paste("res_chr",i,sep=""))$resRaw
 resFilter1=resFilter1+get(paste("res_chr",i,sep=""))$resFilter1
 resFilter2=resFilter2+get(paste("res_chr",i,sep=""))$resFilter2
 }
i=21
res_chr21=chr21_res
 resRaw=resRaw+get(paste("res_chr",i,sep=""))$resRaw
 resFilter1=resFilter1+get(paste("res_chr",i,sep=""))$resFilter1
 resFilter2=resFilter2+get(paste("res_chr",i,sep=""))$resFilter2
 i=22
 file=paste("chr",i,".rda",sep="")
 load(file)
 resRaw=resRaw+get(paste("res_chr",i,sep=""))$resRaw
 resFilter1=resFilter1+get(paste("res_chr",i,sep=""))$resFilter1
 resFilter2=resFilter2+get(paste("res_chr",i,sep=""))$resFilter2

rbind(resRaw[,1:5],resFilter1[,1:5],resFilter2[,1:5])

rbind(rbind(apply(tmpRaw,1,mean),apply(tmpRaw,1,sd)),rbind(apply(tmpFilter1,1,mean),apply(tmpFilter1,1,sd)),rbind(apply(tmpFilter2,1,mean),apply(tmpFilter2,1,sd)))

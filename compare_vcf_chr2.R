#chr21<-read.table("/sc/arion/scratch/wangm08/liver_gtex/outdir/liftover_chr/liver_all_combine.chr21.liftover.vcf",as.is=T)
rnaValid<-function(chr21){
#dim(chr21)
#chr21[1:5,1:%]
#chr21[1:5,1:5]
#chr21[1:5,1:10]
cmd = "grep \"#CHROM\"  /sc/arion/scratch/wangm08/liver_gtex/outdir/liftover_chr/liver_all_combine.chr21.liftover.vcf"
#header <- sysmtem(cmd,intern=T)
header <- system(cmd,intern=T)
#header
#strsplit(header,"\t")[[1]]
samp<-strsplit(header,"\t")[[1]][10:147]
#length(samp)
#dim(chr21)
#dp=matrix(NA,nrow=nrow(chr21),ncol=ncol(chr21))
#colnames(dp)=samp
#col.names(dp)=samp
#colnames(dp) <- samp
#dim(dp)
#length(samp)
dp=matrix(NA,nrow=nrow(chr21),ncol=length(samp))
colnames(dp)=samp
#chr21[1:5,1:5]
tmp=paste(chr21[,1],chr21[,2],chr21[,4],chr21[,5],sep="_")
#length(tmp)
#tmp[1]
#row.names(dp)=samp
row.names(dp)<-tmp
gt<-dp
vr<-dp
nr<-dp
#chr21[10:15,1:5]
#chr21[1:5,10:15]
#chr21[1:5,9:15]
#sapply(chr21[1:5,9],function(z) which(strsplit(z,':')[[1]]=="DP"))
#sapply(chr21[1:5,9],function(z) which(strsplit(z,':')[[1]]=="GT"))
#sapply(chr21[1:5,9],function(z) which(strsplit(z,':')[[1]]=="AD"))
#unique(sapply(chr21[,9],function(z) which(strsplit(z,':')[[1]]=="DP")))
#unique(sapply(chr21[,9],function(z) which(strsplit(z,':')[[1]]=="AD")))
#unique(sapply(chr21[,9],function(z) which(strsplit(z,':')[[1]]=="GT")))
#strsplit("./.:0,0:0:.:.",':')[[1]][1]
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],'/')[[1]]
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],'\/')[[1]]
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],['/'])[[1]]
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],"\/")[[1]]
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],"[||/]")[[1]]
#strsplit("./.:0,0:0:.:.",':')[[1]][1]
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],"[||/]")
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],split="[||/]")
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],"[||/]")
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],'[||/]')
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],'||/')
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],[||/])
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],['||/'])
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],'||/')
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],'|/')
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],'[|/]')
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],'[\|/]')
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],'[||/]')
#split(strsplit("./.:0,0:0:.:.",':')[[1]][1],"[||/]")
#strsplit(strsplit("./.:0,0:0:.:.",':')[[1]][1],"[||/]")
#as.numeric(strsplit(strsplit("./.:0,0:0:.:.",':')[[1]][1],"[||/]"))
#strsplit(strsplit("./.:0,0:0:.:.",':')[[1]][1],"[||/]")
#as.numeric(strsplit(strsplit("./.:0,0:0:.:.",':')[[1]][1],"[||/]"))
#as.numeric(strsplit(strsplit("./.:0,0:0:.:.",':')[[1]][1],"[||/]")[[1]])
#which(as.numeric(strsplit(strsplit("./.:0,0:0:.:.",':')[[1]][1],"[||/]")[[1]])>0)
#sum(as.numeric(strsplit(strsplit("./.:0,0:0:.:.",':')[[1]][1],"[||/]")[[1]]))
#sum(as.numeric(strsplit(strsplit("0/0:0,0:0:.:.",':')[[1]][1],"[||/]")[[1]]))
#sum(as.numeric(strsplit(strsplit("0/1:0,0:0:.:.",':')[[1]][1],"[||/]")[[1]]))
#gt[,1]=sapply(chr21[,10],sum(as.numeric(strsplit(strsplit(z,':')[[1]][1],"[||/]")[[1]])))
gt[,1]=sapply(chr21[,10],function(z) sum(as.numeric(strsplit(strsplit(z,':')[[1]][1],"[||/]")[[1]])))
#gt[1:5,1]
for(i in 2:ncol(gt)){
gt[,i]=sapply(chr21[,9+i],function(z) sum(as.numeric(strsplit(strsplit(z,':')[[1]][1],"[||/]")[[1]])))
}
#dim(gt)
#gt[1:5,1:5]
#length(which(gt[,1]==1))
#length(which(gt[,1]==2))
#length(which(gt[,1]==0))
#chr21[1:5,10:15]
#chr21[1:5,9:15]
#strsplit("./.:0,0:0:.:.",':')[[1]][3]
for(i in 1:ncol(dp)){
dp[,i]=sapply(chr21[,9+i],function(z) strsplit(z,':')[[1]][3])
}
#dp[1:5,1:5]
#strsplit(strsplit("./.:0,0:0:.:.",':')[[1]][2],",")[[1]]
#strsplit(strsplit("./.:0,0:0:.:.",':')[[1]][2],",")[[1]][1]
for(i in 1:ncol(vr)){
vr[,i]=sapply(chr21[,9+i],function(z) strsplit(strsplit(z,':')[[1]][2],",")[[1]][2])
}
#for(i in 1:ncol(vt)){
#vt[,i]=sapply(chr21[,9+i],function(z) strsplit(strsplit(z,':')[[1]][2],",")[[1]][2])
#}
#dim(vr)
#ls()
#nr=vr
#vr=NA
#vr[1:5,1:%]
#vr[1:5,1:5]
#vr=matrix(NA,nrow=nrow(chr21),ncol=length(samp))
#row.names(vr)<-tmp
#colnames(vr)=samp
for(i in 1:ncol(nr)){
nr[,i]=sapply(chr21[,9+i],function(z) strsplit(strsplit(z,':')[[1]][2],",")[[1]][1])
}
#savehistory("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_vcf_chr.Rhistory")

tmpr=row.names(dp)
tmpc=colnames(dp)
dp<-matrix(as.numeric(dp),nrow=nrow(dp))
#dp[1:5,1:5]
#length(which(is.na(dp[,1])))
#length(which(is.na(dp)))
row.names(dp)=tmpr
colnames(dp)=tmpc
nr<-matrix(as.numeric(nr),nrow=nrow(nr))
colnames(nr)=colnames(gt)
row.names(nr)=row.names(gt)
vr<-matrix(as.numeric(vr),nrow=nrow(vr))
colnames(vr)=colnames(gt)
row.names(vr)=row.names(gt)
vcfr=vr/(vr+nr)
##############################################################################################################################################

mat<-read.table("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/genotype/rna_geno_match_samp",as.is=T,sep="\t")
#mat[1:5,]
ind=match(mat[,1],colnames(gt))
gtx=gt[,ind]
#dim(gtx)
chr21geno<-read.table("/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir/bychr/chr21.geno.vcf",as.is=T)
#dim(chr21geno)
cmd = "grep \"#CHROM\"  /sc/arion/scratch/wangm08/liver_gtex/genotype/outdir/bychr/chr21.geno.vcf"
header_geno <- system(cmd,intern=T)
#header_geno
samp_geno<-strsplit(header_geno,"\t")[[1]][10:125]
#length(samp_geno)
#length(which(samp_geno!=mat[,2]))
#chr21geno[1:5,9:15]
gtg = matrix(NA,nrow=nrow(chr21geno),ncol=length(samp_geno))
colnames(gtg)=samp_geno
#chr21geno[1:5,1:5]
chr21geno[,1]=paste("chr",chr21geno[,1],sep="")
#chr21geno[1:5,1:5]
tmpx = paste(chr21geno[,1],chr21geno[,2],chr21geno[,4],chr21geno[,5],sep="_")
#tmpx[1:5]
row.names(gtg) = tmpx
#length(tmp)
#length(intersect(tmp,tmpx))
#14100/19920
#length(tmpx)
#gtg=gt
#gt = matrix(NA,nrow=nrow(chr21),ncol=length(samp))
#colnames(gt)=samp
#row.names(gt)<-tmp
#gt[,1]=sapply(chr21[,10],function(z) sum(as.numeric(strsplit(strsplit(z,':')[[1]][1],"[||/]")[[1]])))
#for(i in 2:ncol(gt)){
#gt[,i]=sapply(chr21[,9+i],function(z) sum(as.numeric(strsplit(strsplit(z,':')[[1]][1],"[||/]")[[1]])))
#}
#gt[1:5,1:5]
#gtg[1:5,1:%]
#gtg[1:5,1:5]
#chr21geo[1:5,9:15]
#chr21geno[1:5,9:15]
#strsplit("0|0",'|')
#strsplit("0|0",'\|')
#strsplit("0|0",'|||')
#strsplit("0|0",'[|||]')
for(i in 1:ncol(gtg)){
gtg[,i]=sapply(chr21geno[,9+i],function(z) sum(as.numeric(strsplit(z,"[|||]")[[1]])))
}
#dim(gat)
#dim(gtg)
#gtg[1:5,1:5]
#gt[1:5,1:%]
#gt[1:5,1:5]
#rowSums(gt[1:5,])
#rowSums(gt[1:5,],na.rm=TRUE)
#apply(gt[1:5,],1,function(z) length(which(z>0)))
#which(gt[1,]>0)
#length(which(apply(gt,1,function(z) length(which(z>0)))>1))






dpx=dp[,ind]
vcfrx=vcfr[,ind]
indx=which(apply(gtx,1,function(z) length(which(z>0)))>0)

gtx=gtx[indx,]
dpx=dpx[indx,]
vcfrx=vcfrx[indx,]


resr=matrix(0,nrow=6,ncol=ncol(gtg))
resr[1,]=sapply(c(1:ncol(gtx)),function(z) length(which(gtx[,z]>0)))
resr[2,]=sapply(c(1:ncol(gtx)),function(z) length(which(gtx[,z]==1)))
resr[3,]=sapply(c(1:ncol(gtx)),function(z) length(which(gtx[,z]==2)))
resr[4,]=sapply(c(1:ncol(gtx)),function(z) length(intersect(row.names(gtx)[which(gtx[,z]>0)],row.names(gtg)[which(gtg[,z]>0)])))
resr[5,]=sapply(c(1:ncol(gtx)),function(z) length(intersect(row.names(gtx)[which(gtx[,z]==1)],row.names(gtg)[which(gtg[,z]==1)])))
resr[6,]=sapply(c(1:ncol(gtx)),function(z) length(intersect(row.names(gtx)[which(gtx[,z]==2)],row.names(gtg)[which(gtg[,z]==2)])))
colnames(resr)=colnames(gtg)
row.names(resr)=c("predict","predit_het","predit_homo","positive","het","homo")



#8308/10413
#dp[1:5,1:5]
#dp < data.matrix(dp)
#dp <- data.matrix(dp)
#dp[1:5,1:5]
#tmpr=row.names(dp)
#tmpc=colnames(dp)
#dp<-matrix(as.numeric(dp),nrow=nrow(dp))
#dp[1:5,1:5]
#length(which(is.na(dp[,1])))
#length(which(is.na(dp)))
#row.names(dp)=tmpr
#colnames(dp)=tmpc
#dp[1:5,1:5]
#max(dp)
#max(dp,na.rm=T)
#nr[1:5,1:5]
#length(which(is.na(nr)))
#max(nr)
#min(nr)
#max(vr)
#min(vr)
#max(as.numeric(vr))
#max(as.numeric(nr))
#nr<-matrix(as.numeric(nr),nrow=nrow(nr))
#colnames(nr)=colnames(gt)
#row.names(nr)=row.names(gt)
#vr<-matrix(as.numeric(vr),nrow=nrow(vr))
#colnames(vr)=colnames(gt)
#row.names(vr)=row.names(gt)
#vcfr=vr/(vr+nr)
#vcfr[1:5,1:%]
#vcfr[1:5,1:5]
tmp1=dpx
tmp1[which(dpx<10)]=NA
tmp2=vcfrx
tmp2[which(vcfrx<0.2)]=NA
#tmp2[1:5,1:%]
#tmp2[1:5,1:5]
tmp2[which(is.nan(tmp2))]=NA
#tmp2[1:5,1:%]
#tmp2[1:5,1:5]
tmp3=tmp1*tmp2


#tmp3[which(!is.na(tmp3))]=1
#length(which(tmp3[,1]==1))
#unique(tmp3)
#unique(unlist(tmp3))
gt1=gtx*tmp3
#dim(gtx)
#gtx[1:5,1:%]
#gtx[1:5,1:5]
#ls()
#dim(gtg)
#dim(gt)
#gt1=gt[,ind]
#gt1[1:5,1:%]
#gt1[1:5,1:5]
#length(which(colnames(gt1)!=mat[,1]))
#apply(gt1[1:5,],1,function(z) length(which(!is.na(z))))
#ind1=which(apply(gt1[1:5,],1,function(z) length(which(!is.na(z))))>0)
#length(ind1)
ind1=which(apply(gt1[,],1,function(z) length(which(!is.na(z))))>0)
#length(ind1)
#dim(gt`)
#dim(gt1)
#dim(gt1)
#dim(tmp3)
#tmp3=tmp3[,ind]
#gt1=gt1*tmp3
#ind1=which(apply(gt1[,],1,function(z) length(which(!is.na(z))))>0)
#length(ind1)
#length(intersect(row.names(gt1)[ind1],row.names(gtg)))
#4040/4459
#dim(gt1)
gt1=gt1[ind1,]
#dim(gt2)
#gt2[1:5,1:%]
#gt2[1:5,1:5]
#row.names(gt2)[which(gt2[,1]>0)][1:5]
#length(intersect(row.names(gt2)[which(gt2[,1]>0)],row.names(gtg)[which(gtg[,1]>0)]))
#length(intersect(row.names(gt2)[which(gt2[,1]>0)],row.names(gtg)[which(gtg[,1]>0)]))/length(which(gt2[,1]>0))
#savehistory("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_vcf_chr2.Rhistory")


resr1=matrix(0,nrow=6,ncol=ncol(gtg))
resr1[1,]=sapply(c(1:ncol(gt1)),function(z) length(which(gt1[,z]>0)))
resr1[2,]=sapply(c(1:ncol(gt1)),function(z) length(which(gt1[,z]==1)))
resr1[3,]=sapply(c(1:ncol(gt1)),function(z) length(which(gt1[,z]==2)))
resr1[4,]=sapply(c(1:ncol(gt1)),function(z) length(intersect(row.names(gt1)[which(gt1[,z]>0)],row.names(gtg)[which(gtg[,z]>0)])))
resr1[5,]=sapply(c(1:ncol(gt1)),function(z) length(intersect(row.names(gt1)[which(gt1[,z]==1)],row.names(gtg)[which(gtg[,z]==1)])))
resr1[6,]=sapply(c(1:ncol(gt1)),function(z) length(intersect(row.names(gt1)[which(gt1[,z]==2)],row.names(gtg)[which(gtg[,z]==2)])))
colnames(resr1)=colnames(gtg)
row.names(resr1)=c("predict","predit_het","predit_homo","positive","het","homo")

ind2=which(apply(gt1,1,function(z) length(which(!is.na(z))))>1)
gt2=gt1[ind2,]
resr2=matrix(0,nrow=6,ncol=ncol(gtg))
resr2[1,]=sapply(c(1:ncol(gt2)),function(z) length(which(gt2[,z]>0)))
resr2[2,]=sapply(c(1:ncol(gt2)),function(z) length(which(gt2[,z]==1)))
resr2[3,]=sapply(c(1:ncol(gt2)),function(z) length(which(gt2[,z]==2)))
resr2[4,]=sapply(c(1:ncol(gt2)),function(z) length(intersect(row.names(gt2)[which(gt2[,z]>0)],row.names(gtg)[which(gtg[,z]>0)])))
resr2[5,]=sapply(c(1:ncol(gt2)),function(z) length(intersect(row.names(gt2)[which(gt2[,z]==1)],row.names(gtg)[which(gtg[,z]==1)])))
resr2[6,]=sapply(c(1:ncol(gt2)),function(z) length(intersect(row.names(gt2)[which(gt2[,z]==2)],row.names(gtg)[which(gtg[,z]==2)])))
colnames(resr2)=colnames(gtg)
row.names(resr2)=c("predict","predit_het","predit_homo","positive","het","homo")


res<-list(npredict=nrow(gtx),ngeno=nrow(gtg),nprepdit_pos=length(intersect(row.names(gtx),row.names(gtg))),resRaw=resr,npredict_filter1=nrow(gt1),nfilter1_pos=length(intersect(row.names(gt1),row.names(gtg))),resFilter1=resr1,npredict_filter2=nrow(gt2),nfilter2_pos=length(intersect(row.names(gt2),row.names(gtg))),resFilter2=resr2)ength(intersect(row.names(gt)[indx],row.names(gtg)))

return(res)

}

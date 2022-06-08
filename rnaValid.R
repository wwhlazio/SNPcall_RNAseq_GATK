#chr21<-read.table("/sc/arion/scratch/wangm08/liver_gtex/outdir/liftover_chr/liver_all_combine.chr21.liftover.vcf",as.is=T)
#chr21geno<-read.table("/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir/bychr/chr21.geno.vcf",as.is=T)
rnaValid<-function(chr21_file,chr21geno_file){
chr21<-read.table(chr21_file,as.is=T)
chr21geno<-read.table(chr21geno_file,as.is=T)
cmd = paste("grep \"#CHROM\"  ",chr21_file,sep="")
header <- system(cmd,intern=T)
samp<-strsplit(header,"\t")[[1]][10:147]
dp=matrix(NA,nrow=nrow(chr21),ncol=length(samp))
colnames(dp)=samp
tmp=paste(chr21[,1],chr21[,2],chr21[,4],chr21[,5],sep="_")
row.names(dp)<-tmp
gt<-dp
vr<-dp
nr<-dp

gt[,1]=sapply(chr21[,10],function(z) sum(as.numeric(strsplit(strsplit(z,':')[[1]][1],"[||/]")[[1]])))
for(i in 2:ncol(gt)){
gt[,i]=sapply(chr21[,9+i],function(z) sum(as.numeric(strsplit(strsplit(z,':')[[1]][1],"[||/]")[[1]])))
}

for(i in 1:ncol(dp)){
dp[,i]=sapply(chr21[,9+i],function(z) strsplit(z,':')[[1]][3])
}

for(i in 1:ncol(vr)){
vr[,i]=sapply(chr21[,9+i],function(z) strsplit(strsplit(z,':')[[1]][2],",")[[1]][2])
}

for(i in 1:ncol(nr)){
nr[,i]=sapply(chr21[,9+i],function(z) strsplit(strsplit(z,':')[[1]][2],",")[[1]][1])
}

tmpr=row.names(dp)
tmpc=colnames(dp)
dp<-matrix(as.numeric(dp),nrow=nrow(dp))
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
ind=match(mat[,1],colnames(gt))
gtx=gt[,ind]
cmd = paste("grep \"#CHROM\"  ",chr21geno_file,sep="")
header_geno <- system(cmd,intern=T)
samp_geno<-strsplit(header_geno,"\t")[[1]][10:125]
gtg = matrix(NA,nrow=nrow(chr21geno),ncol=length(samp_geno))
colnames(gtg)=samp_geno
chr21geno[,1]=paste("chr",chr21geno[,1],sep="")

tmpx = paste(chr21geno[,1],chr21geno[,2],chr21geno[,4],chr21geno[,5],sep="_")
row.names(gtg) = tmpx


for(i in 1:ncol(gtg)){
gtg[,i]=sapply(chr21geno[,9+i],function(z) sum(as.numeric(strsplit(z,"[|||]")[[1]])))
}







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



tmp1=dpx
tmp1[which(dpx<10)]=NA
tmp2=vcfrx
tmp2[which(vcfrx<0.2)]=NA
tmp2[which(is.nan(tmp2))]=NA
tmp3=tmp1*tmp2
tmp3[which(!is.na(tmp3))]=1


gt1=gtx*tmp3
ind1=which(apply(gt1[,],1,function(z) length(which(!is.na(z))))>0)
gt1=gt1[ind1,]


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


res<-list(npredict=nrow(gtx),ngeno=nrow(gtg),nprepdit_pos=length(intersect(row.names(gtx),row.names(gtg))),resRaw=resr,npredict_filter1=nrow(gt1),nfilter1_pos=length(intersect(row.names(gt1),row.names(gtg))),resFilter1=resr1,npredict_filter2=nrow(gt2),nfilter2_pos=length(intersect(row.names(gt2),row.names(gtg))),resFilter2=resr2)
return(res)

}

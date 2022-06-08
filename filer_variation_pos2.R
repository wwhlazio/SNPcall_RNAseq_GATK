chr21_file="/sc/arion/scratch/wangm08/liver_gtex/outdir/liftover_chr/liver_all_combine.chr21.liftover.vcf"
chr21geno_file="/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir/bychr/chr21.geno.vcf"
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
ind1=match(intersect(row.names(gtg),row.names(gtx)),row.names(gtx))
ind2=match(setdiff(row.names(gtx),row.names(gtg)),row.names(gtx))
count1=matrix(0,nrow=length(ind1),ncol=4)
row.names(count1)=row.names(gtx)[ind1]
for(i in 1:length(ind1)){
count1[i,match(names(table(gtx[ind1[i],])),c(0,1,2))]=table(gtx[ind1[i],])
count1[i,4]=116-sum(count1[i,1:3])
}
count2=matrix(0,nrow=length(ind2),ncol=4)
row.names(count2)=row.names(gtx)[ind2]
for(i in 1:length(ind2)){
count2[i,match(names(table(gtx[ind2[i],])),c(0,1,2))]=table(gtx[ind2[i],])
count2[i,4]=116-sum(count2[i,1:3])
}
count2x=round(count2*116)
count1x=round(count1*116)
length(which(rowSums(count1[,2:3])>0.05))
length(which(rowSums(count2[,2:3])>0.05))
count1[1:5,1:5]
count1[1:5,]
count2=count2/116
count1=count1/116
count2x=round(count2*116)
count1x=round(count1*116)
length(which(rowSums(count1[,2:3])>0.05))
length(which(rowSums(count2[,2:3])>0.05))
tmp1=dpx
tmp1[which(dpx<10)]=NA
tmp3=tmp1
tmp3[which(!is.na(tmp3))]=1
gt1=gtx*tmp3
ind1=which(apply(gt1[,],1,function(z) length(which(!is.na(z))))>0)
gt1=gt1[ind1,]
length(intersect(row.names(gt1),row.names(gtg)))
ind1=match(intersect(row.names(gtg),row.names(gt1)),row.names(gt1))
ind2=match(setdiff(row.names(gt1),row.names(gtg)),row.names(gt1))
count11=matrix(0,nrow=length(ind1),ncol=4)
for(i in 1:length(ind1)){
count11[i,match(names(table(gt1[ind1[i],])),c(0,1,2))]=table(gt1[ind1[i],])
count11[i,4]=116-sum(count11[i,1:3])
}
count22=matrix(0,nrow=length(ind2),ncol=4)
for(i in 1:length(ind2)){
count22[i,match(names(table(gt1[ind2[i],])),c(0,1,2))]=table(gt1[ind2[i],])
count22[i,4]=116-sum(count22[i,1:3])
}
count11x=count11/116
count22x=count22/116
length(which(rowSums(count22x[,2:3])>0.05))
length(which(rowSums(count11x[,2:3])>0.05))
count1[1:5,]
length(which(rowSums(count1[,2:3])>0.05))
length(which(rowSums(count1[,2:3])>0.05 & count1[,4]<=0.75))
length(which(rowSums(count2[,2:3])>0.05 & count2[,4]<=0.75))
341/(3050+341)
1-341/(3050+341)
length(which(rowSums(count1[,2:3])>0.05 & count1[,4]<=0.5))
length(which(rowSums(count2[,2:3])>0.05 & count2[,4]<=0.5))
214/(214+2130)
1-214/(214+2130)
3050+341
214+2130
tmp1=dpx
tmp1[which(dpx<10)]=NA
tmp3=tmp1
tmp3[which(!is.na(tmp3))]=1
gt1=gtx*tmp3
ind1=which(apply(gt1[,],1,function(z) length(which(!is.na(z))))>0)
gt1=gt1[ind1,]
dim(gt1)
dim(count11)
dim(count22)
length(which(rowSums(count11x[,2:3])>0.05))
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.25))
length(which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.25))
25/(461+25)
1-25/(461+25)
count11x[1:5,]
length(which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.75))
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.75))
64/(65+892)
1-64/(64+892)
65+892
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.5))
length(which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.5))
1-33/687
13483/13483
13483/18482
length(which(count11x[,4]<=0.75))
length(which(count22x[,4]<=0.75))
177(177+2607)
1-177/(177+2607)
177+2607
length(which(count22x[,4]<=0.5))
length(which(count11x[,4]<=0.5))
1874+85
1-85/(1874+85)
savehistory("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/filer_variation_pos2.Rhistory")

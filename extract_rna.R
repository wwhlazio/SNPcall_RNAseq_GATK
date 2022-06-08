extract_rna<-function(chr21_file){
ns=138
chr21<-read.table(chr21_file,as.is=T)
#chr21geno<-read.table(chr21geno_file,as.is=T)
cmd = paste("grep \"#CHROM\"  ",chr21_file,sep="")
header <- system(cmd,intern=T)
samp<-strsplit(header,"\t")[[1]][10:(ns+10-1)] 
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

tmpr=row.names(dp)
tmpc=colnames(dp)
dp<-matrix(as.numeric(dp),nrow=nrow(dp))
row.names(dp)=tmpr
colnames(dp)=tmpc



tmp1=dp
tmp1[which(dp<10)]=NA
tmp3=tmp1
tmp3[which(!is.na(tmp3))]=1
gt1=gt*tmp3
ind1=which(apply(gt1,1,function(z) length(which(!is.na(z))))>0)
gt1=gt1[ind1,]


count=matrix(0,nrow=length(ind1),ncol=4)
for(i in 1:length(ind1)){
count[i,match(names(table(gt1[i,])),c(0,1,2))]=table(gt1[i,])
count[i,4]=ns-sum(count[i,1:3])
}

count=count/ns

ind11 = which(rowSums(count[,2:3])>0.05 & count[,4]<=0.1)
return(gt1[ind11,])
}

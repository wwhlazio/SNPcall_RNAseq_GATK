setwd("/sc/arion/scratch/wangm08/liver_gtex/outdir/liftover/")
#files=list.files(pattern="\\.liftover.g.vcf$")
#files[1:2]
#files
#files=list.files(pattern="\\.liftover.g.vcf$")
a<-read.table("SRR1069141.liftover.g1.vcf",sep="\t",as.is=T)
tmp=cbind(a[[1]],a[[2]],a[[4]],a[[5]],a[[6]])
tmp=cbind(tmp,sapply(a[[8]],function(z) {ifelse(length(which(grepl("DP=",strsplit(z,";")[[1]])))>0, strsplit(strsplit(z,";")[[1]][which(grepl("DP=",strsplit(z,";")[[1]]))],"DP=")[[1]][2], NA)}))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) length(which(as.numeric(strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("GT",strsplit(a[[9]][z],":")[[1]]))],"[||/]")[[1]])>0))))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("AD",strsplit(a[[9]][z],":")[[1]]))],",")[[1]][1]))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) sum(as.integer(strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("AD",strsplit(a[[9]][z],":")[[1]]))],",")[[1]][-1]))))
colnames(tmp)=c("Chr","Pos","Ref","Alt","QUAL","DP","GT","REFR","ALTR")
rna=tmp
#dim(rna)
gt<-read.table("/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir/sample/GTEX-U8XE.vcf",as.is=T)
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/geno_num.R")
gt[[10]]<-as.character(geno_num(gt[[10]]))
gt[[1]]=paste("chr",gt[[1]],sep="")
tmp1=paste(gt[[1]],gt[[2]],sep="_")
tmp2=paste(rna[,1],rna[,2],sep="_")
#length(intersect(tmp1,tmp2))
#length(intersect(tmp1,tmp2))/length(tmp2)
tmp1x=paste(gt[[1]],gt[[2]],gt[[4]],gt[[10]],sep="_")
tmp2x=paste(rna[,1],rna[,2],rna[,3],rna[,7],sep="_")
#length(intersect(tmp1x,tmp2x))
#length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>10)]))
#length(which(as.numeric(tmp[,6])>10))
#25605/32837
#length(tmp2)
#length(intersect(tmp1,tmp2))
#length(intersect(tmp1x,tmp2x))
#67127/94818
#length(which(as.numeric(tmp[,6])>10))
#length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>10)]))
#cut=c(10,20,30,40)
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>10)]))/length(which(as.numeric(tmp[,6])>10)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>10)])),length(which(as.numeric(tmp[,6])>10)))
#res
#rec
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>20)]))/length(which(as.numeric(tmp[,6])>20)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>20)])),length(which(as.numeric(tmp[,6])>20)))
#rec
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>30)]))/length(which(as.numeric(tmp[,6])>30)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>30)])),length(which(as.numeric(tmp[,6])>30)))
#rec
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>40)]))/length(which(as.numeric(tmp[,6])>40)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>40)])),length(which(as.numeric(tmp[,6])>40)))
#rec
#c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>50)])),length(which(as.numeric(tmp[,5])>50)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>50)]))/length(which(as.numeric(tmp[,5])>50)))
#rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>50)])),length(which(as.numeric(tmp[,5])>50)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>50)]))/length(which(as.numeric(tmp[,5])>50)))
#rec
#tmp[1:5,]
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>50)])),length(which(as.numeric(tmp[,5])>70)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>70)]))/length(which(as.numeric(tmp[,5])>70)))
#rec
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>70)]))/length(which(as.numeric(tmp[,5])>70)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>70)])),length(which(as.numeric(tmp[,5])>70)))
#rec
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>100)]))/length(which(as.numeric(tmp[,5])>100)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>100)])),length(which(as.numeric(tmp[,5])>100)))
#rec
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>150)]))/length(which(as.numeric(tmp[,5])>150)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>150)])),length(which(as.numeric(tmp[,5])>150)))
#rec
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>200)]))/length(which(as.numeric(tmp[,5])>200)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>200)])),length(which(as.numeric(tmp[,5])>200)))
#rec
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.1/0.9)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.1/0.9)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.1/0.9)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.1/0.9)))
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.2/0.8)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.2/0.8)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.2/0.8)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.2/0.8)))
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.3/0.7)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.3/0.7)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.3/0.7)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.3/0.7)))
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.4/0.6)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.4/0.6)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.4/0.6)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.4/0.6)))
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.5/0.5)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.5/0.5)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.5/0.5)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>0.5/0.5)))
setwd("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/newpipe")
pdf("hist_VAF.pdf",height=9,width=16)
par(mfrow=c(1,2))
hist(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),8])),main="hist of VAF for True Positive",xlab="VAF")
hist(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),8])),main="hist of VAF for False Positive",xlab="VAF")
dev.off()
pdf("hist_QUAL.pdf",height=9,width=16)
par(mfrow=c(1,2))
hist(log10(as.numeric(tmp[match(intersect(tmp1,tmp2),tmp2),5])),main="hist of QUAL for True Positive",xlab="log10(QUAL)")
hist(log10(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),5])),main="hist of QUAL for False Positive",xlab="log10(QUAL)")
dev.off()
pdf("hist_dp.pdf",height=9,width=16)
par(mfrow=c(1,2))
hist(log10(as.numeric(tmp[match(intersect(tmp1,tmp2),tmp2),6])),main="hist of DP for True Positive",xlab="log10(DP)")
hist(log10(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),6])),main="hist of DP for False Positive",xlab="log10(DP)")
dev.off()




pdf("hist_dp_4part.pdf",height=9,width=16)
par(mfrow=c(2,2))
hist(log10(as.numeric(tmp[match(intersect(tmp1,tmp2),tmp2),6])),main="hist of DP for True Positive 67127",xlab="log10(DP)")
hist(log10(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),6])),main="hist of DP for False Positive 94818-67127",xlab="log10(DP)")
hist(log10(as.numeric(tmp[match(intersect(tmp1x,tmp2x),tmp2x),6])),main="hist of DP for consistent 60880",xlab="log10(DP)")
hist(log10(as.numeric(tmp[setdiff(match(intersect(tmp1,tmp2),tmp2),match(intersect(tmp1x,tmp2x),tmp2x)),6])),main="hist of DP for inconsistent 67127-60880",xlab="log10(DP)")
dev.off()



pdf("hist_VAF_4part.pdf",height=9,width=16)
par(mfrow=c(2,2))
hist(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),8])),main="hist of VAF for False Positive 94818-67127",xlab="VAF")
hist(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),8])),main="hist of VAF for True Positive 67127",xlab="VAF")
hist(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])/(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])+as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),8])),main="hist of VAF for consistent 60880",xlab="VAF")
hist(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])/(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])+as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),8])),main="hist of VAF for inconsistent 67127-60880",xlab="VAF")
dev.off()



pdf("hist_QUAL_4part.pdf",height=9,width=16)
par(mfrow=c(2,2))
hist(log10(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),5])),main="hist of QUAL for False Positive 94818-67127",xlab="log10(QUAL)")
hist(log10(as.numeric(tmp[match(intersect(tmp1,tmp2),tmp2),5])),main="hist of QUAL for True Positive 67127",xlab="log10(QUAL)")
hist(log10(as.numeric(tmp[match(intersect(tmp1x,tmp2x),tmp2x),5])),main="hist of QUAL for consistent 60880",xlab="log10(QUAL)")
hist(log10(as.numeric(tmp[setdiff(match(intersect(tmp1,tmp2),tmp2),match(intersect(tmp1x,tmp2x),tmp2x)),5])),main="hist of QUAL for inconsistent 67127-60880",xlab="log10(QUAL)")
dev.off()


pdf("gt_vaf.pdf",height=9,width=16)
par(mfrow=c(1,3))
boxplot(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),8]))~as.numeric(gt[[10]][match(intersect(tmp2,tmp1),tmp1)]),main="True Positive 67127")
boxplot(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])/(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])+as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),8])) ~ as.numeric(gt[[10]][match(intersect(tmp2x,tmp1x),tmp1x)]),main="consistent 60880")
boxplot(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])/(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])+as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),8])) ~ as.numeric(gt[[10]][setdiff(match(intersect(tmp2,tmp1),tmp1),match(intersect(tmp2x,tmp1x),tmp1x))]),main="inconsistent 67127 - 60880")
dev.off()

setwd("/sc/arion/scratch/wangm08/liver_gtex/outdir/liftover/")
a<-read.table("SRR1069141.liftover.g.vcf",sep="\t",as.is=T)
tmp=cbind(a[[1]],a[[2]],a[[4]],a[[5]],a[[6]])
tmp=cbind(tmp,sapply(a[[8]],function(z) {ifelse(length(which(grepl("DP=",strsplit(z,";")[[1]])))>0, strsplit(strsplit(z,";")[[1]][which(grepl("DP=",strsplit(z,";")[[1]]))],"DP=")[[1]][2], NA)}))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) length(which(as.numeric(strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("GT",strsplit(a[[9]][z],":")[[1]]))],"[||/]")[[1]])>0))))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("AD",strsplit(a[[9]][z],":")[[1]]))],",")[[1]][1]))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) sum(as.integer(strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("AD",strsplit(a[[9]][z],":")[[1]]))],",")[[1]][-1]))))
colnames(tmp)=c("Chr","Pos","Ref","Alt","QUAL","DP","GT","REFR","ALTR")
rna=tmp
tmp2=paste(rna[,1],rna[,2],sep="_")
tmp2x=paste(rna[,1],rna[,2],rna[,3],rna[,7],sep="_")
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>10)]))/length(which(as.numeric(tmp[,6])>10)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>10)])),length(which(as.numeric(tmp[,6])>10)))
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>9)]))/length(which(as.numeric(tmp[,6])>9)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>9)])),length(which(as.numeric(tmp[,6])>9)))
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>19)]))/length(which(as.numeric(tmp[,6])>19)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>19)])),length(which(as.numeric(tmp[,6])>19)))
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>29)]))/length(which(as.numeric(tmp[,6])>29)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>29)])),length(which(as.numeric(tmp[,6])>29)))
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>39)]))/length(which(as.numeric(tmp[,6])>39)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>39)])),length(which(as.numeric(tmp[,6])>39)))


rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>69)]))/length(which(as.numeric(tmp[,5])>69)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>69)])),length(which(as.numeric(tmp[,5])>69)))
rec
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>99)]))/length(which(as.numeric(tmp[,5])>99)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>99)])),length(which(as.numeric(tmp[,5])>99)))
rec
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>149)]))/length(which(as.numeric(tmp[,5])>149)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>149)])),length(which(as.numeric(tmp[,5])>149)))
rec
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>199)]))/length(which(as.numeric(tmp[,5])>199)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>199)])),length(which(as.numeric(tmp[,5])>199)))


rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.1/0.9)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.1/0.9)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.1/0.9)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.1/0.9)))
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.2/0.8)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.2/0.8)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.2/0.8)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.2/0.8)))
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.3/0.7)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.3/0.7)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.3/0.7)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.3/0.7)))
rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.4/0.6)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.4/0.6)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.4/0.6)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.4/0.6)))
setwd("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/")
pdf("hist_VAF.pdf",height=9,width=16)
par(mfrow=c(1,2))
hist(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),8])),main="hist of VAF for False Positive",xlab="VAF")
hist(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),8])),main="hist of VAF for True Positive",xlab="VAF")
dev.off()
pdf("hist_QUAL.pdf",height=9,width=16)
par(mfrow=c(1,2))
hist(log10(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),5])),main="hist of QUAL for False Positive",xlab="log10(QUAL)")
hist(log10(as.numeric(tmp[match(intersect(tmp1,tmp2),tmp2),5])),main="hist of QUAL for True Positive",xlab="log10(QUAL)")
dev.off()





pdf("hist_dp_4part.pdf",height=9,width=16)
par(mfrow=c(2,2))
hist(log10(as.numeric(tmp[match(intersect(tmp1,tmp2),tmp2),6])),main="hist of DP for True Positive 68704",xlab="log10(DP)")
hist(log10(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),6])),main="hist of DP for False Positive 177574-68704",xlab="log10(DP)")
hist(log10(as.numeric(tmp[match(intersect(tmp1x,tmp2x),tmp2x),6])),main="hist of DP for consistent 61739",xlab="log10(DP)")
hist(log10(as.numeric(tmp[setdiff(match(intersect(tmp1,tmp2),tmp2),match(intersect(tmp1x,tmp2x),tmp2x)),6])),main="hist of DP for inconsistent 68704 - 61739",xlab="log10(DP)")
dev.off()



pdf("hist_VAF_4part.pdf",height=9,width=16)
par(mfrow=c(2,2))
hist(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),8])),main="hist of VAF for False Positive 177574-68704",xlab="VAF")
hist(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),8])),main="hist of VAF for True Positive 68704",xlab="VAF")
hist(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])/(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])+as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),8])),main="hist of VAF for consistent 61739",xlab="VAF")
hist(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])/(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])+as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),8])),main="hist of VAF for inconsistent 68704-61739",xlab="VAF")
dev.off()



pdf("hist_QUAL_4part.pdf",height=9,width=16)
par(mfrow=c(2,2))
hist(log10(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),5])),main="hist of QUAL for False Positive 177574-68704",xlab="log10(QUAL)")
hist(log10(as.numeric(tmp[match(intersect(tmp1,tmp2),tmp2),5])),main="hist of QUAL for True Positive 68704",xlab="log10(QUAL)")
hist(log10(as.numeric(tmp[match(intersect(tmp1x,tmp2x),tmp2x),5])),main="hist of QUAL for consistent 61739",xlab="log10(QUAL)")
hist(log10(as.numeric(tmp[setdiff(match(intersect(tmp1,tmp2),tmp2),match(intersect(tmp1x,tmp2x),tmp2x)),5])),main="hist of QUAL for inconsistent 68704 - 61739",xlab="log10(QUAL)")
dev.off()

pdf("gt_vaf.pdf",height=9,width=16)
par(mfrow=c(1,3))
boxplot(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),8]))~as.numeric(gt[[10]][match(intersect(tmp2,tmp1),tmp1)]),main="True Positive 68704")
boxplot(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])/(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])+as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),8])) ~ as.numeric(gt[[10]][match(intersect(tmp2x,tmp1x),tmp1x)]),main="consistent 61739")
boxplot(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])/(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])+as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),8])) ~ as.numeric(gt[[10]][setdiff(match(intersect(tmp2,tmp1),tmp1),match(intersect(tmp2x,tmp1x),tmp1x))]),main="inconsistent 68704 - 61739")
dev.off()

=======================================================================================================================
setwd("/sc/arion/scratch/wangm08/liver_gtex/outdir/liftover/")
a<-read.table("SRR1069141_3.liftover.full.vcf",sep="\t",as.is=T)
tmp=cbind(a[[1]],a[[2]],a[[4]],a[[5]],a[[6]])
tmp=cbind(tmp,sapply(a[[8]],function(z) {ifelse(length(which(grepl("DP=",strsplit(z,";")[[1]])))>0, strsplit(strsplit(z,";")[[1]][which(grepl("DP=",strsplit(z,";")[[1]]))],"DP=")[[1]][2], NA)}))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) length(which(as.numeric(strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("GT",strsplit(a[[9]][z],":")[[1]]))],"[||/]")[[1]])>0))))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("AD",strsplit(a[[9]][z],":")[[1]]))],",")[[1]][1]))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) sum(as.integer(strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("AD",strsplit(a[[9]][z],":")[[1]]))],",")[[1]][-1]))))
colnames(tmp)=c("Chr","Pos","Ref","Alt","QUAL","DP","GT","REFR","ALTR")
row.names(tmp)=NULL
rna=tmp
gt<-read.table("/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir/sample/GTEX-U8XE.vcf",as.is=T)
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/geno_num.R")
gt[[10]]<-as.character(geno_num(gt[[10]]))
gt[[1]]=paste("chr",gt[[1]],sep="")
tmp1=paste(gt[[1]],gt[[2]],sep="_")
tmp2=paste(rna[,1],rna[,2],sep="_")
tmp1x=paste(gt[[1]],gt[[2]],gt[[4]],gt[[10]],sep="_")
tmp2x=paste(rna[,1],rna[,2],rna[,3],rna[,7],sep="_")

compare_2_vcf<-fucntion(tmp1,tmp2){

rec=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>=10)]))/length(which(as.numeric(tmp[,6])>=10)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>=10)])),length(which(as.numeric(tmp[,6])>=10)))
#res
#rec
rec=cbind(rec,c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>=20)]))/length(which(as.numeric(tmp[,6])>=20)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>=20)])),length(which(as.numeric(tmp[,6])>=20))))
#rec
rec=cbind(rec,c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>=30)]))/length(which(as.numeric(tmp[,6])>=30)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>=30)])),length(which(as.numeric(tmp[,6])>=30))))
#rec
rec=cbind(rec,c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>=40)]))/length(which(as.numeric(tmp[,6])>=40)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,6])>=40)])),length(which(as.numeric(tmp[,6])>=40))))

colnames(rec)=c(10,20,30,40)
row.names(rec)=c("Accuracy","TP","Prediction")


rec2=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>=70)]))/length(which(as.numeric(tmp[,5])>=70)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>=70)])),length(which(as.numeric(tmp[,5])>=70)))
rec2=cbind(rec2,c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>=100)]))/length(which(as.numeric(tmp[,5])>=100)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>=100)])),length(which(as.numeric(tmp[,5])>=100))))
rec2=cbind(rec2,c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>=150)]))/length(which(as.numeric(tmp[,5])>=150)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>=150)])),length(which(as.numeric(tmp[,5])>=150))))
rec2=cbind(rec2,c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>=200)]))/length(which(as.numeric(tmp[,5])>=200)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,5])>=200)])),length(which(as.numeric(tmp[,5])>=200))))


colnames(rec2)=c(70,100,150,200)
row.names(rec2)=c("Accuracy","TP","Prediction")

rec3=c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.1/0.9)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.1/0.9)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.1/0.9)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.1/0.9)))
rec3=cbind(rec3,c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.2/0.8)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.2/0.8)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.2/0.8)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.2/0.8))))
rec3=cbind(rec3,c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.3/0.7)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.3/0.7)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.3/0.7)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.3/0.7))))
rec3=cbind(rec3,c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.4/0.6)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.4/0.6)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.4/0.6)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.4/0.6))))
rec3=cbind(rec3,c(length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.5/0.5)]))/length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.5/0.5)),length(intersect(tmp1,tmp2[which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.5/0.5)])),length(which(as.numeric(tmp[,9])/as.numeric(tmp[,8])>=0.5/0.5))))

colnames(rec3)=c(0.1,0.2,0.3,0.4,0.5)
row.names(rec3)=c("Accuracy","TP","Prediction")

return rbind(rec1,rec2,rec3)
}


source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_2vcf.R")
write.csv(compare_2vcf(tmp1,tmp2),"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk35_result.csv",quote=F)  > write.csv(compare_2vcf(tmp1x,tmp2x),"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk35_result_geno.csv",quote=F)

setwd("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk3_res/")
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_2vcf_plot.R")
compare_2vcf_plot(tmp1,tmp2,tmp1x,tmp2x)

==========================================================================================================================================================



setwd("/sc/arion/scratch/wangm08/liver_gtex/outdir/liftover/")
a<-read.table("SRR1069141_2.liftover.full.vcf",sep="\t",as.is=T)
tmp=cbind(a[[1]],a[[2]],a[[4]],a[[5]],a[[6]])
tmp=cbind(tmp,sapply(a[[8]],function(z) {ifelse(length(which(grepl("DP=",strsplit(z,";")[[1]])))>0, strsplit(strsplit(z,";")[[1]][which(grepl("DP=",strsplit(z,";")[[1]]))],"DP=")[[1]][2], NA)}))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) length(which(as.numeric(strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("GT",strsplit(a[[9]][z],":")[[1]]))],"[||/]")[[1]])>0))))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("AD",strsplit(a[[9]][z],":")[[1]]))],",")[[1]][1]))
tmp=cbind(tmp,sapply(c(1:nrow(a)),function(z) sum(as.integer(strsplit(strsplit(a[[10]][z],":")[[1]][which(grepl("AD",strsplit(a[[9]][z],":")[[1]]))],",")[[1]][-1]))))
colnames(tmp)=c("Chr","Pos","Ref","Alt","QUAL","DP","GT","REFR","ALTR")
row.names(tmp)=NULL
rna=tmp
gt<-read.table("/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir/sample/GTEX-U8XE.vcf",as.is=T)
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/geno_num.R")
gt[[10]]<-as.character(geno_num(gt[[10]]))
gt[[1]]=paste("chr",gt[[1]],sep="")
tmp1=paste(gt[[1]],gt[[2]],sep="_")
tmp2=paste(rna[,1],rna[,2],sep="_")
tmp1x=paste(gt[[1]],gt[[2]],gt[[4]],gt[[10]],sep="_")
tmp2x=paste(rna[,1],rna[,2],rna[,3],rna[,7],sep="_")

source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_2vcf.R")
write.csv(compare_2vcf(tmp1,tmp2),"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk4_result.csv",quote=F)
write.csv(compare_2vcf(tmp1x,tmp2x),"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk4_result_geno.csv",quote=F)

setwd("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk4_res/")
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_2vcf_plot.R")
compare_2vcf_plot(tmp1,tmp2,tmp1x,tmp2x)



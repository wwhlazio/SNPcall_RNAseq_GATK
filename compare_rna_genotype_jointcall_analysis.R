#building compare_rna_genotype_jointcall
length(intersect(row.names(gt1),row.names(gtg)))
length(ind1)
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
length(ind1)
length(ind2)
count11x=count11/116
count22x=count22/116
length(which(rowSums(count22x[,2:3])>0.05))
length(which(rowSums(count11x[,2:3])>0.05))
length(which(rowSums(count1[,2:3])>0.05))
length(which(rowSums(count1[,2:3])>0.05 & count1[,4]<=0.75))
length(which(rowSums(count2[,2:3])>0.05 & count2[,4]<=0.75))
length(which(rowSums(count1[,2:3])>0.05 & count1[,4]<=0.5))
tmp1=dpx
tmp1[which(dpx<10)]=NA
tmp3=tmp1
tmp3[which(!is.na(tmp3))]=1
gt1=gtx*tmp3
ind1=which(apply(gt1[,],1,function(z) length(which(!is.na(z))))>0)
length(ind1)
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.25))
length(which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.25))
length(which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.75))
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.75))
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.5))
length(which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.5))
length(which(count11x[,4]<=0.75))
length(which(count22x[,4]<=0.75))
length(which(count22x[,4]<=0.5))
length(which(count11x[,4]<=0.5))
length(ind1)
ind1=match(intersect(row.names(gtg),row.names(gtx)),row.names(gtx))
ind2=match(setdiff(row.names(gtx),row.names(gtg)),row.names(gtx))
length(idn1)
length(idn1)
length(ind1)
length(which(rowSums(count1[,2:3])>0.05))
length(which(rowSums(count2[,2:3])>0.05))
length(which(rowSums(count22x[,2:3])>0.05))
length(which(rowSums(count11x[,2:3])>0.05))
length(which(rowSums(count22x[,2:3])>0.05))
length(which(rowSums(count11x[,2:3])>0.05))
length(which(rowSums(count1[,2:3])>0.05 & count1[,4]<=0.75))
length(which(rowSums(count2[,2:3])>0.05 & count2[,4]<=0.75))
length(which(rowSums(count1[,2:3])>0.05 & count1[,4]<=0.5))
length(which(rowSums(count2[,2:3])>0.05 & count2[,4]<=0.5))
length(intersect(row.names(gt1),row.names(gtg)))
tmp1=dpx
tmp1[which(dpx<10)]=NA
tmp3=tmp1
tmp3[which(!is.na(tmp3))]=1
gt1=gtx*tmp3
ind1=which(apply(gt1[,],1,function(z) length(which(!is.na(z))))>0)
gt1=gt1[ind1,]
length(intersect(row.names(gt1),row.names(gtg)))
length(idn1)
length(ind1)
tmp1=dpx
tmp1[which(dpx<10)]=NA
tmp3=tmp1
tmp3[which(!is.na(tmp3))]=1
gt1=gtx*tmp3
ind1=which(apply(gt1[,],1,function(z) length(which(!is.na(z))))>0)
gt1=gt1[ind1,]
dim(gt1)
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.25))
length(which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.25))
length(which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.75))
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.75))
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.5))
length(which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.5))
res
length(which(count11x[,4]<=0.75))
length(which(count22x[,4]<=0.75))
length(which(count22x[,4]<=0.5))
length(which(count11x[,4]<=0.5))
length(which(rowSums(count1[,2:3])>0.05 & count1[,4]<=0.5))
dim(count22x)
length(which(count22x[,4]<=0.5))
length(which(count11x[,4]<=0.75))
length(which(count11x[,4]<=0.5))
length(which(count22x[,4]<=0.5))
length(which(count11x[,4]<=0.5))
dim(count11x)
dim(count22x)
count11x[1:5,]
dim(gt1)
dim(gtx)
length(ind1)
dim(count11)
ind1=match(intersect(row.names(gtg),row.names(gt1)),row.names(gt1))
ind2=match(setdiff(row.names(gt1),row.names(gtg)),row.names(gt1))
length(ind11)
length(ind1)
length(ind2)
length(idn1)
length(ind1)
dim(gt1)
dim(count11x)
dim(count22x)
length(idn1)
length(idn1)
length(ind1)
dim(count11)
dim(count11x)
dim(count11x)
ind11=ind1[which(count11x[,4]<=0.5)]
length(ind11)
length(ind2)
length(ind1)
5949+893
ind2=match(setdiff(row.names(gt1),row.names(gtg)),row.names(gt1))
length(ind2)
ind1=match(intersect(row.names(gtg),row.names(gt1)),row.names(gt1))
length(ind1)
dim(g1)
dim(gt1)
893+5949
length(intersect(row.names(gtg),row.names(gt1)))
length(setdiff(row.names(gt1),row.names(gtg)))
dim(gt1)
gt1[1:5,1:5]
ind11=ind1[which(count11x[,4]<=0.5)]
length(ind11)
ind21=ind2[which(count22x[,4]<=0.5)]
length(ind21)
gt11=gt1[c(ind11,ind21),]
dim(gt11)
gt11[1:5,1:5]
gtg[1:5,1:5]
length(intersect(row.names(gt1),row.names(gt11)))
length(intersect(row.names(gt11),row.names(gtg)))
length(which(gt11==1))
length(which(gt11==2))
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
dim(gtg1)
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
length(which(gt111==gtg1))
dim(gt111)
dim(gtg1)
length(which(gt111==1)==which(gtg1==1))
length(intersect(which(gt111==1),which(gtg1==1)))
length(intersect(which(gt111==2),which(gtg1==2)))
dim(gt11)
19597/20778
9093/9650
length(whcih(gtg1==1))
length(which(gtg1==1))
length(which(gtg1==2))
19597/23704
11235/23704
9093/11235
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_rna_genotype_jointcall.R")
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_rna_genotype_jointcall.R")
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_rna_genotype_jointcall.R")
res<-compare_rna_genotype_jointcall(chr21_file,chr21geno_file)
length(which(count11x[,4]<=0.75))
ind11=ind1[which(count11x[,4]<=0.75)]
ind21=ind2[which(count22x[,4]<=0.75)]
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
length(which(gt1==1))
length(which(gt1==2))
length(which(gtg1==1))
length(which(gtg1==2))
length(intersect(which(gt111==1),which(gtg1==1)))
length(intersect(which(gt111==2),which(gtg1==2)))
11069/14163
23429/29786
11069/17374
23429/34477
ind11=ind1[which(count11x[,4]<=0.75)]
ind11=ind1[which(count11x[,4]<=0.5)]
ind21=ind2[which(count22x[,4]<=0.5)]
length(ind11)
length(ind21)
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
length(which(gt11==1))
length(which(gt11==2))
length(which(gt111==1))
length(which(gt111==2))
length(which(gtg1==2))
length(which(gtg1==1))
length(intersect(which(gt111==1),which(gtg1==1)))
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.5))
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.1))
length(which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.15))
length(which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.1))
length(which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.15))
336/(336+14)
379/(379+21)
ind11=ind1[which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.1)]
ind21=ind2[which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.1)]
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
length(gt11==1)
length(gt11==2)
length(which(gt11==2))
length(which(gt11==1))
length(which(gtg1==1))
length(which(gtg1==2))
length(intersect(which(gt111==1),which(gtg1==1)))
10588/10695
length(intersect(which(gt111==2),which(gtg1==2)))
5104/5118
ind11 = which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.15)
ind21 = which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.15)
loci_check<-function(ind11,ind21)
resx=matrix(2,3)
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
resx[1,1]=length(which(gtg1==1))
resx[2,1]=length(which(gtg1==2))
resx[1,2]=length(which(gt11==1))
resx[2,2]=length(which(gt11==2))
resx[1,3]=length(intersect(which(gt111==1),which(gtg1==1)))
resx[2,3]=length(intersect(which(gt111==2),which(gtg1==2)))
return(rex)
}
loci_check<-function(ind11,ind21){
resx=matrix(2,3)
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
resx[1,1]=length(which(gtg1==1))
resx[2,1]=length(which(gtg1==2))
resx[1,2]=length(which(gt11==1))
resx[2,2]=length(which(gt11==2))
resx[1,3]=length(intersect(which(gt111==1),which(gtg1==1)))
resx[2,3]=length(intersect(which(gt111==2),which(gtg1==2)))
return(rex)
}
loci_check<-function(ind11,ind21){
resx=matrix(2,3)
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
resx[1,1]=length(which(gtg1==1))
resx[2,1]=length(which(gtg1==2))
resx[1,2]=length(which(gt11==1))
resx[2,2]=length(which(gt11==2))
resx[1,3]=length(intersect(which(gt111==1),which(gtg1==1)))
resx[2,3]=length(intersect(which(gt111==2),which(gtg1==2)))
return(resx)
}
loci_check(ind11,ind21)
ind11 = which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.15)
ind21 = which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.15)
loci_check<-function(ind11,ind21){
resx=matrix(2,3)
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
resx[1,1]=length(which(gtg1==1))
resx[2,1]=length(which(gtg1==2))
resx[1,2]=length(which(gt11==1))
resx[2,2]=length(which(gt11==2))
resx[1,3]=length(intersect(which(gt111==1),which(gtg1==1)))
resx[2,3]=length(intersect(which(gt111==2),which(gtg1==2)))
return(resx)
}
loci_check(ind11,ind21)
loci_check<-function(ind11,ind21,gtg,gt1){
resx=matrix(2,3)
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
resx[1,1]=length(which(gtg1==1))
resx[2,1]=length(which(gtg1==2))
resx[1,2]=length(which(gt11==1))
resx[2,2]=length(which(gt11==2))
resx[1,3]=length(intersect(which(gt111==1),which(gtg1==1)))
resx[2,3]=length(intersect(which(gt111==2),which(gtg1==2)))
return(resx)
}
loci_check(ind11,ind21)
loci_check(ind11,ind21,gtg,gt1)
ind11 = ind1[which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.15)]
ind21 = ind2[which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.15)]
loci_check<-function(ind11,ind21,gtg,gt1){
resx=matrix(2,3)
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
resx[1,1]=length(which(gtg1==1))
resx[2,1]=length(which(gtg1==2))
resx[1,2]=length(which(gt11==1))
resx[2,2]=length(which(gt11==2))
resx[1,3]=length(intersect(which(gt111==1),which(gtg1==1)))
resx[2,3]=length(intersect(which(gt111==2),which(gtg1==2)))
return(resx)
}
loci_check(ind11,ind21,gtg,gt1)
resx=matrix(2,3)
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
resx[1,1]=length(which(gtg1==1))
resx[2,1]=length(which(gtg1==2))
resx[1,2]=length(which(gt11==1))
resx[2,2]=length(which(gt11==2))
resx[1,3]=length(intersect(which(gt111==1),which(gtg1==1)))
resx[2,3]=length(intersect(which(gt111==2),which(gtg1==2)))
loci_check<-function(ind11,ind21,gtg,gt1){
resx=matrix(0,nrow=2,ncol=3)
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
resx[1,1]=length(which(gtg1==1))
resx[2,1]=length(which(gtg1==2))
resx[1,2]=length(which(gt11==1))
resx[2,2]=length(which(gt11==2))
resx[1,3]=length(intersect(which(gt111==1),which(gtg1==1)))
resx[2,3]=length(intersect(which(gt111==2),which(gtg1==2)))
return(resx)
}
loci_check(ind11,ind21,gtg,gt1)
ind11 = which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.1)
ind11 = ind1[which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.11)]
ind11 = ind1[which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.1)]
ind21 = ind2[which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.1)]
loci_check(ind11,ind21,gtg,gt1)
length(ind1)
length(ind2)
ind1=match(intersect(row.names(gtg),row.names(gt1)),row.names(gt1))
ind2=match(setdiff(row.names(gt1),row.names(gtg)),row.names(gt1))
length(ind1)
length(ind2)
ind11=ind1[which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.1)]
dim(count11x)
ind21=ind2[which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.1)]
length(ind11)
length(ind21)
ind11 = ind1[which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.15)]
ind21 = ind2[which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.15)]
length(ind11)
length(ind21)
ind11 = ind1[which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.1)]
ind21=ind2[which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.1)]
length(ind11)
length(ind21)
ind21 = ind2[which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.15)]
ind11 = ind1[which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.15)]
length(ind11)
length(ind21)
loci_check(ind11,ind21,gtg,gt1)
length(ind11)
length(ind21)
ind11 = ind1[which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.1)]
ind21 = ind2[which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.1)]
loci_check(ind11,ind21,gtg,gt1)
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_rna_genotype_jointcall.R")
res<-compare_rna_genotype_jointcall(chr21_file,chr21geno_file)
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_rna_genotype_jointcall.R")
res<-compare_rna_genotype_jointcall(chr21_file,chr21geno_file)
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_rna_genotype_jointcall.R")
compare_rna_genotype_jointcall
res<-compare_rna_genotype_jointcall(chr21_file,chr21geno_file)
res
loci_check(ind11,ind21,gtg,gt1)
5104/5236
5104/5483
ind11 = ind1[which(rowSums(count11x[,2:3])>0.05 & count11x[,4]<=0.1)]
ind21 = ind2[which(rowSums(count22x[,2:3])>0.05 & count22x[,4]<=0.1)]
gt1=gtx*tmp3
ind1=which(apply(gt1[,],1,function(z) length(which(!is.na(z))))>0)
gt1=gt1[ind1,]
dim(gt1)
ind1=match(intersect(row.names(gtg),row.names(gt1)),row.names(gt1))
ind2=match(setdiff(row.names(gt1),row.names(gtg)),row.names(gt1))
gt11=gt1[c(ind11,ind21),]
gtg1=gtg[match(intersect(row.names(gt11),row.names(gtg)),row.names(gtg)),]
gt111=gt11[match(intersect(row.names(gt11),row.names(gtg)),row.names(gt11)),]
length(which(gt11==1))
length(which(gt11==2))
dim(gt11)
length(ind2)
dim(count22)
res
5788/6280
5788/6064
5104/5236
5104/5483
10588/10944
10588/10849
=============================================================================================================================================================

setwd("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/valid_rna_geno")
list.files()
#i=1
#load(paste("chr",i,".rda",sep=""))
#ls()
#res_chr1
for(i in 1:22){
load(paste("chr",i,".rda",sep=""))
}
#ls()
#res_chr21
resloci=matrix(0,nrow=22,ncol=6)
i=1
#get(paste("res_chr",i,sep=""))
#get(paste("res_chr",i,sep=""))[,11:12]
tmp=get(paste("res_chr",i,sep=""))
resloci[i,1:2]=tmp[,11]
resloci[i,3]=tmp[2,11]/tmp[1,11]
#resloci[i,3:5]=tmp[,12]
resloci[i,4:5]=tmp[,12]
resloci[i,6]=tmp[2,12]/tmp[1,12]
#resloci[1,]
for (i in 2:22){
tmp=get(paste("res_chr",i,sep=""))
resloci[i,1:2]=tmp[,11]
resloci[i,3]=tmp[2,11]/tmp[1,11]
resloci[i,4:5]=tmp[,12]
resloci[i,6]=tmp[2,12]/tmp[1,12]
}
#resloci
#res_chr1
#86905/127612
#22745/26770
#38605/46738
#4432/4913
#5796/6537
#14103/15249
#chr1 <- read.table("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/liftover_chr/liver_all_combine.chr1.liftover.vcf",as.is=T)
#dim(chr1)
#res_chr21
#dim(chr1)
#res_chr1
#chr1geno<-read.table("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/genotype_chr/1.geno.vcf",as.is=T)
#cmd = paste("grep \"#CHROM\"  ","/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/liftover_chr/liver_all_combine.chr1.liftover.vcf",sep="")
#header <- system(cmd,intern=T)
#samp<-strsplit(header,"\t")[[1]][10:147]
#dp=matrix(NA,nrow=nrow(chr1),ncol=length(samp))
#colnames(dp)=samp
#tmp=paste(chr1[,1],chr1[,2],chr1[,4],chr1[,5],sep="_")
#dim(gt)
#row.names(dp)<-tmp
#dim(gt1)
#dim(gtx)
#resloci
#res_chr1
#76493/86704
#48028/51326
#67855/76514
#43139/45759
resloci15=matrix(0,nrow=22,ncol=10)
i=1
#tmp=paste(chr1[,1],chr1[,2],chr1[,4],chr1[,5],sep="_")
tmp=get(paste("res_chr",i,sep=""))
#tmp
resloci15[i,1:3]=tmp[1,14:16]
resloci15[i,4]=tmp[1,16]/tmp[1,15]
resloci15[i,5]=tmp[1,16]/tmp[1,14]
#resloci16[i,6:8]=tmp[2,14:16]
resloci15[i,6:8]=tmp[2,14:16]
resloci15[i,9]=tmp[2,16]/tmp[2,15]
resloci15[i,10]=tmp[2,16]/tmp[2,14]
#resloci15[1,]
#48028/50024
for (i in 2:22){
tmp=get(paste("res_chr",i,sep=""))
resloci15[i,1:3]=tmp[1,14:16]
resloci15[i,4]=tmp[1,16]/tmp[1,15]
resloci15[i,5]=tmp[1,16]/tmp[1,14]
resloci15[i,6:8]=tmp[2,14:16]
resloci15[i,9]=tmp[2,16]/tmp[2,15]
resloci15[i,10]=tmp[2,16]/tmp[2,14]
}
#resloci15
#res_chr1
resloci10=matrix(0,nrow=22,ncol=10)
for (i in 1:22){
tmp=get(paste("res_chr",i,sep=""))
resloci10[i,1:3]=tmp[1,17:19]
resloci10[i,4]=tmp[1,19]/tmp[1,18]
resloci10[i,5]=tmp[1,19]/tmp[1,17]
resloci10[i,6:8]=tmp[2,17:19]
resloci10[i,9]=tmp[2,19]/tmp[2,18]
resloci10[i,10]=tmp[2,19]/tmp[2,17]
}
#resloci10
#5104/5483
#resloci
row.names(resloci)=c(1:22)
colnames(resloci)=c("#COV10VAF05CR85","#TP","ACCURACY","#COV10VAF5CR90","#TP","ACCURACY")
#resloci
#resloci15
colnames(resloci15)=c("#GT1","#GT1_Called","#GT1_TP","Accuracy_GT1","Recall_GT1","#GT2","#GT2_called","#GT2_TP","Accuracy_GT2","Recall_GT2")
colnames(resloci10)=c("#GT1","#GT1_Called","#GT1_TP","Accuracy_GT1","Recall_GT1","#GT2","#GT2_called","#GT2_TP","Accuracy_GT2","Recall_GT2")
#write.table(resloci,"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/valid_rna_geno/COV10_VAF05_CR85_CR90_loci_accuracy",sep="\t",quote=F,row.names=F,col.names=F)
#write.table(resloci15,"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/valid_rna_geno/COV10_VAF05_CR85_loci_geno_accuracy_recall",sep="\t",quote=F,row.names=F,col.names=F)
#write.table(resloci10,"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/valid_rna_geno/COV10_VAF05_CR90_loci_geno_accuracy_recall",sep="\t",quote=F,row.names=F,col.names=F)
savehistory("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_rna_genotype_jointcall_analysis.Rhistory")
write.table(resloci10,"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/valid_rna_geno/COV10_VAF05_CR90_loci_geno_accuracy_recall",sep="\t",quote=F,row.names=T,col.names=T)
write.table(resloci15,"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/valid_rna_geno/COV10_VAF05_CR85_loci_geno_accuracy_recall",sep="\t",quote=F,row.names=T,col.names=T)
write.table(resloci,"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/valid_rna_geno/COV10_VAF05_CR85_CR90_loci_accuracy",sep="\t",quote=F,row.names=T,col.names=T)
savehistory("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_rna_genotype_jointcall_analysis.Rhistory")

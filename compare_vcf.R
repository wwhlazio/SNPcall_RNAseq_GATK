load("/sc/arion/scratch/wangm08/liver_gtex/genotype/compare_rna_genotype.RData")
source("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/geno_num.R")
setwd("/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir/sample/")
maty=mat[match(ind[,2],mat[,2]),]


compare_vcf<-function(gt,rna){
        gt[[10]]<-as.character(geno_num(gt[[10]]))
	gt[[1]]=paste("chr",gt[[1]],sep="")
	res=matrix(0,nrow=1,ncol=4)
	res[1]=nrow(gt)
	res[2]=nrow(rna)
        tmp1=paste(gt[[1]],gt[[2]],sep="_")
	tmp2=paste(rna[,1],rna[,2],sep="_")
	res[3]=length(intersect(tmp1,tmp2))
	tmp1=paste(gt[[1]],gt[[2]],gt[[4]],gt[[10]],sep="_")
        tmp2=paste(rna[,1],rna[,2],rna[,3],rna[,7],sep="_")
	res[4]=length(intersect(tmp1,tmp2))
	return(res)
}

i=1
 file=paste(maty[i,2],".vcf",sep="")
        gt<-read.table(file,as.is=T)
        rna<-get(maty[i,1])
        res<-compare_vcf(gt,rna)



for(i in 2:nrow(maty)){
	file=paste(maty[i,2],".vcf",sep="")
	gt<-read.table(file,as.is=T)
	rna<-get(maty[i,1])
	res<-rbind(res,compare_vcf(gt,rna))
	if(i>1){
		print(cbind(maty[c(1:i),],res[c(1:i),]))
	}
	print(i)
}




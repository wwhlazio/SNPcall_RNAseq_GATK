setwd("/sc/arion/scratch/wangm08/liver_gtex/outdir/vpass")


geno_num<-function(seq){
	res=matrix(0,nrow=1,ncol=length(seq))
	for (i in 1:length(seq)){
		ge<-strsplit(seq[i],'[|]')[[1]]
		if (ge[1]=='.'){
			res[i]=NA
		}else{
			res[i]=as.integer(ge[1])+as.integer(ge[2])
		}
	}
	return(res)
}


extract<-function(startLine,endLine,ind){

res=''
cmd = paste("sed -n \'",endLine,"q;",startLine,",",endLine,"p\' /sc/arion/scratch/wangm08/liver_gtex/genotype/phg000830.v1.GTEx_WGS_additional.genotype-calls-vcf.c1/wgs635_genotype.vcf",sep="")
tmp=system(cmd,intern=TRUE)

for(i in 1:length(tmp) ){
	tmpy <- strsplit(tmp[i],"\t")[[1]]
	tmpx<-geno_num( tmpy[ind])
	res=rbind(res,matrix(c(tmpy[c(1,2,4,5)],tmpx),nrow=1))
}

return(res)
}


startLine=30+10000*(i-1)+1
endLine=30+10000*i+1
lsfile=paste(workdir,"extract_",i,".lsf",sep="")
stderr=paste(workdir,"extract_",i,".stderr",sep="")
stdout=paste(workdir,"extract_",i,".stdout",sep="")
output=paste(outdir,"extract_",i,sep="")

sink(lsfile)
cat("#!/bin/bash\n")
cat("#BSUB -P acc_BD2K\n")
cat("#BSUB -q premium\n")
cat("#BSUB -J extract_")
cat(i)
cat("\n")
cat("#BSUB -R \"rusage[mem=5000]\"\n")
cat("#BSUB -W 00:30\n")
cat("#BSUB -o ")
cat(stdout)
cat("\n")
cat("#BSUB -eo ")
cat(stderr)
cat("\n")
cat("#BSUB -L /bin/bash\n")
cat("module load R\n")
cat("Rscript -e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/geno_num.R\\\")\" \\\n")
cat("-e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/extract.R\\\")\" \\\n")
cat("-e \"ind=read.table(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/livergenosampid\\\",as.is=T,sep=\\\"\\t\\\")\" \\\n")
cat("-e \"ss<-extract(")
cat(startLine)
cat(",")
cat(endLine)
cat(",ind[[1]])\" \\\n")
cat("-e \"colnames(ss)=c(\\\"Chr\\\",\\\"Pos\\\",\\\"Ref\\\",\\\"Alt\\\",ind[[2]])\" \\\n")
cat("-e \"write.table(ss,\\\"")
cat(output)
cat("\\\",col.names=T,row.names=F,quote=F,sep=\\\"\\t\\\")\"\n")
sink()
cat(paste("bsub < ",lsfile,sep=""),file=sh_file,sep="\n")

for( i in 2:4066){
startLine=30+10000*(i-1)+1
endLine=30+10000*i+1
lsfile=paste(workdir,"extract_",i,".lsf",sep="")
stderr=paste(workdir,"extract_",i,".stderr",sep="")
stdout=paste(workdir,"extract_",i,".stdout",sep="")
output=paste(outdir,"extract_",i,sep="")

sink(lsfile)
cat("#!/bin/bash\n")
cat("#BSUB -P acc_BD2K\n")
cat("#BSUB -q premium\n")
cat("#BSUB -J extract_")
cat(i)
cat("\n")
cat("#BSUB -R \"rusage[mem=5000]\"\n")
cat("#BSUB -W 00:30\n")
cat("#BSUB -o ")
cat(stdout)
cat("\n")
cat("#BSUB -eo ")
cat(stderr)
cat("\n")
cat("#BSUB -L /bin/bash\n")
cat("module load R\n")
cat("Rscript -e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/geno_num.R\\\")\" \\\n")
cat("-e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/extract.R\\\")\" \\\n")
cat("-e \"ind=read.table(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/livergenosampid\\\",as.is=T,sep=\\\"\\t\\\")\" \\\n")
cat("-e \"ss<-extract(")
cat(startLine)
cat(",")
cat(endLine)
cat(",ind[[1]])\" \\\n")
cat("-e \"colnames(ss)=c(\\\"Chr\\\",\\\"Pos\\\",\\\"Ref\\\",\\\"Alt\\\",ind[[2]])\" \\\n")
cat("-e \"write.table(ss,\\\"")
cat(output)
cat("\\\",col.names=T,row.names=F,quote=F,sep=\\\"\\t\\\")\"\n")
sink()
cat(paste("bsub < ",lsfile,sep=""),file=sh_file,sep="\n",append=T)
}


i=4067
startLine=30+10000*(i-1)+1
endLine=40669411
lsfile=paste(workdir,"extract_",i,".lsf",sep="")
stderr=paste(workdir,"extract_",i,".stderr",sep="")
stdout=paste(workdir,"extract_",i,".stdout",sep="")
output=paste(outdir,"extract_",i,sep="")

sink(lsfile)
cat("#!/bin/bash\n")
cat("#BSUB -P acc_BD2K\n")
cat("#BSUB -q premium\n")
cat("#BSUB -J extract_")
cat(i)
cat("\n")
cat("#BSUB -R \"rusage[mem=5000]\"\n")
cat("#BSUB -W 00:30\n")
cat("#BSUB -o ")
cat(stdout)
cat("\n")
cat("#BSUB -eo ")
cat(stderr)
cat("\n")
cat("#BSUB -L /bin/bash\n")
cat("module load R\n")
cat("Rscript -e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/geno_num.R\\\")\" \\\n")
cat("-e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/extract.R\\\")\" \\\n")
cat("-e \"ind=read.table(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/livergenosampid\\\",as.is=T,sep=\\\"\\t\\\")\" \\\n")
cat("-e \"ss<-extract(")
cat(startLine)
cat(",")
cat(endLine)
cat(",ind[[1]])\" \\\n")
cat("-e \"colnames(ss)=c(\\\"Chr\\\",\\\"Pos\\\",\\\"Ref\\\",\\\"Alt\\\",ind[[2]])\" \\\n")
cat("-e \"write.table(ss,\\\"")
cat(output)
cat("\\\",col.names=T,row.names=F,quote=F,sep=\\\"\\t\\\")\"\n")
sink()
cat(paste("bsub < ",lsfile,sep=""),file=sh_file,sep="\n",append=T)



workdir = "/sc/arion/scratch/wangm08/liver_gtex/genotype/workdir/sample/"
outdir = "/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir/sample/"
sh_file = "/sc/arion/scratch/wangm08/liver_gtex/genotype/workdir/sample/extract_geno_sample.sh"

for( i in 1:nrow(ind)){
lsfile=paste(workdir,"extract_",ind[i,2],".lsf",sep="")
stderr=paste(workdir,"extract_",ind[i,2],".stderr",sep="")
stdout=paste(workdir,"extract_",ind[i,2],".stdout",sep="")
output=paste(outdir,ind[i,2],".vcf",sep="")

sink(lsfile)
cat("#!/bin/bash\n")
cat("#BSUB -P acc_BD2K\n")
cat("#BSUB -q premium\n")
cat("#BSUB -J extract_")
cat(i)
cat("\n")
cat("#BSUB -R \"rusage[mem=10000]\"\n")
cat("#BSUB -W 01:00\n")
cat("#BSUB -o ")
cat(stdout)
cat("\n")
cat("#BSUB -eo ")
cat(stderr)
cat("\n")
cat("#BSUB -L /bin/bash\n")
cat("module load bcftools\n")
cat("bcftools view -s ")
cat(ind[i,2])
cat(" -c 1 /sc/arion/scratch/wangm08/liver_gtex/genotype/phg000830.v1.GTEx_WGS_additional.genotype-calls-vcf.c1/GTEx_Analysis_20160115_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.PIR.vcf.gz  > ")
cat(output)
cat("\n")
sink()
cat(paste("bsub < ",lsfile,sep=""),file=sh_file,sep="\n",append=T)
}


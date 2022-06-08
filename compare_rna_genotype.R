compare_rna_genotype<-function(rna,genotype){
	
}


setwd("/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir")
files=list.files()


load_geno<-function(start,end,outfile){
	setwd("/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir")
	files=list.files()
	res=''
	for(i in c(start:end)){
		a=read.tables(files[i],header=T,as.is=T,sep="\t")
		tmp=apply(a[,5:120],1,function(z) length(which(z>0)))
		res=rbind(res,a[which(tmp>0)])
	}
	
	write.table(res[-1,],outfile,col.names=T,row.names=F,sep="\t",quote=F)
}
cmd = 

for (i in 2:8){
	start=500*(i-1)+1
	end=500*i
	outfile=paste("\"/sc/arion/scratch/wangm08/liver_gtex/genotype/gt_",i,"\"",sep="")
	cmd = "Rscript -e \"source(\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/load_geno.R\")\""
        cmd1 = paste(paste(" -e \"load_geno(",paste(start,end,outfile,sep=","),sep=""),")\"",sep="")
	cmd = paste(cmd,cmd1,sep="")
	print(cmd)
	
}

start=4001
end=4067
outfile="/sc/arion/scratch/wangm08/liver_gtex/genotype/gt_9"
cmd = "Rscript -e \"source(\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/load_geno.R\")\""
cmd1 = paste(paste(" -e \"load_geno(",paste(start,end,outfile,sep=","),sep=""),")\"",sep="")
cmd = paste(cmd,cmd1,sep="")

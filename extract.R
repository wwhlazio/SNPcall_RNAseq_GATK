extract<-function(startLine,endLine,ind){

res=''
cmd = paste("sed -n \'",endLine,"q;",startLine,",",endLine,"p\' /sc/arion/scratch/wangm08/liver_gtex/genotype/phg000830.v1.GTEx_WGS_additional.genotype-calls-vcf.c1/wgs635_genotype.vcf",sep="")
tmp=system(cmd,intern=TRUE)

for(i in 1:length(tmp) ){
        tmpy <- strsplit(tmp[i],"\t")[[1]]
        tmpx<-geno_num( tmpy[ind])
        res=rbind(res,matrix(c(tmpy[c(1,2,4,5)],tmpx),nrow=1))
}
res=res[-1,]
return(res)
}


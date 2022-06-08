load_geno<-function(start,end,outfile){
        setwd("/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir")
        files=list.files()
        res=''
        for(i in c(start:end)){
                a=read.table(files[i],header=T,as.is=T,sep="\t")
                tmp=apply(a[,5:120],1,function(z) length(which(z>0)))
                res=rbind(res,a[which(tmp>0),])
        }

        write.table(res[-1,],outfile,col.names=T,row.names=F,sep="\t",quote=F)
}


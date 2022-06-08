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


extrac_bychr<-fucntion(chr){
cmd = paste("bcftools view -S genotype/geno_samp -c 1 /sc/arion/scratch/wangm08/liver_gtex/genotype/phg000830.v1.GTEx_WGS_additional.genotype-calls-vcf.c1/GTEx_Analysis_20160115_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.PIR.vcf.gz -r ",chr," > /sc/arion/scratch/wangm08/liver_gtex/genotype/outdir/bychr/chr",chr,".geno.vcf",sep="")
system(cmd)
}

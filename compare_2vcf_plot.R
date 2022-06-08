compare_2vcf_plot<-function(tmp1,tmp2,tmp1x,tmp2x){

pdf("hist_dp_4part.pdf",height=9,width=16)
par(mfrow=c(2,2))
hist(log10(as.numeric(tmp[match(intersect(tmp1,tmp2),tmp2),6])),main="hist of DP for True Positive ",xlab="log10(DP)")
hist(log10(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),6])),main="hist of DP for False Positive",xlab="log10(DP)")
hist(log10(as.numeric(tmp[match(intersect(tmp1x,tmp2x),tmp2x),6])),main="hist of DP for consistent GT in positive",xlab="log10(DP)")
hist(log10(as.numeric(tmp[setdiff(match(intersect(tmp1,tmp2),tmp2),match(intersect(tmp1x,tmp2x),tmp2x)),6])),main="hist of DP for inconsistent GT in positive",xlab="log10(DP)")
dev.off()



pdf("hist_VAF_4part.pdf",height=9,width=16)
par(mfrow=c(2,2))
hist(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),8])),main="hist of VAF for False Positive",xlab="VAF")
hist(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),8])),main="hist of VAF for True Positive",xlab="VAF")
hist(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])/(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])+as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),8])),main="hist of VAF for consistent GT in positive",xlab="VAF")
hist(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])/(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])+as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),8])),main="hist of VAF for inconsistent GT in positive",xlab="VAF")
dev.off()



pdf("hist_QUAL_4part.pdf",height=9,width=16)
par(mfrow=c(2,2))
hist(log10(as.numeric(tmp[match(setdiff(tmp2,tmp1),tmp2),5])),main="hist of QUAL for False Positive",xlab="log10(QUAL)")
hist(log10(as.numeric(tmp[match(intersect(tmp1,tmp2),tmp2),5])),main="hist of QUAL for True Positive",xlab="log10(QUAL)")
hist(log10(as.numeric(tmp[match(intersect(tmp1x,tmp2x),tmp2x),5])),main="hist of QUAL for consistent GT in positive",xlab="log10(QUAL)")
hist(log10(as.numeric(tmp[setdiff(match(intersect(tmp1,tmp2),tmp2),match(intersect(tmp1x,tmp2x),tmp2x)),5])),main="hist of QUAL for inconsistent GT in positive",xlab="log10(QUAL)")
dev.off()

pdf("gt_vaf.pdf",height=9,width=16)
par(mfrow=c(1,3))
boxplot(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])/(as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),9])+as.numeric(tmp[match(intersect(tmp2,tmp1),tmp2),8]))~as.numeric(gt[[10]][match(intersect(tmp2,tmp1),tmp1)]),main="VAF_RNA of GT_WGS for True Positive")
boxplot(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])/(as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),9])+as.numeric(tmp[match(intersect(tmp2x,tmp1x),tmp2x),8])) ~ as.numeric(gt[[10]][match(intersect(tmp2x,tmp1x),tmp1x)]),main="TP and consistent GT")
boxplot(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])/(as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),9])+as.numeric(tmp[setdiff(match(intersect(tmp2,tmp1),tmp2),match(intersect(tmp2x,tmp1x),tmp2x)),8])) ~ as.numeric(gt[[10]][setdiff(match(intersect(tmp2,tmp1),tmp1),match(intersect(tmp2x,tmp1x),tmp1x))]),main="TP but inconsistent GT")
dev.off()
}

compare_2vcf<-function(tmp1,tmp2){

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

return(cbind(rec,rec2,rec3))
}


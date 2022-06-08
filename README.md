# SNPcall_RNAseq_GATK
implement_rnaseq_snpcall_gatk3.py is teh final pipeline. comparing with gatk4 and found that RNAseq with 2 step alignment performs much better.
run single sample snpcall first and genotype together performs better
R code is used to select the best cut off based on comparing with DNA genotype

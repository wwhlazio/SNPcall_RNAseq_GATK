#/usr/bin/python3

import subprocess
import os


def run_fastq_filter(data_dir,outdir,workdir):
    acc = "premium"
    time1 = "4:00"
    mem1 ="2500"
    nthr = 4
    
    cmd = "ls "+data_dir+"*.fastq"
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    samp = []
    for i in range(len(files)):
        tmpx = files[i].split("/")[-1]    
        tmp = tmpx.split("_")[0]
        if tmp not in samp:
            samp.append(tmp)

    shfile = workdir + "filter.sh"
    fx = open(shfile,"w")
    for i in range(len(samp)-1):
        lsfile = workdir + samp[i] + "_filter.lsf"
        inputfile1 = data_dir + samp[i] + "_1.fastq"
        inputfile2 = data_dir + samp[i] + "_2.fastq"
        outputfile1 = outdir + samp[i] + "_1.fq.gz"
        outputfile2 = outdir + samp[i] + "_2.fq.gz"

        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -n %s\n" % nthr)
        f.write("#BSUB -J filter_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %sfastq_filter_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %sfastq_filter_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")
        f.write("module load fastp\n")
        f.write("fastp -w %d -i %s -o %s -I %s -O %s\n" % (nthr,inputfile1,outputfile1,inputfile2,outputfile2))
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

    

def run_star(data_dir,outdir,workdir,genomedir):
     
    acc = "premium"
    time1 = "4:00"
    mem1 ="60000"
    nthr = 4  
    
    cmd = "ls "+data_dir+"*.fq.gz"
    print(cmd)
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    samp = []
    for i in range(len(files)):
        tmpx = files[i].split("/")[-1]
        tmp = tmpx.split("_")[0]
        if tmp not in samp:
            samp.append(tmp)   

    print(files)    

    shfile = workdir + "align.sh"
    fx = open(shfile,"w")
    for i in range(len(samp)-1):
        lsfile = workdir + samp[i] + "_align.lsf"
        inputfile1 = data_dir + samp[i] + "_1.fq.gz"
        inputfile2 = data_dir + samp[i] + "_2.fq.gz"
	outdir_tmp = outdir + samp[i] + "/"
        cmd = "mkdir " + outdir_tmp
        os.system(cmd)
        outputprefix = outdir_tmp + samp[i]
        
        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J align_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %salign_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %salign_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")
        f.write("module load star\n")
        f.write("STAR  --genomeDir %s  --outFilterMultimapNmax 1 --outReadsUnmapped Fastx  --outSAMtype BAM SortedByCoordinate  --twopassMode Basic  --runThreadN %d --readFilesCommand zcat --quantMode GeneCounts --readFilesIn %s %s --outFileNamePrefix %s\n" % (genomedir,nthr,inputfile1,inputfile2,outputprefix))
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

def run_picard(data_dir,workdir,outdir):
    
    acc = "premium"
    time1 = "4:00"
    mem1 ="20000"


    cmd = "find "+data_dir+" -name \"*.bam\""
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    samp = []
    for i in range(len(files)):
        tmpx = files[i].split("/")[-1]
        tmp = tmpx.split("Aligned.sortedByCoord.out.bam")[0]
        if tmp not in samp:
            samp.append(tmp)
    
    shfile = workdir + "picard.sh"
    
    fx = open(shfile,"w")
    for i in range(len(samp)-1):
        lsfile = workdir + samp[i] + "_picard.lsf"
              
 
        workdir_tmp = workdir + samp[i] + "/"
        cmd = "mkdir " + workdir_tmp
        os.system(cmd)
        
        outdir_tmp = outdir + samp[i] + "/"
        cmd = "mkdir " + outdir_tmp
        os.system(cmd)
        
        outfile1 = outdir_tmp + samp[i] + ".mkDup.bam"
        outfile2 = outdir_tmp + samp[i] + ".marked_dup_metrics.txt"
        outfile3 = outdir_tmp + samp[i] +  ".sorted.bam"
        
        
        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J picard_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %spicard_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %spicard_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")
        f.write("module load java\n")
        f.write("module load samtools\n")
        f.write("java -jar -Xmx10g \"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/picard.jar\" MarkDuplicates --TMP_DIR %s --INPUT %s --OUTPUT %s --METRICS_FILE %s\n" % (workdir_tmp,files[i],outfile1,outfile2))
        f.write("java -jar -Xmx10g \"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/picard.jar\" SortSam --TMP_DIR %s --INPUT %s --OUTPUT %s  --SORT_ORDER coordinate \n" % (workdir_tmp,outfile1,outfile3))
        f.write("samtools index %s\n" % outfile3)
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()


def run_gatk(data_dir,workdir,outdir,genome,knownsites):
    
    acc = "premium"
    time1 = "8:00"
    mem1 ="20000"



    cmd = "find "+data_dir+" -name \"*sorted.bam\""
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    samp = []
    for i in range(len(files)):
        tmpx = files[i].split("/")[-1]
        tmp = tmpx.split(".sorted.bam")[0]
        if tmp not in samp:
            samp.append(tmp)

    shfile = workdir + "gatk.sh"
    
    
    
    fx = open(shfile,"w")
    for i in range(len(samp)-1):
        lsfile = workdir + samp[i] + "_gatk.lsf"


        workdir_tmp = workdir + samp[i] + "/"
        cmd = "mkdir " + workdir_tmp
        os.system(cmd)

        outdir_tmp = outdir + samp[i] + "/"
        cmd = "mkdir " + outdir_tmp
        os.system(cmd)

        outfile1 = outdir_tmp + samp[i] + ".addgrp.bam"
        outfile2 = outdir_tmp + samp[i] + ".splited.bam"
        outfile3 = outdir_tmp + samp[i] + ".recal_data.table"
        outfile4 = outdir_tmp + samp[i] + ".bqsr.bam"
        outfile5 = outdir_tmp + samp[i] + ".AnalyzeCovariates.pdf"

        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J gatk_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %sgatk_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %sgatk_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")
         
        f.write("module load java\n")
        f.write("module load R\n")
        f.write("java -jar -Xmx10g \"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/picard.jar\" AddOrReplaceReadGroups I= %s O= %s  RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20   \n" % (files[i],outfile1))
        f.write("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk/gatk-4.2.5.0/./gatk SplitNCigarReads -R %s --input %s --output %s --tmp-dir %s\n" % (genome,outfile1,outfile2,workdir_tmp))
        f.write("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk/gatk-4.2.5.0/./gatk BaseRecalibrator -R %s --known-sites %s --input %s --output %s --tmp-dir %s\n" % (genome,knownsites,outfile2,outfile3,workdir_tmp))
        f.write("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk/gatk-4.2.5.0/./gatk ApplyBQSR -R %s --input %s --bqsr-recal-file %s --output %s --tmp-dir %s\n" % (genome, outfile2, outfile3,outfile4,workdir_tmp))
        f.write("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk/gatk-4.2.5.0/./gatk AnalyzeCovariates -bqsr %s -plots %s\n" % (outfile3,outfile5))
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

def run_variantcall(data_dir,workdir,outdir,genome):
    
    acc = "premium"
    time1 = "8:00"
    mem1 ="20000"
    
    cmd = "find "+data_dir+" -name \"*bqsr.bam\""
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    samp = []
    for i in range(len(files)):
        tmpx = files[i].split("/")[-1]
        tmp = tmpx.split(".bqsr.bam")[0]
        if tmp not in samp:
            samp.append(tmp)

    shfile = workdir + "variantcall.sh"
    
    
    fx = open(shfile,"w")
    for i in range(len(samp)-1):
        lsfile = workdir + samp[i] + "_varcall.lsf"


        workdir_tmp = workdir + samp[i] + "/"
        cmd = "mkdir " + workdir_tmp
        os.system(cmd)

        outfile=outdir+samp[i]+".caller.g.vcf.gz"
        

        

        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J vcall_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %svcall_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %svcall_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")

        f.write("module load java\n")
        f.write("module load R\n")
 
        f.write("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk/gatk-4.2.5.0/./gatk --java-options -Xmx10g HaplotypeCaller -R %s --input %s --output %s  -ERC GVCF --tmp-dir %s\n" % (genome,files[i],outfile,workdir_tmp))
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()
      
def run_variantfilter(data_dir,workdir,outdir,genome):
    acc = "premium"
    time1 = "8:00"
    mem1 ="20000"


    cmd = "find "+data_dir+" -name \"*.caller.g.vcf.gz\""
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    samp = []
    for i in range(len(files)):
        tmpx = files[i].split("/")[-1]
        tmp = tmpx.split(".caller.g.vcf.gz")[0]
        if tmp not in samp:
            samp.append(tmp)


    shfile = workdir + "variantfilter.sh"
    fx = open(shfile,"w")
    for i in range(len(samp)-1):
        lsfile = workdir + samp[i] + "_varfilter.lsf"
        
        workdir_tmp = workdir + samp[i] + "/"
        cmd = "mkdir " + workdir_tmp
        os.system(cmd)

        outfile=outdir+samp[i]+".filter.g.vcf.gz"
        
        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J vfilter_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %svfilter_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %sfilter_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")

        f.write("module load java\n")
        f.write("/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/gatk/gatk-4.2.5.0/./gatk VariantFiltration -R %s --filter-expression \"QD<2.0\" --filter-name \"QD2\" --filter-expression \"QUAL<30.0\" --filter-name \"QUAL30\" --filter-expression \"SOR>3.0\" --filter-name \"SOR3\" --filter-expression \"FS>60.0\" --filter-name \"FS60\"  --filter-expression \"MQ<40.0\" --filter-name \"MQ40\" --filter-expression \"MQRankSum<-12.5\" --filter-name \"MQRankSum-12.5\" --filter-expression \"ReadPosRankSum<-8.0\" --filter-name \"ReadPosRankSum-8\" -V %s -O %s  --tmp-dir  %s --verbosity ERROR \n" % (genome,files[i],outfile,workdir_tmp))
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

    
def run_pass(data_dir,workdir,outdir):
    
    acc = "premium"
    time1 = "8:00"
    mem1 ="20000"


    cmd = "find "+data_dir+" -name \"*.filter.g.vcf.gz\"" 
    
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    samp = []
    for i in range(len(files)):
        tmpx = files[i].split("/")[-1]
        tmp = tmpx.split(".filter.g.vcf.gz")[0]
        if tmp not in samp:
            samp.append(tmp)



    shfile = workdir + "variantpass.sh"
    fx = open(shfile,"w")
    for i in range(len(samp)-1):
        lsfile = workdir + samp[i] + "_varpass.lsf"
        
        outfile=outdir+samp[i]+".pass.g.vcf"
         
        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J vpass_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %svpass_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %spass_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n") 
        
        f.write("module load bcftools\n")

        f.write("bcftools view -f PASS %s > %s\n" % (files[i],outfile))
        f.write("bgzip -c %s > %s.gz\n" % (outfile,outfile))
        f.write("tabix -p vcf %s.gz\n" % outfile)
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()


def run_liftover(data_dir,workdir,outdir):        

    acc = "premium"
    time1 = "1:00"
    mem1 ="10000"
    
    cmd = "find "+data_dir+" -name \"*.pass.g.vcf.gz\""
    
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    samp = []
    for i in range(len(files)):
        tmpx = files[i].split("/")[-1]
        tmp = tmpx.split(".pass.g.vcf.gz")[0]
        if tmp not in samp:
            samp.append(tmp)

    shfile = workdir + "liftover.sh"
    fx = open(shfile,"w")
    for i in range(len(samp)-1):
        lsfile = workdir + samp[i] + "_lover.lsf"
        
        outfile=outdir+samp[i]+".liftover.g.vcf"
        outfile2=outdir+samp[i]+".reject_variant.vcf"
        
        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J lover_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %slover_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %slover_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")

        f.write("module load java\n")
        f.write("java -jar  -Xmx8g \"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/picard.jar\" LiftoverVcf I=%s O=%s CHAIN=/sc/arion/projects/zhuj05a/Wenhui/bladder/liftover/hg38ToHg19.over.chain.gz REJECT=%s R=/sc/arion/projects/zhuj05a/Wenhui/HBV/script/test_new_mosaik/VirusSeq/Mosaik_JumpDb/hg19.fa ALLOW_MISSING_FIELDS_IN_HEADER=true WARN_ON_MISSING_CONTIG=true\n" % (files[i],outfile,outfile2)) 
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()



if __name__=="__main__":
    data_dir = "/sc/arion/scratch/leee19/gtex-liver/ncbi/fastqfiles/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/filter/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/filter/"
    #run_fastq_filter(data_dir,outdir,workdir)


    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/filter/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/align/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/align/"
    genomedir = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/"
    #run_star(data_dir,outdir,workdir,genomedir)
    
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/align/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/picard/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/picard/"
    #run_picard(data_dir,workdir,outdir)
    
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/picard/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/gatk/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/gatk/"
    genome = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/GRCh38.primary_assembly.genome.fa"
    knownsites = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf"
    #run_gatk(data_dir,workdir,outdir,genome,knownsites)

    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/gatk/"
    workdir =  "/sc/arion/scratch/wangm08/liver_gtex/workdir/vcall/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vcall/"
    genome = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/GRCh38.primary_assembly.genome.fa"
    #run_variantcall(data_dir,workdir,outdir,genome)

    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vcall/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/vfilter/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vfilter/"
    genome = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/GRCh38.primary_assembly.genome.fa"
    #run_variantfilter(data_dir,workdir,outdir,genome)

    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vfilter/"
    workdir =  "/sc/arion/scratch/wangm08/liver_gtex/workdir/vpass/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vpass/"
    #run_pass(data_dir,workdir,outdir)
    
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vpass/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/liftover/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/liftover/"
    run_liftover(data_dir,workdir,outdir)

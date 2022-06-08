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
        outdir_tmp = outdir + samp[i] + "/"
        cmd = "mkdir " + outdir_tmp
        os.system(cmd)

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
        f.write("module load trim_galore\n")
        f.write("trim_galore --paired %s %s --no_report_file --length 36 --quality 20 -o %s\n" % (inputfile1,inputfile2,outdir_tmp))
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

    

def run_star(data_dir,outdir,workdir,genomedir):
     
    acc = "premium"
    time1 = "4:00"
    mem1 ="60000"
    nthr = 4  
    
    cmd = "find "+data_dir+" -name \"*.fq\""
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
        inputfile1 = data_dir + samp[i] + "/" + samp[i] + "_1_val_1.fq"
        inputfile2 = data_dir + samp[i] + "/" + samp[i] + "_2_val_2.fq"
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
        f.write("STAR  --runThreadN %d --genomeDir %s \\\n--readFilesIn %s %s \\\n--outFileNamePrefix %s \\\n--outSAMtype BAM SortedByCoordinate --outFilterType BySJout\n" % (nthr,genomedir,inputfile1,inputfile2,outputprefix))
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

def run_regenerate_genome(data_dir,genomedir):
    

    cmd = "find " + data_dir + " -name \"*.tab\""
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    

    stjdir = genomedir + "gatk_liver_star/stj/"
    for i in range(len(files)):
        cmd = "cp " + files[i] + " " + stjdir
        os.system(cmd)
    
    cmd = "cat " + stjdir + "*.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > " + stjdir + "SJ.filtered.tab"
    print(cmd)
    
    cmd = "STAR --runMode genomeGenerate --genomeDir " + genomedir + "gatk_liver_star/" + "--genomeFastaFiles " + genomedir + "GRCh38.primary_assembly.genome.fa --sjdbGTFfile " + genomedir + "gencode.v39.primary_assembly.annotation.gtf --runThreadN 8 --sjdbOverhang 75 --limitSjdbInsertNsj 2050000 --sjdbFileChrStartEnd " + stjdir + "SJ.filtered.tab"
    print(cmd)


def run_star_2pass(data_dir,workdir,outdir,genomedir):
  
    acc = "premium"
    time1 = "4:00"
    mem1 ="60000"
    nthr = 4

    cmd = "find "+data_dir+" -name \"*.fq\""
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

    shfile = workdir + "align2.sh"
    fx = open(shfile,"w")
    for i in range(len(samp)-1):
        lsfile = workdir + samp[i] + "_align2.lsf"
        inputfile1 = data_dir + samp[i] + "/" + samp[i] + "_1_val_1.fq"
        inputfile2 = data_dir + samp[i] + "/" + samp[i] + "_2_val_2.fq"
        outdir_tmp = outdir + samp[i] + "/"
        cmd = "mkdir " + outdir_tmp
        os.system(cmd)
        outputprefix = outdir_tmp + samp[i]

        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J align2_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %salign2_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %salign2_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")
        f.write("module load star\n")
        f.write("STAR  --runThreadN %d --genomeDir %s \\\n--readFilesIn %s %s \\\n--outFileNamePrefix %s \\\n--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM\n" % (nthr,genomedir,inputfile1,inputfile2,outputprefix))
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()
    




    
    

def run_picard(data_dir,workdir,outdir):
    
    acc = "premium"
    time1 = "4:00"
    mem1 ="20000"


    cmd = "find "+data_dir+" -name \"*Aligned.sortedByCoord.out.bam\""
    print(cmd)
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
        
        outfile0 = outdir_tmp + samp[i] + "Aligned.sortedByCoord.cs.bam"
        outfile1 = outdir_tmp + samp[i] + ".addgrp.bam"
        outfile2 = outdir_tmp + samp[i] + ".uniq.bam"
        outfile3 = outdir_tmp + samp[i] +  ".mkDup.bam"
        outfile4 = outdir_tmp + samp[i] + ".marked_dup_metrics.txt"
        
        
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
        f.write("java -jar -Xmx10g \"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/picard.jar\" CleanSam I= %s O= %s\n" % (files[i],outfile0))
        f.write("java -jar -Xmx10g \"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/picard.jar\" AddOrReplaceReadGroups I= %s O= %s RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20\n" % (outfile0,outfile1))
        f.write("rm %s\n" % outfile0)
        f.write("samtools view -h -q 255 %s -bo %s\n" % (outfile1, outfile2))
        f.write("rm %s\n" % outfile1)
        f.write("java -jar -Xmx10g \"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/picard.jar\" MarkDuplicates TMP_DIR= %s CREATE_INDEX= true VALIDATION_STRINGENCY=SILENT I=  %s O= %s M= %s\n" % (workdir_tmp,outfile2,outfile3,outfile4))
        f.write("rm %s\n" % outfile2)
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()


def run_gatk(data_dir,workdir,outdir,genome,knownsites):
    
    acc = "premium"
    time1 = "8:00"
    mem1 ="20000"



    cmd = "find "+data_dir+" -name \"*mkDup.bam\""
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    samp = []
    for i in range(len(files)):
        tmpx = files[i].split("/")[-1]
        tmp = tmpx.split(".mkDup.bam")[0]
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

        #outfile1 = outdir_tmp + samp[i] + ".addgrp.bam"
        outfile2 = outdir_tmp + samp[i] + ".splited.bam"
        outfile21 = outdir_tmp + samp[i] + ".real.intervals"
        outfile22 = outdir_tmp + samp[i] + ".real.bam"
        outfile3 = outdir_tmp + samp[i] + ".recal.table"
        outfile4 = outdir_tmp + samp[i] + ".recal.bam"

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
        #f.write("java -jar -Xmx10g \"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/picard.jar\" AddOrReplaceReadGroups I= %s O= %s  RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20   \n" % (files[i],outfile1))

        f.write("java -Xmx10g -jar /hpc/packages/minerva-common/gatk/3.5-0/src/GenomeAnalysisTK.jar -T SplitNCigarReads \\\n-R %s \\\n-I %s \\\n-o %s \\\n-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS\n\n" % (genome,files[i],outfile2))
        
        f.write("java -Xmx10g -jar /hpc/packages/minerva-common/gatk/3.5-0/src/GenomeAnalysisTK.jar -T RealignerTargetCreator \\\n")
        f.write("-R %s \\\n" % genome)
        f.write("-I %s \\\n" % outfile2)
        f.write("-known %s \\\n" % knownsites[1])
        f.write("-known %s \\\n" % knownsites[2])
        f.write("-o %s\n\n" % outfile21)
        
        f.write("java -Xmx10g -jar /hpc/packages/minerva-common/gatk/3.5-0/src/GenomeAnalysisTK.jar -T IndelRealigner \\\n")
        f.write("-I %s \\\n" % outfile2)
        f.write("-R %s \\\n" % genome)
        f.write("-targetIntervals %s \\\n" % outfile21)
        f.write("-known %s \\\n" % knownsites[1])
        f.write("-known %s \\\n" % knownsites[2])
        f.write("-o %s\n\n" % outfile22)
        
        f.write("java -Xmx10g -jar /hpc/packages/minerva-common/gatk/3.5-0/src/GenomeAnalysisTK.jar -T BaseRecalibrator \\\n")
        f.write("-S SILENT \\\n")
        f.write("-I %s \\\n" % outfile22)
        f.write("-R %s \\\n" % genome)
        f.write("-knownSites %s \\\n" % knownsites[0])
        f.write("-knownSites %s \\\n" % knownsites[1])
        f.write("-knownSites %s \\\n" % knownsites[2])
        f.write("-o %s\n\n" % outfile3)

        f.write("java -Xmx10g -jar /hpc/packages/minerva-common/gatk/3.5-0/src/GenomeAnalysisTK.jar -T PrintReads \\\n")
        f.write("-I %s \\\n" % outfile22)
        f.write("-R %s \\\n" % genome)
        f.write("-BQSR %s \\\n" % outfile3)
        f.write("-o %s\n\n" % outfile4)
        f.write("rm %s %s\n" % (outfile2,outfile22))
        

        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

def run_variantcall(data_dir,workdir,outdir,genome):
    
    acc = "premium"
    time1 = "8:00"
    mem1 ="20000"
    
    cmd = "find "+data_dir+" -name \"*recal.bam\""
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    samp = []
    for i in range(len(files)):
        tmpx = files[i].split("/")[-1]
        tmp = tmpx.split(".recal.bam")[0]
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
        outfile1=outdir+samp[i]+".caller.full.vcf.gz"

        

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
 
        f.write("java -Xmx10g -jar /hpc/packages/minerva-common/gatk/3.5-0/src/GenomeAnalysisTK.jar \\\n")
        f.write("-T HaplotypeCaller \\\n")
        f.write("-I %s \\\n" % files[i])
        f.write("-R %s \\\n" % genome)
        f.write("--min_base_quality_score 10 --min_mapping_quality_score 20 -ERC GVCF -dontUseSoftClippedBases \\\n")
        f.write("-o %s\n\n\n" % outfile)




        #f.write("java -Xmx10g -jar /hpc/packages/minerva-common/gatk/3.5-0/src/GenomeAnalysisTK.jar -T GenotypeGVCFs \\\n")
        #f.write("--variant %s \\\n" % outfile)
        #f.write("-R %s \\\n" % genome)
        #f.write("-stand_call_conf 20.0 \\\n")
        #f.write("-o %s\n" % outfile1)

        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

def run_rename_vcf(data_dir,outdir,workdir):
    acc = "premium"
    time1 = "8:00"
    mem1 ="20000"

    cmd = "find "+data_dir+" -name \"*caller.g.vcf.gz\""
    a = subprocess.check_output(cmd,shell=True).decode("utf-8")
    files = a.split("\n")
    samp = []
    for i in range(len(files)):
        tmpx = files[i].split("/")[-1]
        tmp = tmpx.split(".caller.g.vcf.gz")[0]
        if tmp not in samp:
            samp.append(tmp)

    shfile = workdir + "renamevcf.sh"


    fx = open(shfile,"w")
    for i in range(len(samp)-1):
        lsfile = workdir + samp[i] + "_renamevcf.lsf"
        outfile = outdir + samp[i] + ".caller.new.g.vcf.gz"       
 
        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J rnmv_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %srnmv_%s.stdout\n" % (workdir,samp[i]))
        f.write("#BSUB -eo %srnmv_%s.stderr\n" % (workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")

        f.write("module load java\n")
        f.write("java -jar /sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/picard.jar RenameSampleInVcf \\\n")
        f.write("INPUT=%s \\\n" % files[i])
        f.write("OUTPUT=%s \\\n" % outfile)
        f.write("NEW_SAMPLE_NAME=\"%s\" \\\n" % samp[i]) 
        f.write("CREATE_INDEX=true\n")
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

 

def run_combinegvcf(inputlist,workdir,outdir,genome):
    
    acc = "premium"
    time1 = "24:00"
    mem1 ="20000"
    
    chro = [] 
    for i in range(22):
        chro.append("chr"+str(i+1))
    
    chro.append("chrX")
    chro.append("chrY")
    chro.append("chrM")
    
    shfile = workdir + "combinegvcf.sh"
    fx = open(shfile,"w")

    for i in range(25):
        lsfile = workdir + chro[i] + "_combinegvcf.lsf"
        outfile = outdir + "liver_all_combine." + chro[i] + ".g.vcf.gz"        
 
        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J comb_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %scomb_%s.stdout\n" % (workdir,chro[i]))
        f.write("#BSUB -eo %scomb_%s.stderr\n" % (workdir,chro[i]))
        f.write("#BSUB -L /bin/bash\n")
        
        f.write("module load java\n")
        f.write("java -Xmx10g -jar /hpc/packages/minerva-common/gatk/3.5-0/src/GenomeAnalysisTK.jar -T CombineGVCFs \\\n")
        f.write("--variant %s \\\n" % inputlist)
        f.write("-R %s \\\n" % genome)
        f.write("-L %s \\\n" % chro[i])
        f.write("-o %s\n" % outfile)
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

def run_genotype(datadir,workdir,outdir,genome):
    acc = "premium"
    time1 = "24:00"
    mem1 ="20000"

    chro = []
    for i in range(22):
        chro.append("chr"+str(i+1))

    chro.append("chrX")
    chro.append("chrY")
    chro.append("chrM")

    shfile = workdir + "genotype.sh"
    fx = open(shfile,"w")

    for i in range(25):
        lsfile = workdir + chro[i] + "_genotype.lsf"
        inputfile = datadir + "liver_all_combine." + chro[i] + ".g.vcf.gz"
        outfile = outdir + "liver_all_genotype." + chro[i] + ".full.vcf.gz"

        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J geno_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %sgeno_%s.stdout\n" % (workdir,chro[i]))
        f.write("#BSUB -eo %sgeno_%s.stderr\n" % (workdir,chro[i]))
        f.write("#BSUB -L /bin/bash\n")

        f.write("module load java\n")
        f.write("java -Xmx10g -jar /hpc/packages/minerva-common/gatk/3.5-0/src/GenomeAnalysisTK.jar -T GenotypeGVCFs \\\n")
        f.write("--variant %s \\\n" % inputfile)
        f.write("-R %s \\\n" % genome)
        f.write("-stand_call_conf 20.0 \\\n")
        f.write("-o %s\n" % outfile)
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()
    


def run_variantfilter_chr(data_dir,workdir,outdir,genome):
    acc = "premium"
    time1 = "8:00"
    mem1 ="20000"


    chro = []
    for i in range(22):
        chro.append("chr"+str(i+1))

    chro.append("chrX")
    chro.append("chrY")
    chro.append("chrM")

    shfile = workdir + "filter_chr.sh"
    fx = open(shfile,"w")    

    fx = open(shfile,"w")

    for i in range(25):
        lsfile = workdir + chro[i] + "_fliterchr.lsf"
        outfile = outdir + "liver_all_combine." + chro[i] + ".filter.vcf.gz"
        outfile0 = outdir + "liver_all_combine." + chro[i] + ".select.vcf.gz"
        inputfile = data_dir + "liver_all_genotype." + chro[i] + ".full.vcf.gz"

        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J vfilter_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %svfilter_%s.stdout\n" % (workdir,chro[i]))
        f.write("#BSUB -eo %svfilter_%s.stderr\n" % (workdir,chro[i]))
        f.write("#BSUB -L /bin/bash\n")
        
        f.write("module load java\n")
        f.write("java -Xmx10g -jar /hpc/packages/minerva-common/gatk/3.5-0/src/GenomeAnalysisTK.jar -T SelectVariants \\\n")
        f.write("-R %s \\\n" % genome)
        f.write("-selectType SNP -selectType INDEL --restrictAllelesTo BIALLELIC \\\n")
        f.write("-V %s \\\n" % inputfile)
        f.write("-o %s\n\n\n" % outfile0)
        
        f.write("java -Xmx10g -jar /hpc/packages/minerva-common/gatk/3.5-0/src/GenomeAnalysisTK.jar -T VariantFiltration \\\n")
        f.write("-R %s \\\n" % genome)
        f.write("-window 35 -cluster 3 --filterExpression \"QD<2.0\" --filterName \"QD2\" --filterExpression \"FS>30.0\" --filterName \"FS30\" \\\n")
        f.write("-V %s \\\n" % outfile0)
        f.write("-o %s\n" % outfile)
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

def run_pass_chr(data_dir,workdir,outdir):
    acc = "premium"
    time1 = "8:00"
    mem1 ="20000"


    chro = []
    for i in range(22):
        chro.append("chr"+str(i+1))

    chro.append("chrX")
    chro.append("chrY")
    chro.append("chrM")

    shfile = workdir + "pass_chr.sh"
    fx = open(shfile,"w")

    fx = open(shfile,"w")

    for i in range(25):
        lsfile = workdir + chro[i] + "_passchr.lsf"
        inputfile = data_dir + "liver_all_combine." + chro[i] + ".filter.vcf.gz"
        outfile = outdir + "liver_all_combine." + chro[i] + ".pass.vcf"
    
        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J vpass_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %svpass_%s.stdout\n" % (workdir,chro[i]))
        f.write("#BSUB -eo %svpass_%s.stderr\n" % (workdir,chro[i]))
        f.write("#BSUB -L /bin/bash\n")

        f.write("module load bcftools\n")

        f.write("bcftools view -f PASS %s > %s\n" % (inputfile,outfile))
        f.write("bgzip -c %s > %s.gz\n" % (outfile,outfile))
        f.write("tabix -p vcf %s.gz\n" % outfile)
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

def run_liftover_chr(data_dir,workdir,outdir):
    acc = "premium"
    time1 = "8:00"
    mem1 ="20000"


    chro = []
    for i in range(22):
        chro.append("chr"+str(i+1))

    chro.append("chrX")
    chro.append("chrY")
    chro.append("chrM")

    shfile = workdir + "liftover_chr.sh"
    fx = open(shfile,"w")


    for i in range(25):
        lsfile = workdir + chro[i] + "_liftoverchr.lsf"
        inputfile = data_dir + "liver_all_combine." + chro[i] + ".pass.vcf.gz"
        outfile = outdir + "liver_all_combine." + chro[i] + ".liftover.vcf"
        outfile2 = outdir + "liver_all_combine." + chro[i] + ".rejectvariant.vcf"

        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J vlover_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %svlover_%s.stdout\n" % (workdir,chro[i]))
        f.write("#BSUB -eo %svlover_%s.stderr\n" % (workdir,chro[i]))
        f.write("#BSUB -L /bin/bash\n")

        f.write("module load java\n")
        f.write("java -jar  -Xmx8g \"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/picard.jar\" LiftoverVcf I=%s O=%s CHAIN=/sc/arion/projects/zhuj05a/Wenhui/bladder/liftover/hg38ToHg19.over.chain.gz REJECT=%s R=/sc/arion/projects/zhuj05a/Wenhui/HBV/script/test_new_mosaik/VirusSeq/Mosaik_JumpDb/hg19.fa ALLOW_MISSING_FIELDS_IN_HEADER=true WARN_ON_MISSING_CONTIG=true\n" % (inputfile,outfile,outfile2))
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

def run_rnavalid_chr(data_dir1,data_dir2,workdir,outdir):

    acc = "premium"
    time1 = "8:00"
    mem1 ="40000"

    chro = []
    for i in range(22):
        chro.append("chr"+str(i+1))

    chro.append("chrX")
    chro.append("chrY")
    chro.append("chrM")

    shfile = workdir + "rnavalid_chr.sh"
    fx = open(shfile,"w")
    
    for i in range(25):
        lsfile = workdir + chro[i] + "_rnavalidchr.lsf"
        inputfile1 = data_dir1 + "liver_all_combine." + chro[i] +".liftover.vcf"
        inputfile2 = data_dir2 + chro[i] + ".geno.vcf"
        outfile = outdir + chro[i] + ".rda"
         
        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J rvalid_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %srvalid_%s.stdout\n" % (workdir,chro[i]))
        f.write("#BSUB -eo %srvalid_%s.stderr\n" % (workdir,chro[i]))
        f.write("#BSUB -L /bin/bash\n")
        
        f.write("module load R\n")
        f.write("Rscript -e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/rnaValid.R\\\")\" \\\n")
        f.write("-e \"res_%s<-rnaValid(\\\"%s\\\",\\\"%s\\\")\" \\\n" % (chro[i],inputfile1,inputfile2))
        f.write("-e \"save(res_%s,file=\\\"%s\\\")\"\n" % (chro[i],outfile))
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

def run_rnavalid_chr2(data_dir1,data_dir2,workdir,outdir):

    acc = "premium"
    time1 = "8:00"
    mem1 ="40000"

    chro = []
    for i in range(22):
        chro.append("chr"+str(i+1))

    chro.append("chrX")
    chro.append("chrY")
    chro.append("chrM")

    shfile = workdir + "rnavalid_chr.sh"
    fx = open(shfile,"w")

    for i in range(25):
        lsfile = workdir + chro[i] + "_rnavalidchr.lsf"
        inputfile1 = data_dir1 + "liver_all_combine." + chro[i] +".liftover.vcf"
        inputfile2 = data_dir2 + chro[i].split("chr")[1] + ".geno.vcf"
        outfile = outdir + chro[i] + ".rda"
        
        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J rvalid_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %srvalid_%s.stdout\n" % (workdir,chro[i]))
        f.write("#BSUB -eo %srvalid_%s.stderr\n" % (workdir,chro[i]))
        f.write("#BSUB -L /bin/bash\n")


        f.write("module load R\n")
        f.write("Rscript -e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/compare_rna_genotype_jointcall.R\\\")\" \\\n")
        f.write("-e \"res_%s<-compare_rna_genotype_jointcall(\\\"%s\\\",\\\"%s\\\")\" \\\n" % (chro[i],inputfile1,inputfile2))
        f.write("-e \"save(res_%s,file=\\\"%s\\\")\"\n" % (chro[i],outfile))
        f.close()
        fx.write("bsub < %s\n" % lsfile)
    fx.close()

def extract_loci_chr(data_dir1,data_dir2,workdir,outdir):

    acc = "premium"
    time1 = "8:00"
    mem1 ="40000"

    chro = []
    for i in range(22):
        chro.append("chr"+str(i+1))

    chro.append("chrX")
    chro.append("chrY")
    chro.append("chrM")

    shfile = workdir + "extract_loci_chr.sh"
    fx = open(shfile,"w")

    for i in range(25):
        lsfile = workdir + chro[i] + "_extractlocichr.lsf"
        inputfile1 = data_dir1 + "liver_all_combine." + chro[i] +".liftover.vcf"
        inputfile2 = data_dir2 + chro[i].split("chr")[1] + ".geno.vcf"
        outfile = outdir + chro[i] + ".rda"

        f = open(lsfile,"w")
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -J extrloci_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %sextrloci_%s.stdout\n" % (workdir,chro[i]))
        f.write("#BSUB -eo %sextrloci_%s.stderr\n" % (workdir,chro[i]))
        f.write("#BSUB -L /bin/bash\n")

        f.write("module load R\n")
        f.write("Rscript -e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/bin/extract_loci.R\\\")\" \\\n")
        f.write("-e \"res_%s<-extract_loci(\\\"%s\\\",\\\"%s\\\")\" \\\n" % (chro[i],inputfile1,inputfile2))
        f.write("-e \"save(res_%s,file=\\\"%s\\\")\"\n" % (chro[i],outfile))
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

    #step1 raw reads filtering: cut reads with adapter seq and filter low quality reads
    data_dir = "/sc/arion/scratch/leee19/gtex-liver/ncbi/fastqfiles/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/filter/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/filter/"
    #run_fastq_filter(data_dir,outdir,workdir)

    #step2 first pass star alignment
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/filter/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/align/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/align/"
    genomedir = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/"
    #run_star(data_dir,outdir,workdir,genomedir)
    
    #step3 collect all the *tab files from all samples first pass align results  and filter tab and regenerated new 
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/align/"
    genomedir = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/"
    #run_regenerate_genome(data_dir,genomedir)
     
    #step4 2nd pass alignement
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/filter/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/align2/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/align2/"
    genomedir = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/gatk_liver_star/"
    #run_star_2pass(data_dir,workdir,outdir,genomedir)
    
    #step5 clean bam file, add group, filter multiple mapping reads and markduplicate
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/align2/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/picard/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/picard/"
    #run_picard(data_dir,workdir,outdir)
    
    #step6 SplitNCigarReads, indel realignment and base quality recalibration
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/picard/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/gatk/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/gatk/"
    genome = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/GRCh38.primary_assembly.genome.fa"
    knownsites = ["/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf","/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/Homo_sapiens_assembly38.known_indels.vcf.gz","/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"]
    #run_gatk(data_dir,workdir,outdir,genome,knownsites)

    #step7 haplotype call with gvcf
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/gatk/"
    workdir =  "/sc/arion/scratch/wangm08/liver_gtex/workdir/vcall/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vcall/"
    genome = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/GRCh38.primary_assembly.genome.fa"
    #run_variantcall(data_dir,workdir,outdir,genome)
    
    #step8 rename the vcf file (defualt one with the same sample name)
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vcall/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vrname/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/vrname/"
    #run_rename_vcf(data_dir,outdir,workdir)

    #find /sc/arion/scratch/wangm08/liver_gtex/outdir/vrname -type f -name "*.g.vcf.gz" > input.list
    #step9 combine the gvcfs from samples for each chromosome (cutdown time)
    inputlist = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vrname/input.list"
    workdir =  "/sc/arion/scratch/wangm08/liver_gtex/workdir/combinegvcf/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/combinegvcf/"
    genome = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/GRCh38.primary_assembly.genome.fa"
    #run_combinegvcf(inputlist,workdir,outdir,genome)
    
    #step10 joint genotype for each chrom
    datadir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/combinegvcf/"
    workdir =  "/sc/arion/scratch/wangm08/liver_gtex/workdir/genotype/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/genotype/"
    genome = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/GRCh38.primary_assembly.genome.fa"
    #run_genotype(datadir,workdir,outdir,genome)

    #step11 select SNPs, indels and bialliec SNVs and filter the SNVs
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/genotype/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/filter_chr/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/filter_chr/"
    genome = "/sc/arion/projects/zhuj05a/Wenhui/IBD/data/GRCh38/GRCh38.primary_assembly.genome.fa"
    #run_variantfilter_chr(data_dir,workdir,outdir,genome)

    #step 12 keep only SNVs pass the filter
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/filter_chr/"
    workdir =  "/sc/arion/scratch/wangm08/liver_gtex/workdir/pass_chr/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/pass_chr/"
    #run_pass_chr(data_dir,workdir,outdir)

    #step 13 liftover to hg19
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/pass_chr/"
    workdir =  "/sc/arion/scratch/wangm08/liver_gtex/workdir/liftover_chr/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/liftover_chr/"
    #run_liftover_chr(data_dir,workdir,outdir)
    
    #step 14 validate the RNAseq results with genotype from wgs
    data_dir1 = "/sc/arion/scratch/wangm08/liver_gtex/outdir/liftover_chr/"
    data_dir2 = "/sc/arion/scratch/wangm08/liver_gtex/genotype/outdir/bychr/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/rnavalid_chr/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/rnavalid_chr/"
    #run_rnavalid_chr(data_dir1,data_dir2,workdir,outdir)

    data_dir1 = "/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/liftover_chr/"
    data_dir2 = "/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/genotype_chr/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/rnavalid_chr/"
    outdir = "/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/valid_rna_geno/"
    #run_rnavalid_chr2(data_dir1,data_dir2,workdir,outdir)

    data_dir1 = "/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/liftover_chr/"
    data_dir2 = "/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/genotype_chr/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/extract_loci/"
    outdir = "/sc/arion/projects/zhuj05a/Wenhui/IBD/colon/data/gatk_somatic_snv/liver_GTEx_rnaseq_snp/rna_seq_joint_call/extract_loci/"
    extract_loci_chr(data_dir1,data_dir2,workdir,outdir)
 
#####################################################################################################################################################
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vfilter/"
    workdir =  "/sc/arion/scratch/wangm08/liver_gtex/workdir/vpass/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vpass/"
    #run_pass(data_dir,workdir,outdir)
    
    data_dir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/vpass/"
    workdir = "/sc/arion/scratch/wangm08/liver_gtex/workdir/liftover/"
    outdir = "/sc/arion/scratch/wangm08/liver_gtex/outdir/liftover/"
    #run_liftover(data_dir,workdir,outdir)

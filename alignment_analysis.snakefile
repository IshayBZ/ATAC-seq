# Purpose: align paired-end data to produce bam files and then extract fragment ends into txt files
# Input:
#   fastq
# Output:
#   bam and txt

# TODO
# 1) Specify working directory (with snakefile and python codes)
# 2) Specify bowtie2 params with path to Bowtie2 reference
# 3) Specify "if" statement that defines which refernece to choose for each sample
# UPDATE - python code will make txt files for all coordinates (inc HO locus) so samples aligned to YeastGenome will have HO txt files

# To run:
# 1) cd to the folder where the snakefile and Coord-extraction python codes are
# 2) snakemake --snakefile filename.snakefile -j 25 --cluster "sbatch -c 4 -n 1 -p short --mem 4G -t 0-2:00" --latency-wait 60
#    (this will skip steps that have already been done and their output files exist)
# or: same with -F
#    (this will not skip steps that have already been done and overwrite their output files)

# Working directory (where sankefiles and py codes are)
workdir: "/n/groups/springer/Ishay/O3/fastq"

import re
import glob

ref_base = "/n/groups/springer/Ishay/test_reference/YeastGenome"
ref_HO = "/n/groups/springer/Ishay/test_reference/YeastGenome_HO_SNPJ"
ref_HO_TTDEL = "/n/groups/springer/Ishay/test_reference/YeastGenome_HO_SNPJ_TTDEL"

# list(range(1,n)) creates a list from 1 to n-1
list_base = list(range(17, 29)) + list(range(33, 43))
list_HO = list(range(29, 33)) + list(range(43, 63)) + list(range(71, 76)) + list(range(78, 79)) + list(range(84, 86)) + list(range(87, 95)) + list(range(97, 101)) + list(range(103, 107)) + list(range(109, 113))
# don't comment list_HO, just put irrelevant numbers
# list_HO_TTDEL = list(range(1, 17))        # UNEEDED - USE IF-ELSE INSTEAD

SAMPLES = []
fastq_files = glob.glob("*.fastq.gz")
for i in fastq_files:
    SAMPLES.append(i.split('_R')[0])
# this takes the file name up to '_R' as the sample name - this usually looks like 'LIB051037_GEN00209659_S1' but sometimes just 'S1'


rule all:
# This is the target rule - dictates which files to make and therefore which rules to run
    input:
#	expand("trimmed/{sample}_R1_001_paired.fastq.gz", sample = SAMPLES),
        expand("sorted/{sample}_Yeast_sorted.bam", sample = SAMPLES),
	expand("filtered/{sample}_Yeast_sorted_filtered.bam", sample = SAMPLES),
        "coords_extracted.done"

#rule trimmomatic:
#    input:
#        R1="{sample}_R1_001.fastq.gz",
#        R2="{sample}_R2_001.fastq.gz"
#    params:
#        "DNE.fa"
#    output:
#        R1_pair="trimmed/{sample}_R1_001_paired.fastq.gz",
#        R1_unpair=temp("trimmed/{sample}_R1_001_fail.fastq.gz"),
#        R2_pair="trimmed/{sample}_R2_001_paired.fastq.gz",
#        R2_unpair=temp("trimmed/{sample}_R2_001_fail.fastq.gz")
#    log:
#        "trimmed/log_files/{sample}_trim.log"
#    shell:  
#        """
#        (trimmomatic PE -threads 4 {input.R1} {input.R2} {output.R1_pair} {output.R1_unpair} {output.R2_pair} {output.R2_unpair} \
#            ILLUMINACLIP:{params}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60) 2> {log}
#        """

rule bowtie2:
    input:
        R1_pair="{sample}_R1_001.fastq.gz",
        R2_pair="{sample}_R2_001.fastq.gz"
#  	R1_pair="trimmed/{sample}_R1_001_paired.fastq.gz", # in case we start with a trimming step
#       R2_pair="trimmed/{sample}_R2_001_paired.fastq.gz"
    params:
# this just takes whatever is before 'S' in the sample name
# updated on 22/1/25 to split by 'S' and not '_S' since some data sets had file names that started with S and number wo the long beginning
        ref_genome = lambda wildcards : ref_base if (int(wildcards.sample.split('S')[1]) in list_base) else (ref_HO if (int(wildcards.sample.split('S')[1]) in list_HO) else ref_HO_TTDEL),
#	    ref_genome = lambda wildcards : ref_HO if (int(wildcards.sample.split('S')[1]) in list_base) else ref_HO_TTDEL,
	    ref_log = "sorted/{sample}_which_genome.log"
    output:
        sam=temp("temp/{sample}_Yeast.sam"),
        bam=temp("temp/{sample}_Yeast.bam"),
        sort="sorted/{sample}_Yeast_sorted.bam",
        index="sorted/{sample}_Yeast_sorted.bam.bai",
        idx="sorted/{sample}_Yeast_counts.tsv"
    log:
        "sorted/log_files/{sample}_Align_stats.log",
    shell:
        """
	echo {params.ref_genome} > {params.ref_log}
        (bowtie2 -x {params.ref_genome} --omit-sec-seq --no-unal --very-sensitive --maxins 1000 --no-discordant -p 4 \
            -1 {input.R1_pair} -2 {input.R2_pair} -S {output.sam}) 2> {log}
        samtools view -bS {output.sam} > {output.bam}
        samtools sort {output.bam} -o {output.sort}
        samtools index {output.sort} {output.index}
        samtools idxstats {output.sort} &> {output.idx}
        """

rule filter_bam: 		# Filter proper pairs (aligned concordantly)
    input:
        "sorted/{sample}_Yeast_sorted.bam"
    output:
        "filtered/{sample}_Yeast_sorted_filtered.bam"
    shell:
        """
        samtools view -h -f 3 {input} -o {output}
        samtools index {output}
        """

rule read_coords_extract:
    input:
        expand("filtered/{sample}_Yeast_sorted_filtered.bam", sample = SAMPLES)
    output:
        touch("coords_extracted.done")
    shell:
        """
        mkdir -p filtered/txt # make dir but not overwrite if it exists
	mkdir -p filtered/txt/AllChr
        python3 read_coords_extract.py
        """

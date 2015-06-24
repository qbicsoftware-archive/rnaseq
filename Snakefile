import os
import sys
import subprocess
import tempfile
import uuid
import shutil
import jinja2
from datetime import datetime
from os.path import join as pjoin
from os.path import exists as pexists
import hashlib

workdir: "../../../../var/marius_wf/rna_seq/TopHat2_all"

INPUT_FILES = []
for name in os.listdir('../../../../data'):
    if name.lower().endswith('.fastq'):
        if not name.endswith('.fastq'):
            print("Extension fastq is case sensitive.", file=sys.stderr)
            exit(1)
        INPUT_FILES.append(os.path.basename(name)[:-6])

OUTPUT_FILES = ["Summary/NumReads/Original/{name}.txt".format(name=name) for name in INPUT_FILES],\
               ["Summary/NumReads/PreFilter/{name}.txt".format(name=name) for name in INPUT_FILES],\
               ["FastQC/{name}".format(name=name) for name in INPUT_FILES],\
               ["Summary/NumReads/CutAdaptMerge/{name}.txt".format(name=name) for name in INPUT_FILES]

rule all:
    #input: ["CutAdaptMerge/{name}.fastq".format(name=name) for name in INPUT_FILES], OUTPUT_FILES
    #input: expand("HTSeqCounts/{name}.txt FastQC/{name}".split(), name=INPUT_FILES
    input: ["HTSeqCounts/{name}.txt".format(name=name) for name in INPUT_FILES], OUTPUT_FILES

rule PreFilterReads:
    input: "../../../../data/{name}.fastq"
    output: "PreFilterReads/{name}.fastq"
    shell: 'grep "1:N:" --no-group-separator -A 3 {input} > {output}'
    
rule FastQC:
    input: "PreFilterReads/{name}.fastq"
    output: "FastQC/{name}","FastQC/{name}/{name}_fastqc"
    shell: 'mkdir -p {output} && (fastqc {input} -o {output} || (rm -rf {output} && exit 1))'

rule FastQCcut:
    input: "CutAdaptMerge/{name}.fastq"
    output: "FastQCcut/{name}"
    shell: 'mkdir -p {output} && (fastqc {input} -o {output} || (rm -rf {output} && exit 1))'

rule Overrepresented:
    input: "FastQC/{name}/{name}_fastqc/fastqc_data.txt"
    output: "Overrepresented/{name}.txt"
    run:
        f = open(str(input))
        out = open(str(output), "w")
        sw = False
        for line in f:
            if line.startswith(">>Overrepresented"):
                sw = True
            if not line.startswith(">>END_MODULE") and sw:
                out.write(line)
            else:
                sw = False
        f.close()
        out.close()

rule OverrepTxtFasta:
    input: "Overrepresented/{name}.txt"
    output: "Overrepresented/{name}.fasta"
    run:
        f = open(str(input))
        out = open(str(output), "w")
        for line in f:
            if not (line.startswith("#") or line.startswith(">>")):
                 tokens1 = line.split("\t")
                 print(tokens1)
                 fseq = tokens1[0]
                 fheader = '>' + tokens1[3]
                 out.write(fheader)#+'\n')#it just happens that is the last column
                 out.write(fseq+'\n')
        f.close()
        out.close()

rule MergeAdapters:
    input: expand("Overrepresented/{name}.fasta", name=INPUT_FILES)
    output: "MergeAdapters/merged.fasta"
    shell: "cat {input} > {output}"

rule CutAdapt:
    input: "MergeAdapters/merged.fasta", "PreFilterReads/{name}.fastq" 
    output: "CutAdaptMerge/{name}.fastq"
    shell: 'cutadapt --discard-trimmed -a file:{input[0]} -o {output} {input[1]}'

rule TopHat2:
    input: "CutAdaptMerge/{name}.fastq"
    output: "TopHat2/{name}"
    shell: 'tophat --no-coverage-search -o {output} -p 2 -G ../../../../ref/UCSC_mm10/annotation/genes.gtf ../../../../ref/UCSC_mm10/Sequence/Bowtie2Index/genome {input}'

rule HTSeqCounts:
    input: "TopHat2/{name}"
    output: "HTSeqCounts/{name}.txt"
    shell: "samtools view {input}/accepted_hits.bam | htseq-count -i gene_id -t exon -s yes - ../../../../ref/UCSC_mm10/annotation/genes.gtf > {output}"

rule IndexBAM:
    input: "TopHat2/{name}"
    output: "TopHat2/{name}/accepted_hits.bai"
    shell: "samtools index {input}/accepted_hits.bam  {output}"

rule PerBaseCoverage:
    input: "TopHat2/{name}"
    output: "Statistics/{name}_PerBaseCoverage.txt"
    shell: "samtools depth {input}/accepted_hits.bam > {output}"

rule Numreads:
    input: "PreFilterReads/{name}.fastq"
    output: "Summary/NumReads/PreFilter/{name}.txt"
    shell: 'wc -l {input} > {output}'

rule NumreadsCut:
    input: "CutAdaptMerge/{name}.fastq"
    output: "Summary/NumReads/CutAdaptMerge/{name}.txt"
    shell: 'wc -l {input} > {output}'

rule NumreadsOrig:
    input: "../../../../data/{name}.fastq"
    output: "Summary/NumReads/Original/{name}.txt"
    shell: 'wc -l {input} > {output}'


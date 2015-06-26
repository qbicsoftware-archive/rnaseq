import os
import sys
import subprocess
from os.path import join as pjoin
from os.path import exists as pexists

configfile: "config.json"
workdir: config["var"]


DATA = config['data']
RESULT = config['result']
LOGS = config['logs']
REF = config['ref']
INI_PATH = config['etc']
SNAKEDIR = config['src']
params = config['params']

def data(path):
    return os.path.join(DATA, path)

def ref(path):
    return os.path.join(REF, path)

def log(path):
    return os.path.join(LOGS, path)

def result(path):
    return os.path.join(RESULT, path)


gtf = os.path.join(REF, params["gtf"])
if not gtf:
    raise ValueError("no gtf file supplied.")

indexedGenome = os.path.join(config['ref'], params["indexedGenome"])
if not indexedGenome:
    raise ValueError("no indexedGenome file supplied.")


INPUT_FILES = []
for name in os.listdir(DATA):
    if name.lower().endswith('.fastq'):
        if not name.endswith('.fastq'):
            print("Extension fastq is case sensitive.", file=sys.stderr)
            exit(1)
        INPUT_FILES.append(os.path.basename(name)[:-6])
    elif name.lower().endswith('.fastq.gz'):
        if not name.endswith('.fastq.gz'):
            print("Extension fastq is case sensitive.", file=sys.stderr)
            exit(1)
        INPUT_FILES.append(os.path.basename(name)[:-len('.fastq.gz')])
    else:
        print("Unknown data file: %s" % name)
        exit(1)

if len(set(INPUT_FILES)) != len(INPUT_FILES):
    print("Some input file names are not unique")
    exit(1)

OUTPUT_FILES = []

OUTPUT_FILES.extend(expand("Summary/NumReads/Original/{name}.txt", name=INPUT_FILES, result=RESULT))
OUTPUT_FILES.extend(expand("Summary/NumReads/PreFilter/{name}.txt", name=INPUT_FILES, result=RESULT))
OUTPUT_FILES.extend(expand("{result}/FastQC_{name}.zip", name=INPUT_FILES, result=RESULT))
OUTPUT_FILES.extend(expand("Summary/NumReads/CutAdaptMerge/{name}.txt", name=INPUT_FILES, result=RESULT))
OUTPUT_FILES.extend(expand("{result}/HTSeqCounts_{name}.txt", name=INPUT_FILES, result=RESULT))
OUTPUT_FILES.extend(expand("TopHat2/{name}/accepted_hits.bai", name=INPUT_FILES, result=RESULT))
OUTPUT_FILES.extend(expand("Summary/MappingStats/{name}.txt", name=INPUT_FILES, result=RESULT))


rule all:
    input: OUTPUT_FILES

rule LinkUncompressed:
    input: data("{name}.fastq")
    output: "fastq/{name}.fastq"
    shell: "ln -s {input} {output}"

rule Uncompress:
    input: data("{name}.fastq.gz")
    output: "fastq/{name}.fastq"
    shell: "zcat {input} > {output}"

rule PreFilterReads:
    input: "fastq/{name}.fastq"
    output: "PreFilterReads/{name}.fastq"
    shell: 'grep "1:N:" --no-group-separator -A 3 {input} > {output}'

rule FastQC:
    input: "PreFilterReads/{name}.fastq"
    output: "FastQC/{name}"
    shell: 'mkdir -p {output} && (fastqc {input} -o {output} || (rm -rf {output} && exit 1))'

rule FastQCCpToResult:
    input: "FastQC/{name}"
    output: result("FastQC_{name}.zip")
    shell: "cp {input}/{wildcards.name}_fastqc.zip {output}"

rule FastQCcut:
    input: "CutAdaptMerge/{name}.fastq"
    output: "FastQCcut/{name}"
    shell: 'mkdir -p {output} && (fastqc {input} -o {output} || (rm -rf {output} && exit 1))'

rule Overrepresented:
    input: "FastQC/{name}"
    output: "Overrepresented/{name}.txt"
    run:
        f = open(str(input) + "/" + wildcards['name'] + "_fastqc/fastqc_data.txt")
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
    #shell: 'tophat --no-coverage-search -o {output} -p 2 -G ' + os.path.join(REF, params["gtf"]) + ' ' + os.path.join(REF, params["indexedGenome"]) + ' {input}'
    shell: 'tophat --no-coverage-search -o {output} -p 2 -G ' + gtf + ' ' + indexedGenome + ' {input}'

rule HTSeqCounts:
    input: "TopHat2/{name}"
    output: result("HTSeqCounts_{name}.txt")
    shell: "samtools view {input}/accepted_hits.bam | htseq-count -i gene_id -t exon -s yes - " + gtf + "  > {output}"

rule IndexBAM:
    input: "TopHat2/{name}"
    output: "TopHat2/{name}/accepted_hits.bai"
    shell: "samtools index {input}/accepted_hits.bam && mv -f {input}/accepted_hits.bam.bai  {output}"

rule CpAlignSummary:
    input: "TopHat2/{name}"
    output: "Summary/MappingStats/{name}.txt"
    shell: "cp {input}/align_summary.txt {output}"

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
    input: "fastq/{name}.fastq"
    output: "Summary/NumReads/Original/{name}.txt"
    shell: 'wc -l {input} > {output}'


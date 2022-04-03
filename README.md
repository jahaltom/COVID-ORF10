# COVID-19-Quantification
This pipeline will take Human RNA-Seq data and concurrently quantify the expression on the Human (GencodeV36) and SARS-COV-2(ASM985889v3) transcriptomes as well as human evidence based (EB) gene treanscripts using pyrpipe https://github.com/urmi-21/pyrpipe.

## Prerequisites
* Run: `conda env create -f env.yaml`
* Activate the environment: `conda activate pyrpipe_covid`


## Snakefiles and filters
* Bam-Fastq_Quant.py: Takes in bam/fastq files and will quantify the expression on Human and SARS-COV-2 transcriptome as well as EB transcripts. Outputs TPM and counts at gene and transcript level. Gene level summed up across transcripts. 

  * Specify FileType and Layout in config.yaml.
  * ids.txt must contain directory names for individual sample fastq(s). Below would be (Sample1,Sample2) <-one per line in ids.txt. 
  * Bam/Fastq files should be in this structre:
```
out/Sample1/Sample1.fastq
out/Sample2/Sample2.fastq
```

* SingleStudy_Filter.py: This takes output from Bam-Fastq_Quant.py (results_TPM_gene.tsv and results_Count_gene.tsv) and filters it by removing EB genes where the median TPM < 1.  


## Execution
```
snakemake -j 50 -s Bam-Fastq_Quant.py --cluster "sbatch -t 01:00:00 -c 25"
```

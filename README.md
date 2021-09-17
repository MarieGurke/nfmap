# Nfmap - nextflow mapping pipeline

**UNDER CONSTRUCTION!**

This nextflow pipeline maps paired and merged reads from several samples to a reference genome. 

Dependencies are samtools and bwa-mem2. The steps included in this pipeline are: 

1. Indexing the reference genome using samtools faidx and bwa-mem2 index
2. Mapping two paired end and one single end merged reads file to the reference genome using bwa-mem2
3. Quality filtering at a mapping quality score of 30 (default) and sorting of the reads using samtools view and sort
4. Merging the two filtered and sorted .bam files into only one file per sample using samtools merge
5. Collecting stats about the mapped reads (samtools idxstats and depth) and creating an overview sheet containing information about all samples



## Usage	



```bash
nextflow run main.nf --read_dir ../data/reads/ --ref_dir /data/reference/ --out_dir /data/output/ -profile mfn
```

Parameters: 

* --read_dir 	:		Directory containing subfolders for each samples. In those subfolders must be three \*.fastq.gz files. Two containing paired end reads (\*R{1,2}.fastq.gz) and one containing merged reads (*U.fastq.gz). 
* --ref_dir 	   :	    Directory containing the reference genome in *.fasta format. 
* --out_dir       :        Directory where the pipeline ouptut should be written into. 

* -profile         :        Can be either local (for running it on a local computer) or mfn for running it on the museums cluster. 


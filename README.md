# Nfmap - nextflow mapping pipeline

This nextflow pipeline maps paired and merged reads from several samples to a reference genome. Currently adapted to work with the output of this [nf-polish](https://github.com/MozesBlom/nf-polish) pipeline on the MfN cluster.

Dependencies are samtools and bwa mem. The steps included in this pipeline are:

1. Indexing the reference genome using bwa index
2. Mapping two paired end and one single end file with merged reads to the reference genome using bwa mem
3. Quality filtering at a mapping quality score of 30 (default) and sorting of the reads using samtools view and sort
4. Merging the two filtered and sorted .bam files into only one file per sample using samtools merge
5. Collecting stats about the mapped reads (samtools idxstats and depth) and creating an overview sheet containing information about all samples



## Usage



```bash
nextflow run main.nf --read_dir ../data/reads/ --ref /data/reference/ref_file.fa --out_dir /data/output/ -profile mfn
```

**Parameters:**

* --read_dir 	    :		Directory containing sub folders for each samples. In those sub folders must be three \*.fastq.gz files. Two containing paired end reads (\*R{1,2}.fastq.gz) and one containing merged reads (*U.fastq.gz).
* --ref 	        :	    Path to reference file in .fasta format. 
* --out_dir            :        Directory where the pipeline output should be written into.
* --qual_fil            :        (INT) mapping quality score cut of for the quality filter steps (Default: 30)
* -profile              :        Can be either local (for running it on a local computer) or mfn for running it on the museums cluster.
* -resume            :        Resume analysis if it was interrupted

* --with-report    :        Make nextflow create a report about the resources the job used.

**Output:**

Output folder will be filled with folders for each sample that contains one *.bam file. Additionally, there will be a file in the output folder called Mapping_stats.tsv. It contains the number of reads mapped, total bp mapped, average depth and the coverage of each samples in a tab separated table.

#!/usr/bin/env nextflow
log.info """\
         Fancy nextflow mapping script
         ===================================
         Read directory          : ${params.read_dir}
         Reference directory     : ${params.ref_dir}
         Output directory        : ${params.out_dir}
         Mapping quality cutoff  : ${params.qual_fil}
         """
         .stripIndent()

// holds the reference file, goes into index_ref
Channel
  .fromPath("${params.ref_dir}*.fasta*")
  .into{ref_index_ch;ref_len_ch}

// holds the forward and reverse reads + the sample id, goes into map_PE
Channel
  .fromFilePairs("${params.read_dir}*/*_R{1,2}.fastq.gz")
  .set{readPE_ch}

// holds the merged reads, goes into map_SE
Channel
  .fromPath("${params.read_dir}*/*_U.fastq.gz")
  .set{readSE_ch}

read_ch = readPE_ch.merge(readSE_ch)


process index_ref {

  label 'RAM_high'

  input:
  file(ref) from ref_index_ch

  output:
  tuple file(ref), file('*') into ref_mapPE_ch
  tuple file(ref), file('*') into ref_mapSE_ch

  script:
  """
    bwa-mem2 index $ref
    samtools faidx $ref
  """
}

// combine the ref and read channel to use the same reference for all reads
// (This is basically the replacemnt for "each tuple", which nf doesn't allow)
ref_and_reads = ref_mapPE_ch.combine(read_ch)

process mapping {

  label 'RAM_high'

  input:
  tuple file(ref), file('*'), val(sample_id), file(reads_PE), file(reads_SE) from ref_and_reads

  output:
  tuple val(sample_id), file('*_PE.sam'), file('*_SE.sam') into mapped_ch
  val(sample_id) into id_channel

  script:
  """
    bwa-mem2 mem -t ${task.cpus} $ref $reads_PE > ${sample_id}_PE.sam
    bwa-mem2 mem -t ${task.cpus} $ref $reads_SE > ${reads_SE}_SE.sam
  """

}

process filtering {

  label 'RAM_high'

  input:
  tuple val(sample_id), file(samPE),file(samSE) from mapped_ch

  output:
  tuple val(sample_id), file('*_PE.bam'),file('*_SE.bam') into filtered_ch

  script:
  """
    samtools view -@ ${task.cpus} -q ${params.qual_fil} -Su $samPE | samtools sort -@ ${task.cpus} > ${samPE.baseName}.bam
    samtools view -@ ${task.cpus} -q ${params.qual_fil} -Su $samSE | samtools sort -@ ${task.cpus} > ${samSE.baseName}.bam

  """
  // samtool view command:
  // -q mapping quality cut-off
  // -S sam input file
  // -u uncompressed bam output file
}

process merge_bams {
  publishDir "$params.out_dir/$sample_id", mode: 'copy'

  input:
  tuple val(sample_id), file(pe),file(se) from filtered_ch

  output:
  file('*.bam') into merged_bam

  script:
  """
    samtools merge -@ ${task.cpus} ${sample_id}.bam $se $pe
  """
  // -@ Number of threads
}

process stats_file {
  publishDir "$params.out_dir", mode: 'copy'

  output:
  file('Mapping_stats.tsv') into stats_file_ch

  script:
  """
  echo -e \
   "Sample\tNumber of reads mapped\tTotal bp mapped\tAverage depth\tCoverage" \
   > Mapping_stats.tsv
  """
}

process stats {
  //publishDir "$params.out_dir/${bam.baseName}", pattern: '*.txt', mode: 'copy'
  publishDir "$params.out_dir", pattern: '*.tsv', mode: 'copy'
  input:
  file(stats_file) from stats_file_ch
  each file(bam) from merged_bam

  output:
  file('*_idxstats.txt')
  file('*_depth.txt')
  file(stats_file)

  script:
  """
    samtools index -@ ${task.cpus} $bam
    samtools idxstats $bam > ${bam.baseName}_idxstats.txt
    samtools depth $bam > ${bam.baseName}_depth.txt
    re_map=\$(awk '{if(NR==1) print \$3}' ${bam.baseName}_idxstats.txt)
    tot_bp=\$(awk '{sum+=\$3}END{print sum}' ${bam.baseName}_depth.txt)
    depth=\$(awk '{sum+=\$3;cnt++}END{print sum/cnt}' ${bam.baseName}_depth.txt| sed 's/,/./')
    ref_len=\$(awk '{if(NR==1) print \$2}' ${bam.baseName}_idxstats.txt)
    cov=\$(echo "scale=2 ; \$tot_bp / \$ref_len" | bc)
    echo -e "${bam.baseName}\t\$re_map\t\$tot_bp\t\$depth\t\$cov" >> $stats_file
  """
}

#!/usr/bin/env nextflow
log.info """\
         Fancy nextflow mapping script
         ===================================
         read_dir        : ${params.read_dir}
         ref_dir         : ${params.ref_dir}
         out_dir         : ${params.out_dir}
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

process index_ref {

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
ref_and_reads = ref_mapPE_ch.combine(readPE_ch)

process map_PE {
  input:
  tuple file(ref), file('*'), val(sample_id), file(reads) from ref_and_reads

  output:
  file('*.sam') into mapped_PE
  val(sample_id) into id_channel

  script:
  """
    bwa-mem2 mem $ref $reads > ${sample_id}_PE.sam
  """

}

process map_SE {
  input:
  tuple file(ref), file('*') from ref_mapSE_ch
  each file(reads) from readSE_ch

  output:
  file('*.sam') into mapped_SE

  script:
  """
    bwa-mem2 mem $ref $reads > ${reads}_SE.sam
  """
}

process filter_PE {
  input:
  file(sam) from mapped_PE

  output:
  file('*.bam') into filtered_PE

  script:
  """
    samtools view -q 30 -Su $sam | samtools sort > ${sam.baseName}.bam

  """
  // samtool view command:
  // -q mapping quality cut-off
  // -S sam input file
  // -u uncompressed bam output file
}

process filter_SE {
  input:
  file(sam) from mapped_SE

  output:
  file('*.bam') into filtered_SE

  script:
  """
    samtools view -q 30 -Su $sam | samtools sort > ${sam.baseName}.bam

  """
  // samtool view command:
  // -q mapping quality cut-off
  // -S sam input file
  // -u uncompressed bam output file
}

process merge_bams {
  publishDir "$params.out_dir/$sample_id", mode: 'copy'

  input:
  file(se) from filtered_SE
  file(pe) from filtered_PE
  val(sample_id) from id_channel

  output:
  file('*.bam') into merged_bam

  script:
  """
    samtools merge -@ 2 ${sample_id}.bam $se $pe
  """
  // -@ Number of threads
}

process stats_file {
  publishDir "$params.out_dir", mode: 'copy'

  output:
  file('Mapping_stats.tsv') into stats_file_ch

  script:
  """
  echo \
   "Sample\tNumber of reads mapped\tTotal bp mapped\tAverage depth\tCoverage" \
   > Mapping_stats.tsv
  """
}

process stats {
  publishDir "$params.out_dir/${bam.baseName}", pattern: '*.txt', mode: 'copy'
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
    samtools index $bam
    samtools idxstats $bam > ${bam.baseName}_idxstats.txt
    samtools depth $bam > ${bam.baseName}_depth.txt
    re_map=\$(awk '{if(NR==1) print \$3}' ${bam.baseName}_idxstats.txt)
    tot_bp=\$(awk '{sum+=\$3}END{print sum}' ${bam.baseName}_depth.txt)
    depth=\$(awk '{sum+=\$3;cnt++}END{print sum/cnt}' ${bam.baseName}_depth.txt| sed 's/,/./')
    ref_len=\$(awk '{if(NR==1) print \$2}' ${bam.baseName}_idxstats.txt)
    cov=\$(echo "scale=2 ; \$tot_bp / \$ref_len" | bc)
    echo "${bam.baseName}\t\$re_map\t\$tot_bp\t\$depth\t\$cov" >> $stats_file
  """
}

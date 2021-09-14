#!/usr/bin/env nextflow
log.info """\
         Fancy nextflow mapping script
         ===================================
         read_dir        : ${params.read_dir}
         ref             : ${params.ref_dir}
         """
         .stripIndent()

Channel
  .fromPath("${params.ref_dir}*.fasta*")
  .set{ref_index_ch}

Channel
  .fromFilePairs("${params.read_dir}*/*_R{1,2}.fastq.gz")
  .set{readPE_ch}


process index_ref {

  input:
  file(ref) from ref_index_ch

  output:
  tuple file(ref), file('*') into ref_mapPE_ch

  script:
  """
    bwa-mem2 index $ref
    samtools faidx $ref
  """
}

ref_and_reads = ref_mapPE_ch.combine(readPE_ch)

process map_PE {
  input:
  tuple file(ref), file('*'), val(sample_id), file(reads) from ref_and_reads

  output:
  file('*.sam') into mapped_PE

  script:
  """
    bwa-mem2 mem $ref $reads > ${sample_id}_PE.sam
  """

}

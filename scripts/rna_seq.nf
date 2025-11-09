nextflow.enable.dsl=2

params.reads = "raw_reads/*_{1,2}.fastq.gz"    
params.index = "reference/grch38/genome" 
params.gtf = "reference/Homo_sapiens.GRCh38.115.gtf"               
params.outdir = "results"
params.threads = 8
params.featurecounts_threads = 4
params.paired = true

process hisat2_align {
  tag { sample_id }
  publishDir "${params.outdir}/alignment", mode: 'copy'

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}.sorted.bam")

\
  script:
  """
  R1=${reads[0]}
  R2=${reads[1]}
  hisat2 -x ${params.index} -1 ${R1} -2 ${R2} --rg-id ${sample_id} --rg "SM:${sample_id}" --threads ${task.cpus} \
    | samtools view -bS - \
    | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
  samtools index ${sample_id}.sorted.bam
  """
}

process featurecounts_run {
  tag "featureCounts"
  publishDir "${params.outdir}/counts", mode: 'copy'

  input:
  set val(sample_id), path(bam) from hisat2_align.out

  output:
  path "featureCounts.txt"



  script:
  """
  featureCounts -T ${task.cpus} -a ${params.gtf} -o featureCounts.txt ${bam}

  """
}

workflow {
  // Channel: pair FASTQ files by sample
  reads_ch = Channel.fromFilePairs(params.reads, flat: false)

  // Map to (sample_id, [R1,R2]) tuples
  paired = reads_ch.map { sample, pair -> tuple(sample, pair) }

  // Run HISAT2 alignment producing sorted BAMs
  paired | hisat2_align

  // Collect BAMs and run featureCounts once
  hisat2_align.out.collect().set { bam_list }
  bam_list | featurecounts_run
}

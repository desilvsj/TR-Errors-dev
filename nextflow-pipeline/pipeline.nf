#!/home/desil/.local/bin/nextflow

//rm -rf .nextflow* work/ results/
//nextflow run pipeline.nf -profile docker --r1 inputs/example_trimmed_R1.fastq --r2 inputs/example_trimmed_R2.fastq --n  10000 -with-docker tr-errors:dev

/*
 * Running our consensus
 */

process consensusBuilder {
  publishDir 'output/Step-1', mode: 'copy', overwrite: true
  container 'tr-errors:dev'

  input:
    path r1_input
    path r2_input
    val  max_reads

  output:
    path 'output.fastq.gz'
    path 'metadata.txt.gz'

  script:
  """
  python3 /mnt/c/Janesh/GitHub/TR-Errors-Pipeline-V2/nextflow-pipeline/src/consensus_builder.py $r1_input $r2_input --max-reads $max_reads --fastq-out output.fastq.gz --meta-out metadata.txt.gz
  """
}

// process kallistoRun {
//     publishDir 'output/Step-2', mode: 'copy'

//     input
// }


/*
 * Pipeline parameters
 */
// params.r1 / params.r2 for quick testing; later you can swap to a samplesheet
params.r1 = null
params.r2 = null
params.n  = 10000  // max_reads (static; could pass as val if per-sample later)

workflow {
  if( !params.r1 || !params.r2 ){
    log.error 'Provide --r1 and --r2'; System.exit(1)
  }

  // One channel; one item per sample
  r1_ch = Channel.fromPath(params.r1)
  r2_ch = Channel.fromPath(params.r2)
  n_max_ch = Channel.of(params.n)

  // Step 1
  consensusBuilder(r1_ch, r2_ch, n_max_ch)
}

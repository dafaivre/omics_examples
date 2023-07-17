/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"

log.info """\
	PIPELINE INPUT PARAMS
	=====================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

// read_pairs_ch = Channel.fromFilePairs(params.reads)
Channel
	.fromFilePairs( params.reads, checkIfExists: true )
	.set { read_pairs_ch } 
	
read_pairs_ch.view()
/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"
params.transcript = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "$baseDir/output"

// println "reads: $params.reads"

log.info """
	PIPELINE INPUT PARAMS
	===============
	transcriptome		${params.transcript}
	reads			${params.reads}
	multiqc folder		${params.multiqc}
	outdir			${params.outdir}
	"""
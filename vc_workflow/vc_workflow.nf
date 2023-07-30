/*
 * Define the default parameters
 */ 
params.outdir = 'results'
params.genome = "$baseDir/data/ref_genome/ecoli_rel606.fasta"
params.reads = "$baseDir/data/trimmed_fastq/*_{1,2}.trim.sub.fastq"

log.info """\
PIPELINE INPUT PARAMS
=====================
genome       : ${params.genome}
reads        : ${params.reads}
outdir       : ${params.outdir}
"""
.stripIndent()

/*
 * Align reads to reference genome & create BAM file with FASTQC
 */
process FASTQC {
    tag{"FASTQC ${reads}"}

    publishDir("${params.outdir}/fastqc_trim", mode: 'copy')

    input:
    tuple val( sample_id ), path( reads )

    output:
    path( "*_fastqc*" ), emit: fastqc_out

    script:
    """
    fastqc ${reads}
    """
}

/*
 * Index the reference genome with BWA
 */
process BWA_INDEX {
	tag{"BWA_INDEX ${genome}"}

	publishDir("${params.outdir}/bwa_index", mode: 'copy')

	input:
	path genome

	output:
	tuple path( genome ), path( "*" ), emit: bwa_index

	script:
	"""
	bwa index ${genome}
	"""
}

/*
 * Align reads to reference genome & create BAM file
 */
process BWA_ALIGN {
    tag{"BWA_ALIGN ${sample_id}"}

    publishDir("${params.outdir}/bwa_align", mode: 'copy')

    input:
    tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    tuple val( sample_id ), path( "${sample_id}.aligned.bam" ), emit: aligned_bam

    script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    bwa mem \$INDEX ${reads} > ${sample_id}.aligned.sam
    samtools view -S -b ${sample_id}.aligned.sam > ${sample_id}.aligned.bam
    """
}

/*
 * Convert the format of the alignment to sorted BAM.
 */
process SAMTOOLS_SORT {
	tag{"SAMTOOLS_SORT ${sample_id}"}

	publishDir("${params.outdir}/bam_align", mode: 'copy')

	input:
	tuple val( sample_id ), path( bam )

	output:
	tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam" ), emit: sorted_bam

	script:
	"""
	samtools sort -o "${sample_id}.aligned.sorted.bam" ${bam}
	"""
}


/*
 * Index the BAM file for visualization purpose
 */
process SAMTOOLS_INDEX {
	tag{"SAMTOOLS_INDEX ${sample_id}"}
	
	publishDir("${params.outdir}/samtools_index", mode: 'copy')
	
	input:
	tuple val( sample_id ), path( sorted_bam )

	output:
	tuple path( sorted_bam ), path( "*" ), emit: samtools_index

	script:
	"""
	samtools index ${sorted_bam}
	"""
}

/*
 * Calculate the read coverage of positions in the genome.
 */
process BCFTOOLS_MPILEUP {
	tag{"BCFTOOLS_MPILEUP ${sample_id}"}
	
	publishDir("${params.outdir}/bcf", mode: 'copy')

	input:
	tuple val( sample_id ), path( sorted_bam ), path( genome )

	output:
	tuple val( sample_id ), path( "${sample_id}_raw.bcf" ), emit: raw_bcf

	script:
	"""
	bcftools mpileup -O b -o ${sample_id}_raw.bcf -f ${genome} ${sorted_bam}
	"""
}

/*
 * Detect the single nucleotide variants (SNVs).
 */
process BCFTOOLS_CALL {
	tag{"BCFTOOLS_CALL ${sample_id}"}
	
	publishDir("${params.outdir}/vcf", mode: 'copy')

	input:
	tuple val( sample_id ), path( raw_bcf )

	output:
	tuple val( sample_id ), path( "${sample_id}_variants.vcf" ), emit: variants

	script:
	"""
	bcftools call --ploidy 1 -m -v -o ${sample_id}_variants.vcf $raw_bcf
	"""
}

/*
 * Filter and report the SNVs in VCF (variant calling format).
 */
process VCFUTILS {
	tag{"VCFUTILS ${sample_id}"}
	
	publishDir("${params.outdir}/vcf", mode: 'copy')

	input:
	tuple val( sample_id ), path( variants )

	output:
	tuple val( sample_id ), path( "${sample_id}_final_variants.vcf" ), emit: final_variants

	script:
	"""
	vcfutils.pl varFilter $variants > ${sample_id}_final_variants.vcf
	"""
}

/* 
 * main workflow
 */
workflow {
	ref_ch = Channel.fromPath( params.genome, checkIfExists: true )
	reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )

	// Data preparation
	FASTQC( reads_ch )
	BWA_INDEX( ref_ch )
	BWA_ALIGN( BWA_INDEX.out.bwa_index.combine(reads_ch) )
	SAMTOOLS_SORT( BWA_ALIGN.out.aligned_bam )
	SAMTOOLS_INDEX( SAMTOOLS_SORT.out.sorted_bam )
	BCFTOOLS_MPILEUP( SAMTOOLS_SORT.out.sorted_bam.combine(ref_ch) )
	BCFTOOLS_CALL( BCFTOOLS_MPILEUP.out.raw_bcf )
	VCFUTILS( BCFTOOLS_CALL.out.variants )
}

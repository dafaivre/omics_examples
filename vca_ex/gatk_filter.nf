/*
 * Define the default parameters
 */ 
params.genome     = "$baseDir/data/genome.fa"
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.denylist   = "$baseDir/data/denylist.bed" 
params.reads      = "$baseDir/data/reads/rep2_{1,2}.fq.gz"
params.results    = "results"

log.info """\
PIPELINE INPUT PARAMS
=====================
genome   : $params.genome
reads    : $params.reads
variants : $params.variants
denylist : $params.denylist
results  : $params.results
"""

/*
 * Create a FASTA genome index (.fai) with samtools for GATK
 */

process PREP_GENOME_SAMTOOLS { 
  tag "$genome.baseName"
 
  input: 
    path genome
 
  output: 
    path "${genome}.fai"
  
  script:
  """
  samtools faidx ${genome}
  """
}

/*
 * Create a FASTA genome sequence dictionary with Picard for GATK
 */

process PREP_GENOME_PICARD {
  tag "$genome.baseName"
  label 'mem_xlarge'

  input:
    path genome
  output:
    path "${genome.baseName}.dict"

  script:
  """
  gatk CreateSequenceDictionary -R $genome -O ${genome.baseName}.dict
  """
}


/*
 * Create STAR genome index file.
 */

process PREP_STAR_GENOME_INDEX {
  tag "$genome.baseName"

  input:
    path genome
  output:
    path "genome_dir"

  script:
  """
  mkdir genome_dir

  STAR --runMode genomeGenerate \
       --genomeDir genome_dir \
       --genomeFastaFiles ${genome} \
       --runThreadN ${task.cpus}
  """
}


/*
 * Create a file containing the filtered and recoded set of variants
 */

process PREP_VCF_FILE {
  tag "$variantsFile.baseName"

  input: 
    path variantsFile
    path denylisted

  output:
    tuple \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz"), \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")
  
  script:  
  """
  vcftools --gzvcf $variantsFile -c \
           --exclude-bed ${denylisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz

  tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
  """
}

/*
 * Align reads to the genome with STAR
 */

process MAPPING_STAR {
  tag "$replicateId"

  input: 
    path genome
    path genomeDir
    tuple val(replicateId), path(reads) 

  output: 
    tuple \
      val(replicateId), \
      path('Aligned.sortedByCoord.uniq.bam'), \
      path('Aligned.sortedByCoord.uniq.bam.bai')

  script:
  """
  # ngs-nf-dev Align reads to genome
  STAR --genomeDir $genomeDir \
       --readFilesIn $reads \
       --runThreadN $task.cpus \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999

  # Run 2-pass mapping (improve alignmets using table of splice junctions and create a new index)  
  STAR --genomeDir $genomeDir \
       --readFilesIn $reads \
       --runThreadN $task.cpus \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --sjdbFileChrStartEnd SJ.out.tab \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:$replicateId LB:library PL:illumina PU:machine SM:GM12878

  # Select only unique alignments, no multimaps
  (samtools view -H Aligned.sortedByCoord.out.bam; samtools view Aligned.sortedByCoord.out.bam| grep -w 'NH:i:1') \
  |samtools view -Sb - > Aligned.sortedByCoord.uniq.bam
  
  # Index the BAM file
  samtools index Aligned.sortedByCoord.uniq.bam
  """
}

/*
 * Split reads that contain Ns in their CIGAR string.
 * Creates k+1 new reads (where k is the number of N cigar elements) 
 * that correspond to the segments of the original read beside/between 
 * the splicing events represented by the Ns in the original CIGAR.
 */

process GATK_SPLITNCIGAR {
  tag "$replicateId"
  label 'mem_large'
  
  input: 
    path genome
    path index
    path genome_dict
    tuple val(replicateId), path(bam), path(index)

  output:
    tuple val(replicateId), path('split.bam'), path('split.bai')
  
  script:
  """
  # SplitNCigarReads and reassign mapping qualities
  gatk SplitNCigarReads \
            -R $genome \
            -I $bam \
            --refactor-cigar-string \
            -O split.bam
  """
}

/* 
 * main workflow
 */
workflow {
      reads_ch = Channel.fromFilePairs(params.reads)

      // Data preparation
      PREP_GENOME_SAMTOOLS(params.genome)
      PREP_GENOME_PICARD(params.genome)
      PREP_STAR_GENOME_INDEX(params.genome)
      PREP_VCF_FILE(params.variants, params.denylist)
	  
      // STAR Mapping
      MAPPING_STAR( 
            params.genome, 
            PREP_STAR_GENOME_INDEX.out, 
            reads_ch )

      // GATK Prepare Mapped Reads
      GATK_SPLITNCIGAR(
            params.genome, 
            PREP_GENOME_SAMTOOLS.out, 
            PREP_GENOME_PICARD.out, 
            MAPPING_STAR.out )
}
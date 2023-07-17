/*
 * Define the default parameters
 */ 
params.genome     = "$baseDir/data/genome.fa"
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.denylist   = "$baseDir/data/denylist.bed" 
params.reads      = "$baseDir/data/reads/rep1_{1,2}.fq.gz"
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
 * main workflow
 */
workflow {
      reads_ch = Channel.fromFilePairs(params.reads)

      // Data preparation
      PREP_GENOME_SAMTOOLS(params.genome)
      PREP_GENOME_PICARD(params.genome)
      PREP_STAR_GENOME_INDEX(params.genome)
      PREP_VCF_FILE(params.variants, params.denylist)
}
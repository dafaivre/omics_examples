# omics_examples

## Simple Example
A simple Nextflow/NGS example is in the `simple_ex` folder.

The workflow includes `salmon`, `fastqc`, and `multiqc`.

The `rna-seq.nf` file is the full workflow. The other files were drafts of the steps. 

## Variant Calling Workflow Example 1
A variant calling workflow example is in the `vc_workflow` folder.

The workflow includes `fastqc`, `bwa`, `samtools`, `bcftools`, and `vcfutils.pl`.

## Variant Calling Analysis Example
An in-progress variant calling analysis example is in the `vca_ex` folder.

The workflow currently includes `samtools`, `gatk`, `STAR`, `vcftools`, and `tabix`.
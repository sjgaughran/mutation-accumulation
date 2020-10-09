# mutation-accumulation

## Contents
- [About](#about)
- [SNP Annotation](#snp-annotation)
- [Statistic calculations](#statistic-calculations)
  * [Notes on bootstrapping](#notes-on-bootstrapping)


## About

This mutation accumulation analysis pipeline provides a guide for quantifying differences in mutation accumulation between populations. It follows the statistical approaches proposed by [Simons *et al.* (2014)](https://www.nature.com/articles/ng.2896) (summing alleles) and [Do *et al.* (2015)](https://www.nature.com/articles/ng.3186) (*R<sub>XY</sub>* and related statistics). The Simons *et al.* (2014) method counts derived alleles and is useful for comparing closely related populations. The Do *et al.* (2015) method also takes a counting approach, but normalizes the count for shared derived alleles in a given population pair and is therefore appropriate for phylogeny-level comparisons. In either case, this pipeline assumes that alleles can be polarized by a 
that is outgroup to all other populations of interest. This assumption may be violated by long branches to the outgroup or if there has been gene flow between the outgroup and any of the target populations. 

## SNP Annotation

[SnpEff](https://pcingola.github.io/SnpEff/) is a versatile genomic variant annotation tool. In this pipeline, a SnpEff-annotated VCF is required. The VCF file should containing individuals from the target populations and the outgroup individual. Though not required, we recommend generating the VCF by calling only variants in annotated CDS regions, which will be faster and speed up some downstream analyses. VCFs should be filtered for quality and bi-allelic variants. Additional filtering--such as including only one isoform/transcript per gene--is recommended.

If working with non-model species, you can build an annotation database for the species by following the [SnpEff documentation](https://pcingola.github.io/SnpEff/se_buildingdb/).

## Statistic calculations

The statistical analysis and bootstrapping are done through three major steps:
1) split the SnpEff-annotated VCF into 1000 VCFs of equal size
2) run the Simons *et al.* (2014) method and the first step (*L<sub>XY</sub>* statistic)

*NB: The efficiency of the pipeline could be improved.*




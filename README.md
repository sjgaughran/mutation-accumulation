# mutation-accumulation

## Contents
- [About](#about)
- [SNP Annotation](#snp-annotation)
- [Statistic calculations](#statistic-calculations)
  * [Notes on bootstrapping](#notes-on-bootstrapping)


## About

This mutation accumulation analysis pipeline provides a guide for quantifying differences in mutation accumulation between populations. It follows the statistical approaches proposed by [Simons *et al.* (2014)](https://www.nature.com/articles/ng.2896) (summing alleles) and [Do *et al.* (2015)](https://www.nature.com/articles/ng.3186) (*R<sub>XY</sub>* and related statistics). The Simons *et al.* (2014) method counts derived alleles and is useful for comparing closely related populations. The Do *et al.* (2015) method also takes a counting approach, but normalizes the count for shared derived alleles in a given population pair and is therefore appropriate for phylogeny-level comparisons. In either case, this pipeline assumes that alleles can be polarized by a 
that is outgroup to all other populations of interest. This assumption may be violated by long branches to the outgroup or if there has been gene flow between the outgroup and any of the target populations. 

## Requirements

- [SnpEff](https://pcingola.github.io/SnpEff/)
- Python3

## SNP Annotation

[SnpEff](https://pcingola.github.io/SnpEff/) is a versatile genomic variant annotation tool. In this pipeline, a SnpEff-annotated VCF is required. The VCF file should containing individuals from the target populations and the outgroup individual. Though not required, we recommend generating the VCF by calling only variants in annotated CDS regions, which will be faster and speed up some downstream analyses. VCFs should be filtered for quality and bi-allelic variants. Additional filtering--such as including only one isoform/transcript per gene--is recommended.

If working with non-model species, you can build an annotation database for the species by following the [SnpEff documentation](https://pcingola.github.io/SnpEff/se_buildingdb/).

## Statistic calculations

The statistical analysis and bootstrapping are done through three separate steps:
1) **Split the SnpEff-annotated VCF into 1000 sub-VCFs of equal size.** Run *split_snpeff_vcf.py*, specifying the SnpEff-annotated VCF and the total number of lines (e.g. `wc -l snpeff.vcf`) in the VCF.

```
python split_snpeff_vcf.py snpeff.vcf 1000000
```

2) **Calculate the Simons *et al.* (2014) method and the first step (*L<sub>XY</sub>* statistic) of Do *et al.* (2015) on each individual sub-VCF.** Run *read_snpeff_bs.py* on every sub-VCF (all of which are called *snpeff_bs_1.vcf*, *snpeff_bs_2.vcf*, ..., *snpeff_bs_1000.vcf*). This will output results files for each sub-VCF, which are then processed in the next step. To run this script, you must specify a text file that contains a one population name per line, with sample names listed after a colon (:) and separated by commas. For example:

```
POP1: indiv1, indiv3, indiv5
POP2: indiv2, indiv4, indiv6
POP3: indiv7, indiv8
OUTGROUP: indiv9
```

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;You must also specify the name of the outgroup population as found in the population text file. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The script is then run on each sub-VCF as:

```
python read_snpeff_bs.py snpeff_bs_{number}.py pops.txt outgroup_pop_name
```

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Because the sub-VCFs are numbered, a simple shell script can run the above command on each.

3) **Bootstrap by sampling with replacement from the sub-VCFs 1000 times, and output the average of each statistic with an upper and a lower bound to the 95% CI.**



*NB: The efficiency of the pipeline could be improved.*




# mutation-accumulation

## Contents
- [About](#about)
- [SNP Annotation](#snp-annotation)
- [Statistic calculations](#statistic-calculations)
  * [Notes on bootstrapping](#notes-on-bootstrapping)


## About

This mutation accumulation analysis pipeline provides a guide for quantifying differences in mutation accumulation between populations. It follows the statistical approaches proposed by [Simons *et al.* (2014)](https://www.nature.com/articles/ng.2896) (summing alleles) and [Do *et al.* (2015)](https://www.nature.com/articles/ng.3186) (Rxy and related statistics). The Simons *et al.* (2014) method counts derived alleles and is useful for comparing closely related populations. The Do *et al.* (2015) method also takes a counting approach, but normalizes the count for shared derived alleles in a given population pair and is therefore appropriate for phylogeny-level comparisons. In either case, this pipeline assumes that alleles can be polarized by a species that is outgroup to all other species of interest. This assumption may be violated by long branches to the outgroup or if there has been gene flow between the outgroup and any of the target species. 

## SNP Annotation



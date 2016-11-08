.. _output:

CSV Output
==========

This document describes the columns in the CSV output


sex_check
---------

Sex check performs a comparison between the sex reported in the ped
file and that inferred from the genotypes on the non-PAR regions of
the X chromosome.

1 row per sample with columns of:

+ sample_id: sample from ped.
+ error: boolean indicating wether there is a mismatch between X genotypes and ped sex.
+ het_count: number of heterozygote calls
+ hom_alt_count: number of homozygous-alternate calls
+ hom_ref_count: number of homozygous-reference calls
+ het_ratio: ratio of `het_count` / `hom_alt_count`. Low for males, high for females
+ ped_sex: sex from .ped file 
+ predicted_sex: sex predicted from rate of hets on chrX.


het_check
---------

Het check does general QC including rate of het calls, allele-balance at het calls,
mean and median depth, and a PCA projection onto thousand genomes.

1 row per sample with columns of:

+ sample_id: sample from ped.
+ sampled_sites: number of sites sampled (sufficient call-rate across samples and depth in this sample)
+ mean/median_depth: mean/median depths for the sites tested.
+ depth_outlier: boolean indicating that this sample's depth is considered an outlier relative to the other samples.
+ het_count: number of heterozygote calls in sampled sites. 
+ het_ratio: proportion of sites that were heterozygous.
+ ratio_outlier: boolean indicating that the het_ratio was outside what is normally seen.
+ idr_baf: inter-decile range (90th percentile - 10th percentile) of b-allele frequency. We make a distribution of all sites of
  alts / (ref + alts) and then report the difference between the 90th and the 10th percentile. Large values indicated
  likely sample contamination.
+ p10/p90: the numbers used to calculate idr_baf.

And the PCA columns:

+ PC1/PC2/PC3/PC4: the first 4 values after this sample was projected onto the thousand genomes principle components.
+ ancestry-prediction: one of `AFR AMR EAS EUR SAS UNKNOWN` where it is unknown if `ancestry-prob` < 0.65 for the
  highest population 
+ ancestry-prob: the highest probability from the SVM for any ancestry (between 0 and 1).


ped_check
---------

Ped check compares the relatedness of 2 samples as reported in a .ped file to the
relatedness inferred from the genotypes and ~25K sites in the genome.

This contains 1 row per sample-pair: (n_samples * n_samples) / 2 rows.

+ sample_a/sample_b: the samples indicating the pair in question.
+ n: the number of sites that was used to predict the relatedness.
+ rel: the relatedness calculated from the genotypes.
+ pedigree_relatedness: the relatedness reported in the ped file.
+ rel_difference: difference between the preceding 2 colummns.
+ ibs0: the number of sites at which the 2 samples shared no alleles (should approach 0 for parent-child pairs).
+ ibs2: the number of sites and which the 2 samples where both hom-ref, both het, or both hom-alt.
+ shared_hets: the number of sites at which both samples were hets.
+ hets_a/b: the number of sites at which sample_a/b was het.
+ pedigree_parents: boolean indicating that this pair is a parent-child pair according to the ped file.
+ predicted_parents: boolean indicating that this pair is expected to be a parent-child pair according to the ibs0 (< 0.012) calculated from the genotypes.
+ parent_error: boolean indicating that the preceding 2 columns don't match
+ sample_duplication_error: boolean indicating that rel > 0.75 and ibs0 < 0.012



ancestry
--------

The ancestry check is included in the `het_check` output. Here, we describe its implementation in more detail.
The ancestries of the thousand-genomes samples, along with their genotypes at selected sites are distributed with peddy.
The size of the genotypes for the 2504 samples is about 10MB. When running peddy, we sample those same sites on
the VCF sent in by the user; however, some sites might be excluded due to poor quality or lack of coverage in then
cohort. So, we subset our thousand-genomes sites to those that are also sample in the cohort. (Usually, the intersection
is very high). This gives a matrix of genotypes for the 1KG samples that we know are also represented in the current
query cohort. **We want to train classifier to predict ancestry from genotypes**. To do this efficiently, we first do a
dimensionality reduction using randomized PCA to change the matrix of 2504 samples by ~23,000 sites into a matrix of
2504 samples by 4 principal components. We then train an SVM on that reduced matrix given the known ancestries of the
1KG samples. Now we have a classifier.

To classify the query cohort, we do a dimensionality reduction by projecting the genotype matrix onto the principal
components of the 1KG cohort. We take the resulting, reduced matrix and use the SVM to predict the ancestry of each
sample. The SVM reports a confidence in the prediction. We assign the most likely ancestry if the prediction is greater
than 0.65.

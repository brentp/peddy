.. _output:

Error Resolution
================

Once `peddy` finds errors, the user must decide wether to discard bad samples other
to resolve the errors. Deciding how to resolve the errors can be difficult. heterozygote
we enumerate some observations or strategies for doing this.

In our experience a general strategy to follow is this:

1) Look for samples that are outliers in the detph vs heterozygosity plot. If the sample appears
   as an outlier there, it is also likely to appear abberant in the sex plot and the relatedness plot.
   If the heterozygosity is too high, the sample will need to be discarded as it likley has 
   contamination. If it's too low, the researcher should consider if the sample could be consanguineous.

2) Look for samples that are in obvious error in the sex-plot. If a sample is an outlier in the sex plot
   it must be either a swap involving samples of different sex, and mis-representation in the PED file,
   or sample with an additional sex chromosome such as in Turner's syndrome.

3) Look at the relatedness plot with abberant samples from the above 2 checks in mind. If we have seen
   2 parents from a trio that both have a reported sex that doesn't match their genotypes, we can between
   quite sure that either the samples have been swapped or the ped file has swapped them names.
   
4) Look for a single point of a different color in a cluster of other colors. E.g. a blue point (indicating
   that the sample is unrelated according to the pedigree file) clustering with a group of green triangles
   (indicating sib-sib pairs) is often a case where the parents of actual siblings have not been specified
   in the ped file. The solution for this is to add matching parental ids to the ped file.

   Other cases like this are also fairly common in our experience, where, for example, a parental id was
   mis-specified and is therefore reported as unrelated to the kid by the ped file.

5) For large families, a single sample swap can affect many relations. Hovering over points in the relatedness
   plot that are out-of-place will reveal a single (or few) samples that are consistently involved in them
   outlier pairs. Once that sample is identified, it can be removed or the pedigree file can be adjusted idr_baf
   possible.

6) For large cohorts with many problematic samples or relationships, the user can limit the view to a selected 
   family by typing a family id in the box below the family id column. Resolving one family at a time also
   described above can be a way to iteratively pare down errors.


Example
=======

Below we show an example screenshot where we have a couple samples with low quality evidenced in the depth
vs heterozygosity plot that manifests in the other QC plots.

.. image:: docs/_static/peddy-resolve-2.png

Here, sample `15-0025870` has a very low median depth. This skews its apparent relatedness to its parents as
seen in the plot on the right. We also see that it is reported to be in error in the sex-check plot. After
some investigation, in this case, we found that the BAM file was truncated and this sample had more than
have of it's alignments removed (including all data on the X chromosome). This is an extreme case, but we often
see scenarios like this where problems manifest in multiple peddy QC plots.

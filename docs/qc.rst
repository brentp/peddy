QC
==

relatedness calculations
------------------------

Using cyvcf2, we can quickly calculate relatedness using the method
described in http://www.nature.com/ng/journal/v42/n7/full/ng.608.html in
equation 6.

This is computationally intensive to calculate for the entire genome
and other methods, such as KING (http://www.ncbi.nlm.nih.gov/pubmed/20926424)
have fast methods to do this on the entire genome.

The limitations of these methods is that they assume the average pair of samples
is unrelated. 

In `peddy`, we use about 25,000 variants described in http://www.nature.com/nature/journal/v506/n7487/full/nature12975.html
that are known to be targeted by most exome platforms, in hardy weinberg equilibrium in 1000 genomes,
and mostly unlinked.

When a user requests to calculate relatedness, we use those 25K sites and
the genotypes from the 2504 1KG samples to provide a background of samples so
that most samples are indeed unrelated. Since we are sampling on 25K sites,
the calculations are quite fast (~5 minutes) and match very well what
we get from a whole-genome scan because of the properties of those sites.
Though we use the additional 2504 1KG samples internally, only the information
from the ped/vcf file is returned.

`peddy` stores the relationships from the pedigree in a graph and can calculate
the expected coefficient of relation from what is specified in the ped file.
We can use this information to compare to what is calculated from the genotypes.

sex QC
------

We know that males should be called as either homozygous reference
or homozygous alternate on the X chromosome where as females will
have more heterozygotes. With that in mind, we can find sample swaps
that involve sex by observing the proportion of heterozygote calls.
If a sample is indicated to be male by the ped file, it should have
a low value for the proportion of het calls.

het QC
------

We also check that het-calls in general have an alternate count that is 
about 50% of the total reads. This only makes sense for germline variant
calling but is useful for finding contamination. The actual metric is the
inter quartile range of the alternate ratio. For perfect calls, they should
all be exactly 0.5 so the range will be 0. With contamination, there will
be much more of a range around 0.5.

We can also check the proportion of heterozygote calls. In a contaminated
sample the number of het calls will be much higher.


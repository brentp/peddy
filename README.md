Fast Pedigree::VCF QC
---------------------

<img src="https://cloud.githubusercontent.com/assets/1739/21734388/5a927926-d420-11e6-8ff9-534046df3305.png" style="width:287px;height:168px">

peddy compares familial-relationships and sexes as reported in a [PED/FAM file](https://www.cog-genomics.org/plink2/formats#fam)
with those inferred from a VCF.

It samples the VCF at about 25000 sites (plus chrX) to accurately estimate **relatedness**, **IBS0**, **heterozygosity**, **sex** and **ancestry**. It uses 2504 thousand genome samples as backgrounds to calibrate the relatedness calculation and to make ancestry predictions.

It does this very quickly by sampling, by using C for computationally intensive parts, and by parallelization.

If you use peddy, please cite [Pedersen and Quinlan, Who’s Who? Detecting and Resolving Sample Anomalies in Human DNA
Sequencing Studies with Peddy, The American Journal of Human Genetics (2017),
http://dx.doi.org/10.1016/j.ajhg.2017.01.017](http://www.cell.com/action/showFullTextImages?pii=S0002-9297(17)30017-4)


[![Anaconda-Server Badge](https://anaconda.org/bioconda/peddy/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![PyPI version](https://badge.fury.io/py/peddy.svg)](http://badge.fury.io/py/peddy)
[![Documentation Status](https://readthedocs.org/projects/peddy/badge/?version=latest)](http://peddy.readthedocs.org/en/latest/?badge=latest)
<!--
[![Build Status](https://travis-ci.org/brentp/peddy.svg?branch=master)](https://travis-ci.org/brentp/peddy)
-->

Note that [somalier](https://github.com/brentp/somalier) is a more scalable, faster, replacement for peddy
that uses some of the same methods as peddy along with some new ones.

Quickstart
----------

See installation below.

Most users will only need to run as a command-line tool with a ped and VCF, e.g:

```
python -m peddy -p 4 --plot --prefix ceph-1463 data/ceph1463.peddy.vcf.gz data/ceph1463.ped
```

This will use 4 cpus to run various checks and create **ceph-1463.html** which
you can open in any browser to interactively explore your data.

It will also create create 4 csv files and 4 QC plots.
These will indicate:

+ discrepancies between ped-reported and genotype-inferred relations
+ discrepancies between ped-reported and genotype-inferred sex
+ samples with higher levels of HET calls, lower depth, or more variance in b-allele-frequency (ref / (ref + alt )) for het calls.
+ an ancestry prediction based on projection onto the thousand genomes principal components

Finally, it will create a new file ped files `ceph1463.peddy.ped` that also lists
the most useful columns from the `het-check` and `sex-check`. Users can **first
look at this extended ped file for an overview of likely problems**.

See [the docs](http://peddy.readthedocs.io/en/latest/) for a walk-through and thorough explanation of each plot.

hg38 or custom sites
--------------------

By default, peddy uses hg19/GRCh37. It can be forced to use sites for hg38 by passing `--sites hg38`.
To create custom sites, have a look at [the sites files included with peddy](https://github.com/brentp/peddy/blob/master/peddy/GRCH37.sites) along
with the corresponding .bin.gz which is just the raw binary alternate counts (gt_types) from thousand-genomes that have been written as uint8
and gzipped.


Speed
-----

Because of the sampling approach and parallelization, `peddy` is very fast.
With 4 CPUs, on the 17-member *CEPH1643* pedigree whole-genome VCF, peddy can run the het-check and PCA in ~ 8 seconds. The pedigree check comparing all vs.
all samples run in 3.6 seconds.
It finishes the full set of checks in about 20 seconds.

In comparison [KING](http://people.virginia.edu/~wc9c/KING/manual.html) runs
in 14 seconds (it is **extremely fast**); the time including the conversion
from VCF to binary ped is 85 seconds.

On larger datasets, with hundreds or thousands of samples, it can be beneficial to add as many
cores as possible; for smaller datasets with dozens of samples about 4 processors reduces the
computation time as much as 8 or more would.


Validation
----------

The results between peddy and KING are comparable, but peddy does better on
cohorts where most samples are related. See the figure below where the peddy
relatedness estimate is closer to the actual than KING which over-estimates relatedness.

![Peddy Vs KING](https://raw.githubusercontent.com/brentp/peddy/master/docs/_static/peddy-v-king.png "Comparison with KING")

Peddy uses the KING algorithm for calculating relatedness and so they
match quite well.
Peddy also runs PCA on the 2504 samples from 1000 genomes,
then fitting an SVM and predicting ancestry in addition to
calculating relatedness among all pairwise combinations of the 17 samples.

Warnings and Checks
-------------------

On creating a pedigree object (via Ped('some.ped'). `peddy` will print warnings to STDERR as appropriate like:

```
pedigree warning: '101811-101811' is dad but has female sex
pedigree warning: '101897-101897' is dad but has female sex
pedigree warning: '101896-101896' is mom of self
pedigree warning: '102110-102110' is mom but has male sex
pedigree warning: '102110-102110' is mom of self
pedigree warning: '101381-101381' is dad but has female sex
pedigree warning: '101393-101393' is mom but has male sex

unknown sample: 102498-102498 in family: K34175
unknown sample: 11509-11509 in family: K567331
unknown sample: 5180-5180 in family: K8565
```

Installation
------------

## Conda
Nearly all users should install using conda in the anaconda python distribution. This means
have your own version of python easily installed via:

```
INSTALL_PATH=~/anaconda
wget http://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
# or wget http://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
bash Miniconda2-latest* -fbp $INSTALL_PATH
PATH=$INSTALL_PATH/bin:$PATH

conda update -y conda
conda config --add channels bioconda

conda install -y peddy
```

This should install all dependencies so you can then run peddy with 4 processes as:

```
python -m peddy --plot -p 4 --prefix mystudy $VCF $PED
```

## Github

```
git clone https://github.com/brentp/peddy
cd peddy
pip install -r requirements.txt
pip install --editable .
```

run with 

```
peddy --plot -p 4 --prefix mystudy $VCF $PED
```

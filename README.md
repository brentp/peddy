Fast Pedigree::VCF QC
---------------------

peddy compares familial-relationships and sexes as reported in a [PED file](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped)
with those inferred from a VCF.

It samples the VCF at about 25000 sites (plus chrX) to accurately estimate **relatedness**, **IBS0**, **heterozygosity**, **sex** and **ancestry**. It uses 2504 thousand genome samples as backgrounds to calibrate the relatedness calculation and to make ancestry predictions.

It does this very quickly by sampling, by using C for computationally intensive parts, and by parallelization.


[![PyPI version](https://badge.fury.io/py/peddy.svg)](http://badge.fury.io/py/peddy)
[![Documentation Status](https://readthedocs.org/projects/peddy/badge/?version=latest)](http://peddy.readthedocs.org/en/latest/?badge=latest)
<!--
[![Build Status](https://travis-ci.org/brentp/peddy.svg?branch=master)](https://travis-ci.org/brentp/peddy)
-->


Quickstart
----------

See installation below.

Most users will only need to run as a command-line tool with a ped and VCF, e.g:

```
python -m peddy -p 12 --plot --prefix ceph-1463 data/ceph1463.vcf.gz data/ceph1463.ped
```

This will use 12 cpus to run various checks and create **ceph-1463.html** which
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

See [the docs](TODO) for a walk-through and thorough explanation of each plot.

Speed
-----

Because of the sampling approach and parallelization, `peddy` is very fast.
With 20 CPUs, on the 17-member *CEPH1643* pedigree whole-genome VCF, peddy can run
the het-check and PCA in < 20 seconds. The pedigree check including all vs.
all against the 2504 thousand genomes samples run in 85 seconds.
It finishes the full set of checks in under two minutes.

In comparison [KING](http://people.virginia.edu/~wc9c/KING/manual.html) runs
in 14 seconds (it is **extremely fast**); the time including the conversion
from VCF to binary ped is 85 seconds.


Validation
----------

The results between peddy and KING are comparable, but peddy does better on
cohorts where most samples are related. See the figure below where the peddy
relatedness estimate is closer to the actual than KING which over-estimates relatedness.

![Peddy Vs KING](https://raw.githubusercontent.com/brentp/peddy/master/docs/_static/peddy-v-king.png "Comparison with KING")

Note that the peddy analysis is well-calibrated as it runs with the thousand genomes samples
as background. It also includes running a PCA on the 2504 samples from 1000 genomes,
then fitting an SVM and predicting ancestry in addition to calculating relatedness
among all pairwise combinations of the 2504+17 samples.

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

Nearly all users should install using conda in the anaconda python distribution. This means
have your own version of python easily installed via:

```
INSTALL_PATH=~/anaconda
wget http://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
# or wget http://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
bash Miniconda2-latest\* -bp $INSTALL_PATH
PATH=$INSTALL_PATH/bin:$PATH

conda update -y conda
conda config --add channels bioconda

conda install -y peddy
```

This should install all dependencies so you can then run peddy with 8 processes as:

```
python -m peddy --plot -p 8 --prefix mystudy $VCF $PED
```

To get the development versions of peddy (and cyvcf2), you can follow the above steps and then do:

```
git clone https://github.com/brentp/cyvcf2
cd cyvcf2 && python setup.py install
cd ..
git clone https://github.com/brentp/peddy
cd peddy && python setup.py install
```

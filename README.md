tools for pedigree files
------------------------

[![PyPI version](https://badge.fury.io/py/peddy.svg)](http://badge.fury.io/py/peddy)
[![Build Status](https://travis-ci.org/brentp/peddy.svg?branch=master)](https://travis-ci.org/brentp/peddy)
[![Documentation Status](https://readthedocs.org/projects/peddy/badge/?version=latest)](http://peddy.readthedocs.org/en/latest/?badge=latest)


Quickstart
----------

See installation below.

Most users will only need to run as a command-line tool with a ped and VCF, e.g:

```
python -m peddy --plot --prefix ceph-1463 ceph1463.vcf.gz ceph1463.ped
```

This create **ceph-1463.html** which you can open in any browser to
interactively explore your data.

It will also create create 4 csv files and 4 QC plots.
These will indicate:

+ discrepancies between ped-reported and genotype-inferred relations
+ discrepancies between ped-reported and genotype-inferred sex
+ higher levels of HET calls or more variance in ref / (ref + alt read) for het calls.
+ an ancestry prediction based on to projection onto the thousand genomes principal components

Finally, it will create a new file ped files `ceph1463.peddy.ped` that also lists
the most useful columns from the `het-check` and `sex-check`. Users can **first
look at this extended ped file for an overview of likely problems** and then refer
to the plots and .csv files for more detailed information.

Overview
--------


**NOTE** this module used to be named to "pedagree".

`peddy` is a python library for querying, QC'ing, and manipulating pedigree files.

It currently makes it simple to extract things like:

 + parent-child pairs
 + trios
 + sibs
 + stats on number of males/females/trios/affecteds/unaffecteds
 + families.
 + families with at least N members
 + families with at least N children
 + [not yet] families with at least N generations
 + coefficient of relatedness given relation defined in the pedigree.

Also, given a pedigree file and a VCF file peddy provides tools to:

 + find likely sample mixups (or PED errors)
   - sex mixups on X-Chrom
   - family mixups by inferring relatedness with VCF

 + find mendelian errors


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

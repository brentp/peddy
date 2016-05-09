tools for pedigree files
------------------------

[![PyPI version](https://badge.fury.io/py/peddy.svg)](http://badge.fury.io/py/peddy)
[![Build Status](https://travis-ci.org/brentp/peddy.svg?branch=master)](https://travis-ci.org/brentp/peddy)
[![Documentation Status](https://readthedocs.org/projects/peddy/badge/?version=latest)](http://peddy.readthedocs.org/en/latest/?badge=latest)


Quickstart
----------

Most users will only need to run as a command-line tool with a ped and VCF, e.g:

```
python -m peddy --plot --prefix ceph-1463 ceph1463.vcf.gz ceph1463.ped
```

This create **ceph-1463.html** which you can open in any browser to
interactively explore your data.

It will also create create 3 csv files and 3 QC plots.
These will indicate:

+ discrepancies between ped-reported and genotype-inferred relations
+ discrepancies between ped-reported and genotype-inferred sex
+ higher levels of HET calls or more variance in ref / (ref + alt read) for het calls.

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


Usage
-----

```Python
>>> from peddy import Ped, SEX, PHENOTYPE

>>> p = Ped('my.ped')
# not yet.
#>>> p.dot() # draw the pedigree with graphviz

# not yet
# find any obvious issues (3 parents, mom as male, etc).
>>> p.validate()

# number of affecteds, un, males, females, etc. (contingency table?)
>>> p.summary()

# iterable
>>> p.samples()

>>> p.samples(phenotype=PHENOTYPE.AFFECTED, sex=SEX.MALE)

# sample object
>>> s = next(p.samples())

>>> s.phenotype

>>> s.sex

>>> s.mom

>>> s.dad

>>> s.siblings

>>> s.kids
```

Quality Control
---------------

If cyvcf2 is installed, then, given a ped-file and a VCF, we can look for cases where the relationships
defined in the ped file do not match the relationships derived from the genotypes in the VCF.

```Python
>>> from peddy import Ped
>>> p = Ped('cohort.ped')
>>> df = p.ped_check('cohort.vcf.gz')
>>> df[df.error] # show pairs of samples where the inferred differs from the reported.

```

[![relplot](http://peddy.readthedocs.org/en/latest/_images/ped-check.png)](http://github.com/brentp/cyvcf2/)

We don't see any obvious errors in this pedigree. An obvious error would be when a red colored dot clusters with blue dots. 
The *outlined dots* have a very low IBS0 rate, indicating that they are likely parent-child pairs.

By looking for the frequency of heterozygotes in the non-PAR regions of
the X chromosome, we can determine sex from a VCF:

```Python
>>> from peddy import Ped
>>> p = Ped('cohort.ped')
>>> p.sex_check('cohort.vcf.gz', plot=True)
... List of all samples with number of HETs, HOMREF, HOMALT on X
```
This will also create an image like this one where we can
see a clear sample mixup.

[![sex_plot](https://raw.githubusercontent.com/brentp/peddy/master/images/sex_check.png)](http://github.com/brentp/cyvcf2/)


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

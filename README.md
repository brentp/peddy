tools for pedigree files
------------------------

[![PyPI version](https://badge.fury.io/py/peddy.svg)](http://badge.fury.io/py/peddy)
[![Build Status](https://travis-ci.org/brentp/peddy.svg?branch=master)](https://travis-ci.org/brentp/peddy)

This is currently under development. It already has a nice API for dealing with pedigree files.

It currently makes it simple to extract things like:

 + parent-child pairs
 + trios
 + sibs
 + stats on number of males/females/trios/affecteds/unaffecteds
 + families.
 + families with at least N members
 + families with at least N children
 + [not yet] families with at least N generations


We are working on adding support for the following:

Given a pedigree file and a VCF file:

 + find likely sample mixups (or PED errors)
   - gender mixups on X-Chrom
   - family mixups by inferring relatedness with VCF

 + mendelian errors


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
>>> p.validate('cohort.vcf.gz')
... LIST of QUESTIONABLE SAMPLES
```

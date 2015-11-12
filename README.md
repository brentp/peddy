tools for pedigree files
------------------------

[![PyPI version](https://badge.fury.io/py/pedagree.svg)](http://badge.fury.io/py/pedagree)
[![Build Status](https://travis-ci.org/brentp/pedagree.svg?branch=master)](https://travis-ci.org/brentp/pedagree)

**NOTE** this module used to be named to "peddy".

`pedagree` is a python library for querying, QC'ing, and manipulating pedigree files.

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

Also, given a pedigree file and a VCF file pedagree provides tools to:

 + find likely sample mixups (or PED errors)
   - sex mixups on X-Chrom
   - family mixups by inferring relatedness with VCF

 + find mendelian errors


Usage
-----

```Python
>>> from pedagree import Ped, SEX, PHENOTYPE

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
>>> from pedagree import Ped
>>> p = Ped('cohort.ped')
>>> df = p.ped_check('cohort.vcf.gz')
>>> df[df.error] # show pairs of samples where the inferred differs from the reported.

```

[![relplot](https://raw.githubusercontent.com/brentp/pedagree/master/images/t.png)](http://github.com/brentp/cyvcf2/)


By looking for the frequency of heterozygotes in the non-PAR regions of
the X chromosome, we can determine sex from a VCF:

```Python
>>> from pedagree import Ped
>>> p = Ped('cohort.ped')
>>> p.sex_check('cohort.vcf.gz', plot=True)
... List of all samples with number of HETs, HOMREF, HOMALT on X
```
This will also create an image like this one where we can
see a clear sample mixup.

[![sex_plot](https://raw.githubusercontent.com/brentp/pedagree/master/images/sex_check.png)](http://github.com/brentp/cyvcf2/)


On creating a pedigree object (via Ped('some.ped'). `pedagree` will print warnings to STDERR as appropriate like:

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

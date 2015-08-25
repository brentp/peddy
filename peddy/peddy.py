# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
from collections import OrderedDict

if sys.version_info[0] == 3:
    basestring = str


class UNKNOWN(object):
    def __str__(self):
        return "-9"

class PHENOTYPE(object):
    AFFECTED = True
    UNAFFECTED = False
    UNKNOWN = UNKNOWN()

    @classmethod
    def lookup(cls, i):
        return {'2': cls.AFFECTED, '1': cls.UNAFFECTED}.get(i, cls.UNKNOWN)

    @classmethod
    def rlookup(cls, i):
        return {cls.AFFECTED: 'affected', cls.UNAFFECTED: 'unaffected'}.get(i, cls.UNKNOWN)


class SEX(object):
    MALE = 'male'
    FEMALE = 'female'
    UNKNOWN = UNKNOWN()

    @classmethod
    def lookup(cls, i):
        return {'1': cls.MALE, '2': cls.FEMALE}.get(i, cls.UNKNOWN)

    @classmethod
    def rlookup(cls, i):
        return {cls.MALE: 'male', cls.FEMALE: 'female'}.get(i, cls.UNKNOWN)

UNKNOWN = UNKNOWN()


class Sample(object):

    def __init__(self, family_id, sample_id, paternal_id, maternal_id, sex,
                 phenotype, extra_attrs=None, header=None):
        self.family_id = family_id
        self.sample_id = sample_id
        self.dad = None
        self.mom = None
        self.paternal_id = paternal_id if paternal_id != '-9' else UNKNOWN
        self.maternal_id = maternal_id if maternal_id != '-9' else UNKNOWN
        if self.paternal_id != UNKNOWN:
            self.dad = self.paternal_id
        if self.maternal_id != UNKNOWN:
            self.mom = self.maternal_id
        self._sex = sex
        self._phenotype = phenotype
        self.sex = SEX.lookup(sex)
        self.affected = PHENOTYPE.lookup(phenotype)
        self.kids = []
        self.header = header
        self.attrs = []
        for a in (extra_attrs or []):
            self.attrs.append(a)

    def __eq__(self, other):
        if self is None or other is None:
            return False
        return (self.sample_id == other.sample_id) and (self.family_id ==
                                                        other.family_id)
    def __repr__(self):
        v = "%s('%s', '%s', '%s', '%s', '%s', '%s'" % (self.__class__.__name__,
                                                 self.family_id, self.sample_id,
                                                 self.paternal_id,
                                                 self.maternal_id,
                                                 SEX.rlookup(self.sex),
                                                 PHENOTYPE.rlookup(self.affected))
        if self.attrs:
            v += ", " + str(self.attrs)
        v += ")"
        return v


    @property
    def siblings(self):
        sibs = []
        for parent in (self.mom, self.dad):
            if parent is UNKNOWN: continue
            sibs.extend([x for x in parent.kids if not x in sibs and x != self])
        return sibs

    @property
    def full_siblings(self):
        sibs = []

        if self.mom is None or self.dad is None:
            return sibs
        return [s for s in self.mom.kids if s in self.dad.kids and s != self]

    @classmethod
    def from_row(cls, row, header=None):
        if isinstance(row, basestring):
            row = row.strip("\n").split()
        return cls(row[0], row[1], row[2], row[3], row[4], row[5], row[6:], header=header)

    def __str__(self):
        v = "%s %s %s %s %s %s" % (self.family_id, self.sample_id,
                                   self.paternal_id,
                                   self.maternal_id,
                                   self._sex, self._phenotype)
        if self.attrs:
            v += " " + " ".join(self.attrs)
        return v


class Family(object):
    def __init__(self, samples):
        assert len(set(s.family_id for s in samples)) == 1
        self.unknown_samples = []
        self.samples = samples
        self._build()
        for u in self.unknown_samples:
            print("unknown sample: %s in family: %s" % (u,
                samples[0].family_id), file=sys.stderr)

    def __iter__(self):
        self._index = 0
        return self

    def next(self):
        try:
            sample = self.samples[self._index]
        except IndexError:
            raise StopIteration
        self._index += 1
        return sample

    __next__ = next

    def summary(self):
        affection = {PHENOTYPE.AFFECTED: 0, PHENOTYPE.UNAFFECTED: 0,
                     PHENOTYPE.UNKNOWN: 0}
        sex = {SEX.MALE: 0, SEX.FEMALE: 0, SEX.UNKNOWN: 0}

        trios, quads = 0, 0
        generations = 1

        for s in self.samples:
            affection[s.affected] += 1
            sex[s.sex] += 1
            if s.mom and s.dad:
                trios += 1
            if list(s.full_siblings):
                quads += 1

        return affection, sex, trios, quads

    def _build(self):
        by_id = OrderedDict()
        for s in self.samples:
            by_id[s.sample_id] = s
        by_id[None] = None
        for s in self.samples:
            if s.paternal_id != UNKNOWN:
                try:
                    s.dad = by_id[s.paternal_id]
                    s.dad.kids.append(s)
                except KeyError:
                    s.dad = None
                    self.unknown_samples.append(s.paternal_id)
            if s.maternal_id != UNKNOWN:
                try:
                    s.mom = by_id[s.maternal_id]
                    s.mom.kids.append(s)
                except KeyError:
                    s.mom = None
                    self.unknown_samples.append(s.maternal_id)

    @property
    def affecteds(self):
        for s in self.samples:
            if s.affected:
                yield s

    @property
    def unaffecteds(self):
        for s in self.samples:
            if s.affected is False:
                yield s

    def trios(self, affected=True):
        for s in self.samples:
            if affected is not None:
                if affected != s.affected:
                    continue
            if s.mom and s.dad:
                yield (s, s.mom, s.dad)

class Ped(object):
    """Handle pedigree files

    >>> p = Ped('peddy/tests/a.ped')
    >>> p
    Ped('peddy/tests/a.ped')

    >>> s = next(p.samples())
    >>> s
    Sample('family_1', 'child_1', 'dad_1', 'mom_1', 'male', 'affected', ['caucasian'])

    >>> next(p.samples(phenotype=PHENOTYPE.UNAFFECTED))
    Sample('family_1', 'dad_1', '-9', '-9', 'male', 'unaffected', ['caucasian'])

    >>> next(p.samples(phenotype=PHENOTYPE.UNAFFECTED, sex=SEX.FEMALE))
    Sample('family_1', 'mom_1', '-9', '-9', 'female', 'unaffected', ['caucasian'])

    >>> [x.sample_id for x in p.samples(phenotype=PHENOTYPE.AFFECTED, sex=SEX.FEMALE)]
    ['kid_4', 'kid_4.1']

    #>>> p.summary()

    """

    def __init__(self, ped):
        if isinstance(ped, basestring):
            self.filename = ped
            ped = open(ped)
        else:
            self.filename = ""
        self._parse(ped)

    def _parse(self, fh):
        #header = None
        families = OrderedDict()

        for i, toks in enumerate(l.rstrip().split() for l in fh):
            if i == 0 and toks[0][0] == "#":
                header = toks
                continue
            sample = Sample.from_row(toks, header=header)
            if not sample.family_id in families:
                families[sample.family_id] = []

            families[sample.family_id].append(sample)

        self.families = OrderedDict()

        for family_id, samples in families.items():
            self.families[family_id] = Family(samples)

    def samples(self, phenotype=None, sex=None):

        if phenotype is None:
            samps = (x for fam in self.families.values() for x in fam)
        else:
            samps = (x for fam in self.families.values() for x in fam if x.affected == phenotype)

        if sex is not None:
            return (x for x in samps if x.sex == sex)

        return samps

    def __repr__(self):
        return "%s('%s')" % (self.__class__.__name__, self.filename)

    def relation(self, sample_a, sample_b):
        a = [x for x in self.samples() if x.sample_id == sample_a]
        assert len(a) == 1, (sample_a, "not found in ped file")
        b = [x for x in self.samples() if x.sample_id == sample_b]
        assert len(b) == 1, (sample_b, "not found in ped file")
        a, b = a[0], b[0]

        if a.mom == b or b.mom == a or a.dad == b or b.dad == a:
            return 'parent-child'

        if a.family_id != b.family_id:
            return 'unrelated'

        if a.mom is None and a.dad is None:
            return 'unrelated'

        if a.mom == b.mom and a.dad == b.dad:
            return 'full siblings'

        if (a.mom is not None and a.mom == b.mom) or (a.dad is not None and a.dad == b.dad):
            return 'siblings'

        else:
            return 'related level 2'

    def validate(self, vcf_path, plot=False):
        from cyvcf2 import VCF
        vcf = VCF(vcf_path, gts012=True, lazy=True)
        rels = list(vcf.relatedness(n_variants=9000, gap=25000))

        print("sample_1\tsample_2\tped_relation\tvcf_relation\trel\tIBS")
        for rel in rels:
            sample_a, sample_b = rel['pair']
            ped_rel = self.relation(sample_a, sample_b)
            out_line = "%s\t%s\t%s\t%s\t%.2f\t%.3f" % (sample_a, sample_b,
                    ped_rel, "|".join(rel['tags']), rel['rel'], rel['ibs'])
            if rel['rel'] < 0.04: # likely unrelated
                if ped_rel not in ('related level 2', 'unrelated'):
                    print(out_line)
                continue
            if rel['rel'] < 0.15:
                if ped_rel not in ('unrelated', 'related level 2', 'distant relations'):
                    print(out_line)
                continue
            if 0.26 < rel['rel'] < 0.78:
                if ped_rel not in ('parent-child', 'full siblings'):
                    print(out_line)
                continue

            if 0.15 < rel['rel'] < 0.3:
                if ped_rel not in ('related level 2', 'unrelated'):
                    print(out_line)
                continue

            if ped_rel > 0.78:
                if ped_rel not in ('identical twins', 'self'):
                    print(out_line)
                continue

    def summary(self):
        atrios, aquads = 0, 0
        for family in self.families:
            aff, sex, trios, quads = self.families[family].summary()
            atrios += trios
            aquads += quads
            d = locals()
            d['affected'] = aff[PHENOTYPE.AFFECTED]
            d['unaffected'] = aff[PHENOTYPE.AFFECTED]
            d['unknown'] = aff[PHENOTYPE.UNKNOWN]
            d['sex_unknown'] = sex[SEX.UNKNOWN]
            d['male'] = sex[SEX.MALE]
            d['female'] = sex[SEX.FEMALE]
            print("""\
===============================================
family: '{family}' :: {trios} trios and {quads} quads
===============================================
---------
phenotype
---------
affected: {affected}
unaffected: {unaffected}
unknown: {unknown}

---
sex
---
male: {male}
female: {female}
unknown: {sex_unknown}""".format(**d))

if __name__ == "__main__":
    import doctest
    doctest.testmod()

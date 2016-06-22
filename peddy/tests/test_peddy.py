from __future__ import print_function
import os
import os.path as op
import sys
from peddy import Ped, Family, Sample, PHENOTYPE, SEX

HERE = op.dirname(op.dirname(os.path.abspath(os.path.dirname(__file__))))


def test_sample():
    s = Sample('fam1', 'sample1', '-9', '-9', '2', '2')

    assert s.sex == SEX.FEMALE, (s.sex)
    assert s.affected == PHENOTYPE.AFFECTED
    assert s.kids == []

def test_sample_str_and_from_row():
    s = Sample('fam1', 'sample1', '-9', '-9', '2', '2')
    assert str(s) == "fam1 sample1 -9 -9 2 2", str(s)

    s2 = Sample.from_row(str(s))
    assert s2.sample_id == s.sample_id
    assert s2.sex == s.sex
    assert s2.family_id == s.family_id

def test_sex_check():
    if sys.version_info[0] == 3:
        return

    p = Ped(op.join(HERE, 'peddy/tests/test.mendel.ped'))
    df = p.sex_check(op.join(HERE, 'peddy/tests/test.mendel.vcf.gz'))

    assert "predicted_sex" in df.columns
    assert "ped_sex", df.columns
    assert "error" in df.columns

def test_dict():
    s = Sample('fam1', 'sample1', '-9', '-9', '2', '2')
    d = s.dict()
    assert d == {'maternal_id': '-9', 'paternal_id': '-9', 'sex': 'female',
            'family_id': 'fam1', 'phenotype': 'affected', 'sample_id': 'sample1'}, d

    s = Sample('fam1', 'sample1', 'dad', 'mom', '1', '1')
    d = s.dict()
    assert d == {'maternal_id': 'mom', 'paternal_id': 'dad', 'sex': 'male', 'family_id':
            'fam1', 'phenotype': 'unaffected', 'sample_id': 'sample1'}

    s = Sample('fam1', 'sample1', 'dad', 'mom', '-1', '-1')
    d = s.dict()
    assert d == {'maternal_id': 'mom', 'paternal_id': 'dad', 'sex': '-9',
            'family_id': 'fam1', 'phenotype': 'affected', 'sample_id':
            'sample1'}, d

def test_json():
    p = Ped(op.join(HERE, 'peddy/tests/test.mendel.ped'))
    json = p.to_json()
    #expected = '[{"maternal_id": "-9", "paternal_id": "-9", "sex": "male", "family_id": "CEPH1463", "phenotype": "affected", "sample_id": "NA12889"}, {"maternal_id": "-9", "paternal_id": "-9", "sex": "female", "family_id": "CEPH1463", "phenotype": "affected", "sample_id": "NA12890"}, {"maternal_id": "NA12890", "paternal_id": "NA12889", "sex": "male", "family_id": "CEPH1463", "phenotype": "affected", "sample_id": "NA12877"}]'
    # this test may fail if order of dicts is changed
    assert "CEPH1463" in json, json



def t_ped_check():
    try:
        import pandas as pd
        import cyvcf2
        cyvcf2
    except ImportError:
        return
    p = Ped(op.join(HERE, 'peddy/tests/test.mendel.ped'))
    v = p.ped_check(op.join(HERE, b'peddy/tests/test.mendel.vcf.gz'))
    assert isinstance(v, pd.DataFrame), v

    # remove samples
    f = list(p.families.values())[0]
    l = len(f.samples)
    s = f.samples[-1]
    f.samples = f.samples[:-1]
    assert l -1 == len(f.samples)
    v = p.ped_check(op.join(HERE, b'peddy/tests/test.mendel.vcf.gz'))
    assert isinstance(v, pd.DataFrame), v
    assert "ibs0" in v.columns

    # changed the sample id of a sample
    s.sample_id = "XDFSDFX"
    f.samples.append(s)
    v = p.ped_check(op.join(HERE, b'peddy/tests/test.mendel.vcf.gz'))
    assert isinstance(v, pd.DataFrame), v



def test_relation():
    kid = Sample('fam1', 'kid', 'dad', 'mom', '2', '2')
    dad = Sample('fam1', 'dad', '-9', '-9', '1', '2')
    mom = Sample('fam1', 'mom', '-9', '-9', '2', '2')
    kid.mom = mom
    kid.dad = dad

    from io import StringIO
    p = Ped(StringIO())
    p.families['fam1'] = Family([kid, mom, dad])
    assert p.relation("mom", "dad") == "mom-dad"

def test_relatedness_coefficient_missing_gparent():
    p = Ped(open(os.path.join(HERE, "peddy/tests/test.fam.ped")))
    # uncle
    v = p.relatedness_coefficient('101806-101806', '101811-101811')
    assert v == 0.25, v
    v = p.relatedness_coefficient('101806-101806', '101809-101809')
    assert v == 0.25, v
    # parent-child
    v = p.relatedness_coefficient('101806-101806', '101653-101653')
    assert v == 0.5, v

    p = Ped(open(os.path.join(HERE, "peddy/tests/test.fam2.ped")))
    v = p.relatedness_coefficient('101806-101806', '101811-101811')
    assert v == 0.25, v
    v = p.relatedness_coefficient('101806-101806', '101809-101809')
    assert v == 0.25, v

    # parent-child
    v = p.relatedness_coefficient('101806-101806', '101653-101653')
    assert v == 0.5, v



def test_relatedness_coefficient_missing_parent():

    gma = Sample('X28935', 'gma', '-9', '-9', '2', '1')
    mom = Sample('X28935', 'mom', '-9', 'gma', '2', '1')
    dad = Sample('X28935', 'dad', '-9', '-9', '1', '1')

    kid1 = Sample('X28935', 'kid1', '-9', 'mom', '1', '1')
    kid2 = Sample('X28935', 'kid2', '-9', 'mom', '2', '1')

    kid1 = Sample('X28935', 'kid1', 'dad', 'mom', '1', '1')
    kid2 = Sample('X28935', 'kid2', 'dad', 'mom', '2', '1')

    kid1.mom = mom
    kid2.mom = mom
    mom.mom = gma
    kid1.dad = dad
    kid2.dad = dad

    from io import StringIO
    p = Ped(StringIO())
    p.families['X28935'] = Family([kid1, kid2, mom, gma])#, dad])

    assert "siblings" in p.relation('kid1', 'kid2'), p.relation('kid1', 'kid2')

    v = p.relatedness_coefficient('kid1', 'kid2')
    assert v == 0.5, v

    v = p.relatedness_coefficient('gma', 'kid2')
    assert v == 0.25, v

    v = p.relatedness_coefficient('gma', 'kid1')
    assert v == 0.25, v

    v = p.relatedness_coefficient('gma', 'mom')
    assert v == 0.5, v

def test_relatedness_coefficient():
    kid = Sample('fam1', 'kid', 'dad', 'mom', '2', '2')
    dad = Sample('fam1', 'dad', '-9', '-9', '1', '2')
    mom = Sample('fam1', 'mom', '-9', '-9', '2', '2')
    gma = Sample('fam1', 'gma', '-9', '-9', '2', '2')
    ggma = Sample('fam1', 'ggma', '-9', '-9', '2', '2')
    kid.mom = mom
    kid.dad = dad
    mom.mom = gma
    gma.mom = ggma

    unrelated = Sample('fam1', 'un', '-9', '-9', '2', '2')

    from io import StringIO
    p = Ped(StringIO())
    p.families['fam1'] = Family([kid, mom, dad, gma, ggma, unrelated])
    rel = p.relatedness_coefficient("mom", "dad")
    assert rel == 0.0, rel
    d = p.relatedness_coefficient("mom", "kid")
    assert d == 0.5, d
    d = p.relatedness_coefficient("dad", "gma")
    assert d == 0.0, d

    d = p.relatedness_coefficient("mom", "gma")
    assert d == 0.5, d

    d = p.relatedness_coefficient("kid", "gma")
    assert d == 0.25, d

    d = p.relatedness_coefficient("kid", "ggma")
    assert d == 0.125, d

    assert p.relatedness_coefficient("mom", "mom") == 1.0

    #assert p.relatedness_coefficient("mom", "un") == 0.0


from contextlib import contextmanager

@contextmanager
def redirect_err(new_target=None):
    if new_target is None:
        try:
            from StringIO import StringIO
        except ImportError:
            from io import StringIO
        new_target = StringIO()

    old_target, sys.stderr = sys.stderr, new_target # replace sys.stdout
    try:
        yield new_target # run some code with the replaced stdout
    finally:
        sys.stdout = old_target # restore to the previous value

def test_warnings():

    with redirect_err() as out:
        kid = Sample('fam1', 'kid', 'dad', 'mom', '2', '2')
        mom = Sample('fam1', 'mom', '-9', '-9', '1', '2')
        dad = Sample('fam1', 'dad', '-9', '-9', '2', '2')
        kid.mom = mom
        kid.dad = dad

        v = out.getvalue()
        assert "'dad' is dad but has female sex" in v, v
        assert "'mom' is mom but has male sex" in v, v

    with redirect_err() as out:
        kid = Sample('fam1', 'kid', 'dad', 'mom', '2', '2')
        mom = Sample('fam1', 'mom', '-9', '-9', '-9', '2')
        kid.mom = mom
        v = out.getvalue()
        assert "'mom' is mom but has unknown sex. Setting to female" in v

    with redirect_err() as out:
        kid = Sample('fam1', 'kid', 'dad', 'mom', '2', '2')
        dad = Sample('fam1', 'dad', '-9', '-9', '-9', '2')
        kid.dad = dad
        v = out.getvalue()
        assert "'dad' is dad but has unknown sex. Setting to male" in v

    with redirect_err() as out:
        kid = Sample('fam1', 'kid', 'dad', 'mom', '2', '2')
        kid.dad = kid
        v = out.getvalue()
        assert "'kid' is dad of self" in v, v


def test_family():
    kid = Sample('fam1', 'kid', 'dad', 'mom', '2', '2')
    mom = Sample('fam1', 'mom', '-9', '-9', '2', '2')
    dad = Sample('fam1', 'dad', '-9', '-9', '1', '2')

    f = Family([kid, mom, dad])

    assert mom.kids == [kid]
    assert dad.kids == [kid]

    assert kid.dad == dad
    assert kid.mom == mom

    assert list(f.affecteds) == [kid, mom, dad], list(f.affecteds)
    assert list(f.unaffecteds) == []

    assert list(f) == [kid, mom, dad]


def test_trios():
    p = Ped(op.join(HERE, 'peddy/tests/a.ped'))
    f = p.families['family_4']
    trios = list(f.trios())
    assert len(trios) == 3

    assert [t[0] for t in trios] == list(f.affecteds)



def test_ped():

    p = Ped(op.join(HERE, 'peddy/tests/a.ped'))
    assert len(p.families) == 4

    assert len(list(p.samples())) == 14


def test_getattr():
    p = Ped(op.join(HERE, 'peddy/tests/a.ped'))
    li = list(p.samples(ethnicity='caucasianNEuropean'))
    assert len(li) == 5
    for item in li:
        assert item.ethnicity == 'caucasianNEuropean'

def test_6():
    p = Ped(op.join(HERE, 'peddy/tests/a6.ped'))
    assert len(list(p.samples())) == 14
    for sam in p.samples():
        assert sam.family_id[:3] == "fam"

def test_attrs():
    kid = Sample('fam1', 'kid', 'dad', 'mom', '2', '2', ['asdf', 'hello'])
    assert str(kid) == "fam1 kid dad mom 2 2 asdf hello", str(kid)
    assert repr(kid) == "Sample('fam1', 'kid', 'dad', 'mom', 'female', 'affected', ['asdf', 'hello'])", repr(kid)


def test_distant():

    p = Ped(op.join(HERE, 'peddy/tests/test-unknown-gma.ped'))

    d = p.relatedness_coefficient('kid1', 'cousin1')
    assert d == 0.125, d
    d = p.relatedness_coefficient('kid1', 'aunt')
    assert d == 0.25, d
    d = p.relatedness_coefficient('cousin1', 'aunt')
    assert d == 0.5, d
    d = p.relatedness_coefficient('mom', 'aunt')
    assert d == 0.5, d

    r = p.relation('kid1', 'cousin1')
    assert r == 'cousins', r

    r = p.relation('kid1', 'grandma')
    assert r == 'grandchild', r

    r = p.relation('kid1', 'aunt')
    assert r == 'niece/nephew', r

    # because we don't know that the uncle is related
    r = p.relation('kid1', 'uncle')
    assert r == 'related at unknown level', r

    r = p.relation('cousin1', 'mom')
    assert r == 'niece/nephew', r
    r = p.relation('cousin1', 'dad')
    # because we don't know that the dad is related
    assert r == 'related at unknown level', r

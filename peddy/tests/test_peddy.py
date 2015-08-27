from peddy import Ped, Family, Sample, PHENOTYPE, SEX

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


import sys
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
        assert "'dad' is dad but has female gender" in v, v
        assert "'mom' is mom but has male gender" in v, v

    with redirect_err() as out:
        kid = Sample('fam1', 'kid', 'dad', 'mom', '2', '2')
        mom = Sample('fam1', 'mom', '-9', '-9', '-9', '2')
        kid.mom = mom
        v = out.getvalue()
        assert "'mom' is mom but has unknown gender. Setting to female" in v

    with redirect_err() as out:
        kid = Sample('fam1', 'kid', 'dad', 'mom', '2', '2')
        dad = Sample('fam1', 'dad', '-9', '-9', '-9', '2')
        kid.dad = dad
        v = out.getvalue()
        assert "'dad' is dad but has unknown gender. Setting to male" in v

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
    p = Ped('peddy/tests/a.ped')
    f = p.families['family_4']
    trios = list(f.trios())
    assert len(trios) == 3

    assert [t[0] for t in trios] == list(f.affecteds)



def test_ped():

    p = Ped('peddy/tests/a.ped')
    assert len(p.families) == 4

    assert len(list(p.samples())) == 14


def test_getattr():
    p = Ped('peddy/tests/a.ped')
    li = list(p.samples(ethnicity='caucasianNEuropean'))
    assert len(li) == 5
    for item in li:
        assert item.ethnicity == 'caucasianNEuropean'


def test_attrs():
    kid = Sample('fam1', 'kid', 'dad', 'mom', '2', '2', ['asdf', 'hello'])
    assert str(kid) == "fam1 kid dad mom 2 2 asdf hello", str(kid)
    assert repr(kid) == "Sample('fam1', 'kid', 'dad', 'mom', 'female', 'affected', ['asdf', 'hello'])", repr(kid)


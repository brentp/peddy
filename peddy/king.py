from .utils import which
import os
import subprocess

def run_king(vcf_path, ped_obj):

    king = which('king')
    if king is None:
        raise Exception("requested king analysis, but king executable not found")
    plink = which('plink')
    if plink is None:
        raise Exception("requested king analysis, but plink executable not found")
    out = os.path.basename(vcf_path)
    if out.endswith(".gz"):
        out = out[:-3]
    if out.endswith(".vcf"): out = out[:-4]

    elog = open(out + ".plink.err", "w")
    olog = open(out + ".plink.log", "w")
    subprocess.check_call([plink, "--vcf", vcf_path, "--make-bed",
                           "--noweb",
                           "--out", out, "--biallelic-only", "--geno",
                           "0.95", "--vcf-half-call", "m"], stderr=elog,
                           stdout=olog)
    elog.close() and olog.close()

    elog = open(out + ".plink.err", "w")
    olog = open(out + ".plink.log", "w")
    subprocess.check_call([king, "-b", out + ".bed", "--kinship",
                           "--prefix", out], stderr=elog, stdout=olog)
    elog.close() and olog.close()


    sib_pairs, parent_child = [], []
    for fid, fam in ped_obj.families.items():
        sib_pairs.extend((x.sample_id, y.sample_id) for x, y in fam.sib_pairs)
        parent_child.extend((x.sample_id, y.sample_id) for x, y in fam.parent_child)

    king_pairs = read_king(out + ".kin0")
    kingped(ped_obj, king_pairs, sib_pairs, parent_child)


def read_king(king_file):
    pairs = {}
    import toolshed as ts
    for d in ts.reader(king_file):
        pairs[(d['ID1'], d['ID2'])] = float(d['Kinship'])
        pairs[(d['ID2'], d['ID1'])] = float(d['Kinship'])
    return pairs


def kingped(ped_obj, king_pairs, sib_pairs, parent_child, cutoff=0.13):

    seen = set()
    print("sample_a\tsample_b\terror\tkinship")
    fmt = "%s\t%s\t%s\t%.3f"

    for s in sib_pairs:
        try:
            king_pairs[s]
            seen.add(s)
        except KeyError:
            continue

        v = king_pairs[s]
        if v < cutoff:
            print(fmt % (s[0], s[1], "sibs low kingship", v))

    for kp in parent_child:
        try:
            king_pairs[kp]
            seen.add(kp)
        except KeyError:
            continue
        v = king_pairs[kp]
        if v < cutoff:
            print(fmt % (kp[0], kp[1], "parent-offspring low kinship", v))

    high = [((a, b), v) for (a, b), v in king_pairs.items() if v > cutoff and not ((a, b)
           in seen  or (b, a) in seen)]
    #print "high kinship (> %.3f), but not listed as such in ped file):" % cutoff
    pair_seen = {}
    for (a, b), v in high:
        if (a, b) in pair_seen: continue
        if (b, a) in pair_seen: continue
        try:
            ao, bo = ped_obj[a], ped_obj[b]
        except KeyError:
            continue
        if ao.family_id != bo.family_id:
            print(fmt % (a, b, "high kinship", v))
        pair_seen[(b, a)] = True

    pair_seen = {}
    for k, v in king_pairs.items():
        if k in pair_seen: continue
        if v > 0.42:
            print(fmt % (k[0], k[1], "possible duplicates", v))
            pair_seen[(k[1], k[0])] = True

if __name__ == "__main__":
    import sys
    kingped(sys.argv[1], sys.argv[2])


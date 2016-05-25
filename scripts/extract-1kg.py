import numpy as np
from cyvcf2 import VCF
from scipy.stats import chi2_contingency


def get_hwe_likelihood(obs_hom_ref, obs_het, obs_hom_alt, aaf):
    """
    Compute the likelihood of deviation from HWE using X^2,
    as well as the inbreeding coefficient.
    """
    # Bail out if aaf is undefined. This occurs
    # when there are multiple alternate alleles
    q = aaf
    p = 1 - aaf
    obs = [obs_hom_ref, obs_het, obs_hom_alt]
    n = sum(obs)
    exp = [p ** 2, 2 * p * q, q ** 2]
    exp = [n * v for v in exp]

    return chi2_contingency(np.array([exp, obs]))



if __name__ == "__main__":
    import sys

    #sites = [x.strip().split(":") for x in open('/uufs/chpc.utah.edu/common/home/u6000771/purcell5k.intervals')]
    sites = [x.strip().split(":") for x in open(sys.argv[1])]
    for s in sites:
        s[1] = int(s[1])

    hwes = []

    a = {}
    vcf = VCF('/uufs/chpc.utah.edu/common/home/u6000294/isilon/1kg_p3/up/ALL.phase3.autosome.vcf.gz', gts012=True)
    for i, (chrom, pos) in enumerate(sites):
        k = 0
        for v in vcf('%s:%s-%s' % (chrom, pos, pos)):

            if len(v.ALT) != 1: continue
            if len(v.REF) != 1: continue
            if '>' in v.ALT[0]: continue
            if v.call_rate < 0.8: continue
            if not (0.04 < v.aaf < 0.99): continue

            key = (v.CHROM, v.POS, v.REF, v.ALT[0])
            x = get_hwe_likelihood(v.num_hom_ref, v.num_het, v.num_hom_alt, v.aaf)
            if x[1] < 0.04: continue

            sites[i].extend((v.REF, v.ALT[0]))
            a[key] = np.array(v.gt_types)
            if k > 0:
                print sites[i]
            k += 1
            hwes.append(x[1])
            if len(hwes) % 500 == 0:
                print sites[i], len(hwes)

    print sites[i], len(hwes)

    from matplotlib import pyplot as plt
    plt.hist(hwes, 50)
    plt.show()


    a = np.array([a[tuple(s)] for s in sites if len(s) == 4], dtype=np.uint8)
    sites = [tuple(s) for s in sites if len(s) == 4]
    print(len(sites))
    print(a.shape)

    with open('1kg.sites.%d.%d.bin' % a.shape, 'w') as fh:
        fh.write(a.tostring())
    with open('1kg.sites', 'w') as fh:
        fh.write("\n".join("%s:%d:%s:%s" % tuple(s) for s in sites))
        fh.write("\n")


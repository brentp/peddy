
def main(vcf, ped, prefix, plot=False):
    from .pedagree import Ped

    ped = Ped(ped)

    siter = ped.sex_check(vcf, plot=("%s-sex-check.png" % prefix) if plot else False)
    with open("%s.sex-check.txt" % prefix, "w") as fh:
        for i, d in enumerate(siter):
            if i == 0:
                fh.write("\t".join(d.keys()) + "\n")
            d['het_ratio'] = "%.3f" % d['het_ratio']
            fh.write("\t".join(map(str, d.values())) + "\n")

    df = ped.ped_check(vcf)
    df.to_csv(open("%s.ped-check.txt" % prefix, "w"), sep="\t", index=False,
              float_format="%.3g")


if __name__ == "__main__":
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument("--plot", default=False, action="store_true")
    p.add_argument("--prefix", help="prefix for output files (default is basename of vcf)")
    p.add_argument("vcf", help="bgzipped and tabix'ed VCF")
    p.add_argument("ped", help="pedigree file")
    a = p.parse_args()
    prefix = a.prefix or a.vcf[:-6]
    main(a.vcf, a.ped, prefix, a.plot)

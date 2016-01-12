
def main(vcf, ped, prefix, plot=False):
    from .pedagree import Ped

    ped = Ped(ped)
    prefix = prefix.rstrip(".-")

    if plot:
        plot = prefix + ".het-check.png"
    ped.het_check(vcf, plot=plot).to_csv(prefix + ".het-check.csv", sep=",", index=False)
    if plot:
        plot = prefix + ".ped-check.png"
    ped.ped_check(vcf, plot=plot).to_csv(prefix + ".ped-check.csv", sep=",", index=False)
    if plot:
        plot = prefix + ".sex-check.png"
    ped.sex_check(vcf, plot=plot).to_csv(prefix + ".sex-check.csv", sep=",", index=False)


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

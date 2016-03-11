import multiprocessing as mp
from .pedagree import Ped

def run(args):
    check, pedf, vcf, plot, prefix = args
    p = Ped(pedf)
    print(check)
    if plot:
        plot = prefix + "." + check + "-check.png"

    getattr(p, check)(vcf, plot=plot).to_csv(prefix + "%s-check.csv" % check, sep=",", index=False)


def main(vcf, pedf, prefix, plot=False):

    ped = Ped(pedf)
    prefix = prefix.rstrip(".-")
    print("")

    p = mp.Pool(4)
    list(p.imap(run, [(check, pedf, vcf, plot, prefix) for check in ("het_check", "ped_check", "sex_check")]))

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

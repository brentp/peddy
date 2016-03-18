import multiprocessing as mp
import sys
from .pedagree import Ped

def run(args):
    check, pedf, vcf, plot, prefix = args
    # only print warnings for het_check
    p = Ped(pedf, warn=False)
    print(check)
    sys.stdout.flush()
    if plot:
        plot = prefix + "." + check + "-check.png"

    df = getattr(p, check)(vcf, plot=plot)
    df.to_csv(prefix + (".%s-check.csv" % check), sep=",", index=False)
    return (check, df.to_json(orient='records'))

def main(vcf, pedf, prefix, plot=False):

    ped = Ped(pedf)
    prefix = prefix.rstrip(".-")
    print("")

    p = mp.Pool(4)
    vals = {'pedigree': ped.to_json()}
    for check, json in p.imap(run, [(check, pedf, vcf, plot, prefix) for check
                                    in ("het_check", "ped_check", "sex_check")]):
        vals[check] = json
    sys.stdout.flush()
    p.close()
    print(vals)


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

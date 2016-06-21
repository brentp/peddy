import sys
import os.path as op
import string
from .peddy import Ped
import pandas as pd
from cyvcf2 import VCF
import io
import time

def run(args):
    check, pedf, vcf, plot, prefix, each, ncpus, sites = args
    # only print warnings for het_check
    p = Ped(pedf, warn=False)
    print("\033[1;31m%s\033[0m" % check)

    t0 = time.time()
    sys.stdout.flush()
    background_df = None
    if plot:
        plot = prefix + "." + check + ".png"

    if check in ("ped_check", "het_check"):
        kwargs = {'sites': sites} if check == 'ped_check' else {}
        df = getattr(p, check)(vcf, plot=plot, each=each, ncpus=ncpus,
                               prefix=prefix, **kwargs)
        if check == "het_check":
            df, background_df = df
    else:
        df = getattr(p, check)(vcf, plot=plot)

    df.to_csv(prefix + (".%s.csv" % check), sep=",", index=False, float_format="%.4g")
    if 'keep' in df.columns:
        df = df.ix[df['keep'], :]
        df.drop(['keep'], axis=1, inplace=True)
    d, unit = time.time() - t0, "seconds"
    if d > 100:
        d /= 60
        unit = "minutes"
    print "ran in %.1f %s" % (d, unit)
    if df.shape[0] > 50000 and check == "ped_check":
        # remove unknown relationships that aren't in error.
        df = df[((df.pedigree_relatedness != -1) &
                 (~df.parent_error) &
                 (~df.sample_duplication_error))]

    if check == "ped_check":
        # makes the plot nicer
        df.sort_values(inplace=True, by=["pedigree_relatedness"])

    return (check, df, background_df) #

def main(vcf, pedf, prefix, plot=False, each=1, ncpus=3, sites=None):

    tmpl = string.Template(open(op.join(op.dirname(__file__), "tmpl.html")).read())

    ped = Ped(pedf)
    prefix = prefix.rstrip(".-")
    print("")

    samples = VCF(vcf).samples

    ped_df = pd.read_table(pedf,
                           header=None,
                           names=ped.header or ['family_id', 'sample_id',
                                          'paternal_id', 'maternal_id',
                                          'sex', 'phenotype'],
                           # if there's a header, we skip it as it's inidcated
                           # above.
                           skiprows=1 if ped.header else None,
                           sep="\t")
    ped_df.family_id = ped_df.family_id.astype(basestring)
    ped_df.sample_id = ped_df.sample_id.astype(basestring)
    ped_df.index = ped_df.sample_id

    samples = set(samples)
    samples = [s for s in ped_df['sample_id'] if s in samples]

    ped_df = ped_df.ix[samples, :]


    keep_cols = {"ped_check": [],
                 "het_check": ["het_ratio", "mean_depth", "iqr_baf", "ancestry-prediction", "PC1", "PC2", "PC3"],
                 "sex_check": ["het_ratio", "error"]}

    vals = {'title': op.splitext(op.basename(pedf))[0], 'each': each}
    # background_df only present for het-check. It's the PC's from 1000G for
    # plotting.
    for check, df, background_df in map(run, [(check, pedf, vcf, plot, prefix, each, ncpus, sites) for check
                                 in ("ped_check", "het_check", "sex_check")]):
        vals[check] = df.to_json(orient='split' if check == "ped_check" else 'records', double_precision=3)

        if check != "ped_check":
            df.index = df['sample_id'].astype(basestring)
            for col in keep_cols[check]:
                c = check.split("_")[0] + "_"
                col_name = col if col.startswith(("PC", c)) else c + col
                ped_df[col_name] = list(df[col].ix[samples])
        if background_df is not None:
            vals["background_pca"] = background_df.to_json(orient='records', double_precision=3)

    new_pedf = prefix + ".peddy.ped"
    cols = list(ped_df.columns)
    cols[0] = '#' + cols[0]
    ped_df.columns = cols
    ped_df.to_csv(new_pedf, sep="\t", index=False, float_format="%.4g")

    # output the new version with the extra columns.
    # avoid extra stuff to stderr
    vals['pedigree'] = Ped(new_pedf, warn=False).to_json(samples)

    sys.stdout.flush()
    with open("%s.html" % prefix, "w") as fh:
        fh.write(tmpl.substitute(**vals))


if __name__ == "__main__":
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument("--plot", default=False, action="store_true")
    p.add_argument("-p", "--procs", type=int, default=3, help="number of processors to use")
    p.add_argument("--prefix", help="prefix for output files (default is basename of vcf)")
    p.add_argument("--each", help="sample every nth value from the selected sites instead of every value to speed processing.",
                   type=int, default=1)
    p.add_argument("--sites", help=r"This is rarely used. The path to a file with alternative sites to use for calculating relatedness in format 1:234234\n1:45345345\n..." +
                                " with chrom:pos[:ref:alt] on each line",
                   default=op.join(op.dirname(__file__), '1kg.sites'))
    p.add_argument("vcf", help="bgzipped and tabix'ed VCF")
    p.add_argument("ped", help="pedigree file")
    a = p.parse_args()
    prefix = a.prefix or a.vcf[:-6]
    main(a.vcf, a.ped, prefix, a.plot, a.each, a.procs, a.sites)

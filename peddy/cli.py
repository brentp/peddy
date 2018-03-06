import sys
import os.path as op
import string
import io
import time
from collections import defaultdict
import logging
from . import __version__

import coloredlogs
import click

import numpy as np
import pandas as pd

from .peddy import Ped
from cyvcf2 import VCF

from . import __version__

if sys.version_info[0] == 3:
    basestring = str

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
log = logging.getLogger(__name__)

def run(args):
    check, pedf, vcf, plot, prefix, each, ncpus, sites = args
    # only print warnings for het_check
    p = Ped(pedf, warn=False)
    log.info("\033[1;31m%s\033[0m" % check)

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

    if df is None:
        return check, None, None

    df.to_csv(prefix + (".%s.csv" % check), sep=",", index=False, float_format="%.4g")
    if 'keep' in df.columns:
        df = df.ix[df['keep'], :]
        df.drop(['keep'], axis=1, inplace=True)
    d, unit = time.time() - t0, "seconds"
    if d > 100:
        d /= 60
        unit = "minutes"
    log.info("ran in %.1f %s" % (d, unit))
    if df.shape[0] > 50000 and check == "ped_check":
        # remove unknown relationships that aren't in error.
        df = df[((df.pedigree_relatedness != -1) &
                 (~df.parent_error) &
                 (~df.sample_duplication_error))]

    if check == "ped_check":
        # makes the plot nicer
        df.sort_values(inplace=True, by=["pedigree_relatedness"])

    return (check, df, background_df) #


def correct_sex_errors(ped_df):
    excl = []
    if not 'sex_error' in ped_df.columns:
        return excl
    if np.any(ped_df['sex_error']):
        try:
            ped_df.pop(u'sex_fixed')
        except KeyError:
            pass
        ped_df.rename(columns={'sex_error': 'sex_fixed'}, inplace=True)

        sex_col = ped_df.columns[4]
        # now adjust the sex column as needed.
        sc = np.array(ped_df[sex_col])
        osc = np.array(ped_df[sex_col])
        sf = np.asarray(ped_df['sex_fixed'])
        gt = np.asarray(ped_df['sex_het_ratio'] > 0)
        if sc.dtype.char == 'S':
            sel = (gt & sf & (sc == '1'))
            sc[sel] = '2'
            if sel.sum():
                log.info("set sex of samples: %s to female in peddy.ped"
                         % ",".join(map(str, ped_df.ix[1, sel])))

            sel = (gt & sf & (sc == '2'))
            sc[sel] = '1'
            if sel.sum():
                log.info("set sex of samples: %s to male in peddy.ped"
                         % ",".join(map(str, ped_df.ix[1, sel])))
        else:
            for (ifrom, ito) in ((1, 2), (2, 1)):
                sel = (gt & sf & (sc == ifrom))
                osc[sel] = ito
                if sel.sum():
                    # Should this be a warning?
                    log.info("changed sex of samples: %s to %s in peddy.ped" % (
                          ",".join(map(str, ped_df.ix[sel, 1])),
                          ["", "male", "female"][ito]))


        ped_df[sex_col] = osc
    else:
        excl = ['sex_error']
    return excl


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('vcf')
@click.argument('ped')
@click.option("--plot",
    is_flag=True
)
@click.option("-p", "--procs",
    default=3,
    help="number of processors to use"
)
@click.option("--prefix",
    help="prefix for output files (default is basename of vcf)"
)
@click.option("--each",
    help="sample every nth value from the selected sites instead of every value"\
         " to speed processing.",
    default=1
)
@click.option("--sites",
    help=r"This is rarely used. The path to a file with alternative sites to"\
          " use for calculating relatedness in format 1:234234\n1:45345345..."\
          " with chrom:pos[:ref:alt] on each line",
    default=op.join(op.dirname(__file__), '1kg.sites')
)
@click.option('--loglevel',
    default='INFO',
    type=click.Choice(LOG_LEVELS),
    help="Set the level of log output.",
    show_default=True,
)
@click.version_option(version=__version__, prog_name="peddy")
def peddy(vcf, ped, plot, procs, prefix, each, sites, loglevel):
    """pleasingly pythonic pedigree manipulation"""
    coloredlogs.install(log_level=loglevel)
    log.info("Running Peddy version %s", __version__)
    prefix = prefix or vcf[:-6]

    tmpl = string.Template(open(op.join(op.dirname(__file__), "tmpl.html")).read())

    ped_obj = Ped(ped)
    prefix = prefix.rstrip(".-")

    samples = VCF(vcf).samples

    ped_df = pd.read_table(ped,
                           header=None,
                           names=ped_obj.header or ['family_id', 'sample_id',
                                          'paternal_id', 'maternal_id',
                                          'sex', 'phenotype'],
                           # if there's a header, we skip it as it's inidcated
                           # above.
                           index_col=False,
                           skiprows=1 if ped_obj.header else None,
                           sep="\t")

    ped_df.family_id = [str(x) for x in ped_df.family_id]
    ped_df.index = ped_df.sample_id = [str(x) for x in ped_df.sample_id]

    samples = set(samples)

    # exhaust list explicitly to get str from bytes.
    in_vcf_not_in_ped = samples - set(ped_df.index)
    in_ped_not_in_vcf = set(ped_df.index) - samples
    if in_vcf_not_in_ped:
        log.warning("%d samples in vcf not in ped:\n%s\n" %
                (len(in_vcf_not_in_ped), ",".join(in_vcf_not_in_ped)))
    if in_ped_not_in_vcf:
        log.warning("%d samples in ped not in vcf:\n%s\n" %
               (len(in_ped_not_in_vcf), ",".join(in_ped_not_in_vcf)))

    # keep order.
    samples = [s for s in ped_df['sample_id'] if s in samples]

    ped_df = ped_df.ix[samples, :]


    keep_cols = {"ped_check": [],
                 "het_check": ["call_rate", "het_ratio", "mean_depth", "idr_baf", "ancestry-prediction", "PC1", "PC2", "PC3"],
                 "sex_check": ["het_ratio", "error"]}

    vals = {'title': op.splitext(op.basename(ped))[0], 'each': each}
    # background_df only present for het-check. It's the PC's from 1000G for
    # plotting.
    for check, df, background_df in map(run, [(check, ped, vcf, plot, prefix, each, procs, sites) for check
                                 in ("ped_check", "het_check", "sex_check")]):
        if df is None:
            vals[check] = []
            continue
        columns = df.columns
        if check == "sex_check" and not np.any(df["error"]):
            columns = [c for c in columns if not "error" == c]
        vals[check] = df[columns].to_json(orient='split' if check == "ped_check" else 'records', double_precision=3)

        if check != "ped_check":
            df.index = [str(s) for s in df['sample_id']]
            for col in keep_cols[check]:
                c = check.split("_")[0] + "_"
                col_name = col if col.startswith(("PC", c, "ancestry")) else c + col
                ped_df[col_name] = list(df[col].ix[samples])
        elif any(df['sample_duplication_error']):
            da = df.ix[df['sample_duplication_error'], 'sample_a']
            db = df.ix[df['sample_duplication_error'], 'sample_b']
            d = defaultdict(list)
            for a, b in zip(da, db):
                d[a].append(b)
                d[b].append(a)
            d = dict(d)
            for k in d:
                d[k] = ",".join(d[k])

            ped_df['duplicates'] = [d.get(s, "") for s in samples]
        if check == "ped_check":
            obs_vs_expected(prefix, df[['sample_a', 'sample_b', 'pedigree_relatedness', 'rel', 'n']])

        if background_df is not None:
            vals["background_pca"] = background_df.to_json(orient='records', double_precision=3)
            with open("%s.background_pca.json" % prefix, "w") as fh:
                fh.write(vals["background_pca"])

    ped_df.paternal_id = ped_df.paternal_id.astype(str)
    ped_df.maternal_id = ped_df.maternal_id.astype(str)
    ped_df.loc[ped_df.paternal_id == "", "paternal_id"] = "-9"
    ped_df.loc[ped_df.maternal_id == "", "maternal_id"] = "-9"
    ped_df.loc[ped_df.paternal_id == "nan", "paternal_id"] = "-9"
    ped_df.loc[ped_df.maternal_id == "nan", "maternal_id"] = "-9"

    ped_df[ped_df.paternal_id == ""] = "-9"
    ped_df[ped_df.maternal_id == ""] = "-9"

    new_pedf = prefix + ".peddy.ped"
    cols = list(ped_df.columns)
    cols[0] = '#' + cols[0]
    ped_df.columns = cols

    ped_df.to_csv(new_pedf, sep="\t", index=False, mode='w', float_format="%.4g")

    # output the new version with the extra columns.
    # avoid extra stuff to stderr
    vals['pedigree'] = Ped(new_pedf, warn=False).to_json(samples, exclude=('PC1', 'PC2', 'PC3'))
    import os
    vals['prefix'] = os.path.basename(prefix)

    excl = correct_sex_errors(ped_df)

    ped_df[[c for c in ped_df.columns if not c in excl]].to_csv(new_pedf, sep="\t", index=False, mode='w', float_format="%.4g")


    sys.stdout.flush()
    with open("%s.html" % prefix, "w") as fh:
        fh.write(tmpl.substitute(**vals))

def obs_vs_expected(prefix, df):
    tmpl = string.Template(open(op.join(op.dirname(__file__), "tmpl.vs.html")).read())
    vals = {'title': 'observed vs expected relatedness'}
    vals['data'] = df.to_json(orient='records', double_precision=3)
    with open("%s.vs.html" % prefix, "w") as fh:
        fh.write(tmpl.substitute(**vals))

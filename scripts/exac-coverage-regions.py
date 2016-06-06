import gzip
import sys

def main(cov, pct, directory, patt="Panel.chr%s.coverage.txt.gz"):
    for chrom in range(1, 23) + ["X", "Y"]:
        run_chrom(cov, pct, patt % chrom)
        sys.stdout.flush()

def run_chrom(cov, pct, file):

    fh = gzip.open(file)
    header = next(fh).strip()[1:].split("\t")
    cov = str(cov)

    start = None
    end = -1
    for d in (dict(zip(header, l.rstrip().split("\t"))) for l in fh):
        pos = int(d['pos'])
        d[cov] = float(d[cov])
        # end a region.
        if d[cov] < pct or pos - end > 1:
            if start is not None:
                print "{chrom}\t{start}\t{end}".format(start=start, chrom=d['chrom'], end=end)
            start = None

        # start a new region.
        elif start is None:
            start = pos - 1

        # always increment the end.
        end = pos

    if start is not None:
        print "{chrom}\t{start}\t{end}".format(start=start, chrom=d['chrom'], end=end)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--coverage", choices=(1, 5, 10, 15, 20, 25, 30, 50, 100), required=True, type=int)
    p.add_argument("--percent", type=float, default=80, help="percent of samples having --coverage")

    p.add_argument("-d", "--directory", default="./",
                   help="path to directory with ExAC coverage files")

    a = p.parse_args()

    if a.percent > 1: a.percent /= 100.0
    main(a.coverage, a.percent, a.directory)

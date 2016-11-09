import sys
import toolshed as ts
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def shared(fs):
    sets = []
    for f in fs:
        s = set()
        for d in ts.reader(f, sep=","):
            x = float(d['rel'])
            y = float(d['pedigree_relatedness'])
            #if abs(x - y) > 0.2: continue
            s.add((d['sample_a'], d['sample_b']))
        sets.append(s)
    sall = sets[0]
    for i, s in enumerate(sets):
        if i == 0: continue
        sall &= s
    return sall


def plot(f, axs, shared):

    diffs = []

    xs, ys = [], []
    for d in ts.reader(f, sep=","):
        if not (d['sample_a'], d['sample_b']) in shared: continue

        x = float(d['rel'])
        y = float(d['pedigree_relatedness'])
        #if abs(x - y) > 0.25: continue
        diffs.append(x - y)
        xs.append(x)
        ys.append(y)
    """
    ax.scatter(xs, ys)
    ax.set_xlabel('relatedness by genotype')
    ax.set_ylabel('relatedness by ped file')
    ax.set_title(f)
    """
    p5, p95 = np.percentile(diffs, [2.5, 97.5])
    m, std = np.mean(diffs), np.std(diffs)

    ax2 = axs
    ax2.set_title(convert(f))
    ax2.hist(diffs, 40)
    ax2.text(0.6, 0.8, "95%% range: %.3f - %.3f\nmean: %.3f std: %-3f" % (p5, p95, m, std),
                    transform=ax2.transAxes)
    ax2.set_xlabel("genotype - expected")
    ax2.set_ylabel("count")

def convert(f):
    if ".geo." in f:
        return "geometric mean"
    if ".min." in f:
        return "minimum"
    return "arithmetic mean"

files = sys.argv[1:]
fig, axes = plt.subplots(len(files), 1, figsize=(8, 12))
try:
    axes[0]
except TypeError:
    axes = (axes,)

pairs = shared(files)

for i, f in enumerate(files):
    plot(f, axes[i], pairs)

plt.tight_layout()
plt.show()

import sys

import toolshed as ts
import numpy as np
import statsmodels.api as sm


king = list(ts.reader(1))
peddy = list(ts.reader(2, sep=","))


king_ibs = {tuple(sorted([d['ID1'], d['ID2']])): float(d['N_IBS0']) for d in king}
king_rel = {tuple(sorted([d['ID1'], d['ID2']])): 2.0 * float(d['Kinship']) for d in king}

try:
    actual_rel = {tuple(sorted([d['sample_a'], d['sample_b']])): float(d['pedigree_relatedness']) for d in peddy}
except KeyError:
    actual_rel = None

peddy_ibs = {tuple(sorted([d['sample_a'], d['sample_b']])): float(d['ibs0']) for d in peddy}
peddy_rel = {tuple(sorted([d['sample_a'], d['sample_b']])): float(d['rel']) for d in peddy}

keys = king_ibs.keys()
if set(king_ibs.keys()).symmetric_difference(peddy_ibs.keys()):
    sys.stderr.write("WARNING: not all the same samples")
    keys = list(set(keys).intersection(peddy_ibs.keys()))


from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('ticks')



values = {}
values['king'], values['peddy'] = [king_rel[k] for k in keys], [peddy_rel[k] for k in keys]
if not actual_rel:

    plt.scatter(values['king'], values['peddy'])
    plt.xlabel('king')
    plt.ylabel('other')
    xs, ys = plt.xlim(), plt.ylim()
    lims = [np.min([xs, ys]), np.max([xs, ys])]
    plt.plot(lims, lims, ls='-', c='0.75', lw=3, zorder=0)
    plt.show()
    sys.exit()




fig, axes = plt.subplots(ncols=4, figsize=(15, 5))

values['actual'] = [actual_rel[k] for k in keys]

xys = [('actual', 'king'), ('actual', 'peddy')]

for i in range(2):
    ax = axes[i]

    xn, yn = xys[i]

    ax.scatter(values[xn], values[yn])
    ax.set_xlabel('%s relatedness' % xn)
    ax.set_ylabel('%s relatedness' % yn)
    plt.tight_layout()
    xs, ys = ax.get_xlim(), ax.get_ylim()
    lims = [np.min([xs, ys]), np.max([xs, ys])]
    ax.plot(lims, lims, ls='-', c='0.75', lw=3, zorder=0)
    ax.set_xlim(xs)
    ax.set_ylim(ys)
    fmt = "R-squared: %.3f\np-value: %.3g"

    fit =  sm.OLS(values[yn], values[xn]).fit()
    print(fit.summary())
    print "\n\n\n"
    slope = fit.params[0]
    #ax.text(0.2, 0.8, fmt % (fit.rsquared, fit.f_pvalue),
    #        transform=ax.transAxes)




axes[2].scatter([king_ibs[k] for k in keys], [peddy_ibs[k] for k in keys])
axes[2].set_xlabel('king ibs0')
axes[2].set_ylabel('peddy ibs0')
fit = sm.OLS([king_ibs[k] for k in keys], [peddy_ibs[k] for k in keys]).fit()
axes[2].text(0.2, 0.8, fmt % (fit.rsquared, fit.f_pvalue),
        transform=axes[2].transAxes)

king_ibs = {tuple(sorted([d['ID1'], d['ID2']])): float(d['N_IBS2']) for d in king}
peddy_ibs = {tuple(sorted([d['sample_a'], d['sample_b']])): float(d['ibs2']) for d in peddy}
axes[3].scatter([king_ibs[k] for k in keys], [peddy_ibs[k] for k in keys])
axes[3].set_xlabel('king ibs2')
axes[3].set_ylabel('peddy ibs2')
fit = sm.OLS([king_ibs[k] for k in keys], [peddy_ibs[k] for k in keys]).fit()
axes[3].text(0.2, 0.8, fmt % (fit.rsquared, fit.f_pvalue),
        transform=axes[3].transAxes)

plt.tight_layout()
plt.show()

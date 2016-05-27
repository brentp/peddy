import toolshed as ts
import numpy as np
import statsmodels.api as sm


king = list(ts.reader(1))
peddy = list(ts.reader(2, sep=","))


king_ibs = {tuple(sorted([d['ID1'], d['ID2']])): float(d['N_IBS0']) for d in king}
king_rel = {tuple(sorted([d['ID1'], d['ID2']])): 2.0 * float(d['Kinship']) for d in king}

actual_rel = {tuple(sorted([d['sample_a'], d['sample_b']])): float(d['pedigree_relatedness']) for d in peddy}

peddy_ibs = {tuple(sorted([d['sample_a'], d['sample_b']])): float(d['ibs0']) for d in peddy}
peddy_rel = {tuple(sorted([d['sample_a'], d['sample_b']])): float(d['rel']) for d in peddy}

assert not set(king_ibs.keys()).symmetric_difference(peddy_ibs.keys())

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('ticks')

keys = king_ibs.keys()

fig, axes = plt.subplots(ncols=3, figsize=(15, 5))

values = {}
values['king'], values['peddy'] = [king_rel[k] for k in keys], [peddy_rel[k] for k in keys]
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
axes[2].set_xlabel('king ibs')
axes[2].set_ylabel('peddy ibs')
fit = sm.OLS([king_ibs[k] for k in keys], [peddy_ibs[k] for k in keys]).fit()
axes[2].text(0.2, 0.8, fmt % (fit.rsquared, fit.f_pvalue),
        transform=axes[2].transAxes)
plt.tight_layout()
plt.show()

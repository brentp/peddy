import toolshed as ts
from scipy.stats import linregress


king = list(ts.reader(1))
peddy = list(ts.reader(2, sep=","))


king_ibs = {tuple(sorted([d['ID1'], d['ID2']])): float(d['N_IBS0']) for d in king}
king_rel = {tuple(sorted([d['ID1'], d['ID2']])): 2.0 * float(d['Kinship']) for d in king}

peddy_ibs = {tuple(sorted([d['sample_a'], d['sample_b']])): float(d['ibs0']) for d in peddy}
peddy_rel = {tuple(sorted([d['sample_a'], d['sample_b']])): float(d['rel']) for d in peddy}

assert not set(king_ibs.keys()).symmetric_difference(peddy_ibs.keys())

from matplotlib import pyplot as plt
import seaborn as sns

keys = king_ibs.keys()

fig, axes = plt.subplots(2)

axes[0].scatter([king_rel[k] for k in keys], [peddy_rel[k] for k in keys])
axes[0].set_xlabel('king relatedness')
axes[0].set_ylabel('peddy relatedness')

res = linregress([king_rel[k] for k in keys], [peddy_rel[k] for k in keys])
axes[0].set_title("r^2: %.3f  p-value: %.3g" % (res[2]**2, res[3]))

axes[1].scatter([king_ibs[k] for k in keys], [peddy_ibs[k] for k in keys])
axes[1].set_xlabel('king ibs')
axes[1].set_ylabel('peddy ibs')
res = linregress([king_ibs[k] for k in keys], [peddy_ibs[k] for k in keys])
axes[1].set_title("r^2: %.3f  p-value: %.3g" % (res[2]**2, res[3]))
plt.tight_layout()
plt.show()

"""
plot the first 2 principal components of the 23356 sites that we use in peddy
and show that we nicely recover the population structure even with this small
number of markers.
The PCA runs in ~20 seconds and the randomized PCA runs in ~4 seconds.
"""

import gzip
import numpy as np
import toolshed as ts

import seaborn as sns
pal = sns.color_palette('Set1')

from sklearn.decomposition import PCA, RandomizedPCA

f = "../cyvcf2/cyvcf2/1kg.sites.bin.gz"

tmp = np.fromstring(gzip.open(f).read(), dtype=np.uint8).astype(np.int32)

genos = tmp.reshape((23556, len(tmp) / 23556)).T

clf = RandomizedPCA(n_components=2, whiten=True)

import time
t0 = time.time()
gf = clf.fit_transform(genos)
print time.time() - t0

pops = "AFR AMR EAS EUR SAS".split()

# https://hangouts.google.com/webchat/u/0/frame?v=1463175081&hl=en-US&pvt=AMP3uWbeJzucAfSaxyKZreeU3oew-CyaV-gwRMAvHpeHN9VKU_EhC1Fj75C7UKU-cpIk6HvlXK6K&prop=aChromeExtension#zSoyz
colors = np.array([x['super_pop'] for x in ts.reader('integrated_call_samples_v3.20130502.ALL.panel')])

print gf.shape

from matplotlib import pyplot as plt

for i, p in enumerate(pops):
    subset = colors == p

    plt.scatter(gf[subset, 0], gf[subset, 1], c=pal[i], label=p)
plt.legend()
plt.show()


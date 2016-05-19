"""
plot the first 2 principal components of the 23356 sites that we use in peddy
and show that we nicely recover the population structure even with this small
number of markers.
The PCA runs in ~20 seconds and the randomized PCA runs in ~4 seconds.
"""

import gzip
import numpy as np
import toolshed as ts
from sklearn import svm

import seaborn as sns
pal = sns.color_palette('Set1')

from sklearn.decomposition import PCA, RandomizedPCA
from sklearn.pipeline import make_pipeline

f = "../cyvcf2/cyvcf2/1kg.sites.bin.gz"

tmp = np.fromstring(gzip.open(f).read(), dtype=np.uint8).astype(np.int32)

genos = tmp.reshape((23556, len(tmp) / 23556)).T

clf = make_pipeline(RandomizedPCA(n_components=4, whiten=True, copy=False),
                    svm.SVC(C=2, probability=True))


ipops = "AFR AMR EAS EUR SAS".split()

# https://hangouts.google.com/webchat/u/0/frame?v=1463175081&hl=en-US&pvt=AMP3uWbeJzucAfSaxyKZreeU3oew-CyaV-gwRMAvHpeHN9VKU_EhC1Fj75C7UKU-cpIk6HvlXK6K&prop=aChromeExtension#zSoyz
pops = np.array([x['super_pop'] for x in ts.reader('integrated_call_samples_v3.20130502.ALL.panel')])
target = np.array([ipops.index(x) for x in pops])

#print "|".join(map(str, target))
#1/0
import time

t0 = time.time()
clf.fit(genos, target)
print "time to fit:", time.time() - t0


import cPickle

t0 = time.time()
cPickle.dump(clf, gzip.open('clf.pkl.gz', 'w'))
print time.time() - t0


del clf

clf = cPickle.load(gzip.open('clf.pkl.gz'))
rpca = clf.named_steps['randomizedpca']
gf = rpca.transform(genos)


"""
this code was used to validate the classifier.

from sklearn import cross_validation
cv = cross_validation.ShuffleSplit(len(target), n_iter=20, test_size=0.3, random_state=0)
print cross_validation.cross_val_score(clf, genos, target, cv=cv)
1/0

"""


print gf.shape

from matplotlib import pyplot as plt
fig, axes = plt.subplots(2)

for i, p in enumerate(ipops):
    subset = p == pops

    axes[0].scatter(gf[subset, 0], gf[subset, 1], c=pal[i], label=p)
    axes[1].scatter(gf[subset, 0], gf[subset, 2], c=pal[i], label=p)

axes[0].set_xlabel("PC1")
axes[0].set_ylabel("PC2")
axes[1].set_xlabel("PC1")
axes[1].set_ylabel("PC3")
plt.legend()
plt.show()


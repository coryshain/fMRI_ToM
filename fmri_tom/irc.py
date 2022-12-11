import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off'
)


def bootstrap(group1, group2=None, niter=10000, fisher=True):
    group1 = np.array(group1)
    if fisher:
        group1 = np.arctanh(group1)
    hits = np.zeros(niter)
    if group2 is None:  # Bootstrap test against 0
        val = group1.mean()
        for s in range(niter):
            hits[s] = np.random.choice(group1, size=len(group1), replace=True).mean() <= 0
    else:  # Permutation test
        group2 = np.array(group2)
        if fisher:
            group2 = np.arctanh(group2)
        val = group1.mean() - group2.mean()
        absval = np.abs(val)
        splitix = len(group1)
        v = np.concatenate([group1, group2], axis=0)
        for s in range(niter):
            samp = np.random.permutation(v)
            a = samp[:splitix]
            b = samp[splitix:]
            hits[s] = absval <= np.abs(a.mean() - b.mean())

    p = (hits.sum() + 1) / (niter + 1)

    return val, p


lang_rois = ['IFGorb', 'IFG', 'MFG', 'AntTemp', 'PostTemp']

irc_rest = pd.read_csv('../../data/fMRI_alicemultiling/AliceData_Table4_IRCs_Rest_data.csv', index_col=0)
irc_rest['Expt'] = 'story'
irc_stories = pd.read_csv('../../data/fMRI_alicemultiling/AliceData_Table3_Alice_IRCs_Alice_data.csv')
irc_stories['Expt'] = 'rest'

irc = pd.concat([irc_rest, irc_stories], axis=0).reset_index(drop=True)
irc = irc[irc.Parcel1.isin(lang_rois + ['AngG']) & irc.Parcel2.isin(lang_rois + ['AngG'])]

p1, p2, h1, h2 = irc.Parcel1, irc.Parcel2, irc.Hemi1, irc.Hemi2

comparison = np.zeros_like(irc.Parcel1)
comparison[p1.isin(lang_rois) & p2.isin(lang_rois) & (h1 == 'LH') & (h2 == 'LH')] = 'LangLH-LangLH'
comparison[p1.isin(lang_rois) & p2.isin(lang_rois) & (h1 == 'RH') & (h2 == 'RH')] = 'LangRH-LangRH'
comparison[p1.isin(lang_rois) & p2.isin(lang_rois) & (h1 == 'LH') & (h2 == 'RH')] = 'LangLH-LangRH'
comparison[p1.isin(lang_rois) & p2.isin(lang_rois) & (h1 == 'RH') & (h2 == 'LH')] = 'LangLH-LangRH'
comparison[(p1 == 'AngG') & (p2 == 'AngG')] = 'AngG-AngG'
comparison[p1.isin(lang_rois) & (p2 == 'AngG') & (h1 == 'LH')] = 'LangLH-AngG'
comparison[(p1 == 'AngG') & p2.isin(lang_rois) & (h2 == 'LH')] = 'LangLH-AngG'
comparison[p1.isin(lang_rois) & (p2 == 'AngG') & (h1 == 'RH')] = 'LangRH-AngG'
comparison[(p1 == 'AngG') & p2.isin(lang_rois) & (h2 == 'RH')] = 'LangRH-AngG'

irc['comparison'] = comparison

# Plot correlation heatmap

rois = [
    ('IFGorb', 'LH'),
    ('IFG', 'LH'),
    ('MFG', 'LH'),
    ('AntTemp', 'LH'),
    ('PostTemp', 'LH'),
    ('IFGorb', 'RH'),
    ('IFG', 'RH'),
    ('MFG', 'RH'),
    ('AntTemp', 'RH'),
    ('PostTemp', 'RH'),
    ('AngG', 'LH'),
    ('AngG', 'RH'),
]

R = np.zeros((12, 12))
for key, df in irc.groupby(['Parcel1', 'Hemi1', 'Parcel2', 'Hemi2']):
    key1 = key[:2]
    key2 = key[2:]
    i = rois.index(key1)
    j = rois.index(key2)
    if j > i:
        i, j = j, i
    r = np.tanh(np.arctanh(df.Corr.values).mean())
    # r = df.Corr.mean()
    R[i, j] = r
    R[j, i] = r

for i in range(12):
    R[i,i] = 1

f, ax = plt.subplots(figsize=(6, 5))

# cmap = sns.diverging_palette(230, 20, s=90, as_cmap=True)
# cmap = sns.color_palette("Spectral_r", as_cmap=True)
# cmap = 'YlGnBu'
# cmap = 'hot_r'
cmap = 'bone_r'

# mask = np.triu(np.ones_like(R, dtype=bool), k=1)
mask = np.zeros_like(R, dtype=bool)
sns.heatmap(R, mask=mask, cmap=cmap, vmin=0, vmax=1., center=0.5,
            square=True, linewidths=3, cbar_kws={"shrink": .75})
ax.set_xticks([])
ax.set_yticks([])
ax.tick_params(left=False, bottom=False)

plt.tight_layout()
plt.savefig('alice_corr.png', dpi=300)

plt.close('all')



# Plot correlation barplot

plt.rcParams["font.family"] = "Arial"

means = []
sems = []
names = [
    'LangLH-LangLH',
    'LangRH-LangRH',
    'AngG-AngG',
    'LangLH-LangRH',
    'LangLH-AngG',
    'LangRH-AngG',
]
ticks = [
    'LH-LH',
    'RH-RH',
    'AG-AG',
    'LH-RH',
    'LH-AG',
    'RH-AG'
]

for name in names:
    df = irc[irc.comparison == name]
    means.append(df.Corr.mean())
    sems.append(df.Corr.sem())

color = [
    (0, 0, 0),
    (20, 20, 20),
    (40, 40, 40),
    (60, 60, 60),
    (80, 80, 80),
    (100, 100, 100),
]

color = [(np.array(c) + 50) / 255 for c in color]
print(color)

f, ax = plt.subplots(figsize=(3, 5))

for i in range(len(names)):
    ax.bar(i, means[i], yerr=sems[i], color=color[i], ecolor=color[i], linewidth=2, capsize=4)


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(top='off', bottom='off', left=False, right='off', labelleft='on', labelbottom='on')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, which='major', axis='y', ls='--', lw=.5, c='k', alpha=.3)
ax.axhline(y=0, lw=1, c='gray', alpha=1)
ax.set_xticks(np.arange(len(names)), labels=ticks, rotation=90)

plt.tight_layout()
plt.savefig('alice_barr.png', dpi=300)
plt.close('all')



# Test networks against 0

b, p = bootstrap(irc[irc.comparison == 'LangLH-LangLH'].Corr)
print('LangLH-LangLH: mean = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangRH-LangRH'].Corr)
print('LangRH-LangRH: mean = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'AngG-AngG'].Corr)
print('AngG-AngG: mean = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangLH-LangRH'].Corr)
print('LangLH-LangRH: mean = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangLH-AngG'].Corr)
print('LangLH-AngG: mean = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangRH-AngG'].Corr)
print('LangRH-AngG: mean = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangLH-LangLH'].Corr, group2=irc[irc.comparison == 'LangRH-LangRH'].Corr)
print('LangLH-LangLH vs. LangRH-LangRH: diff = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangLH-LangLH'].Corr, group2=irc[irc.comparison == 'AngG-AngG'].Corr)
print('LangLH-LangLH vs. AngG-AngG: diff = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangRH-LangRH'].Corr, group2=irc[irc.comparison == 'AngG-AngG'].Corr)
print('LangRH-LangRH vs. AngG-AngG: diff = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangLH-LangLH'].Corr, group2=irc[irc.comparison == 'LangLH-LangRH'].Corr)
print('LangLH-LangLH vs. LangLH-LangRH: diff = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangLH-LangLH'].Corr, group2=irc[irc.comparison == 'LangLH-AngG'].Corr)
print('LangLH-LangLH vs. LangLH-AngG: diff = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangRH-LangRH'].Corr, group2=irc[irc.comparison == 'LangLH-LangRH'].Corr)
print('LangRH-LangRH vs. LangLH-LangRH: diff = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangRH-LangRH'].Corr, group2=irc[irc.comparison == 'LangRH-AngG'].Corr)
print('LangRH-LangRH vs. LangRH-AngG: diff = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'AngG-AngG'].Corr, group2=irc[irc.comparison == 'LangLH-AngG'].Corr)
print('AngG-AngG vs. LangLH-AngG: diff = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'AngG-AngG'].Corr, group2=irc[irc.comparison == 'LangRH-AngG'].Corr)
print('AngG-AngG vs. LangRH-AngG: diff = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangLH-LangRH'].Corr, group2=irc[irc.comparison == 'LangLH-AngG'].Corr)
print('LangLH-LangRH vs. LangLH-Ang: diff = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangLH-LangRH'].Corr, group2=irc[irc.comparison == 'LangRH-AngG'].Corr)
print('LangLH-LangRH vs. LangLH-Ang: diff = %f | p = %f' % (b, p))

b, p = bootstrap(irc[irc.comparison == 'LangLH-AngG'].Corr, group2=irc[irc.comparison == 'LangRH-AngG'].Corr)
print('LangLH-AngG vs. LangRH-Ang: diff = %f | p = %f' % (b, p))


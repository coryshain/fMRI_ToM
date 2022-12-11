import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams["font.family"] = "Arial"

lang_lh_core = [
    'LIFGorb',
    'LIFG',
    'LMFG',
    'LAntTemp',
    'LPostTemp',
]

lang_rh = [
    'RIFGorb',
    'RIFG',
    'RMFG',
    'RAntTemp',
    'RPostTemp',
]

lang_angg = [
    'LAngG',
    'RAngG',
]

lang_order = lang_lh_core + lang_rh + lang_angg

lang_fROIs = {
    1: 'LIFGorb',
    2: 'LIFG',
    3: 'LMFG',
    4: 'LAntTemp',
    5: 'LPostTemp',
    6: 'LAngG',
    7: 'RIFGorb',
    8: 'RIFG',
    9: 'RMFG',
    10: 'RAntTemp',
    11: 'RPostTemp',
    12: 'RAngG',
}

lang = pd.read_csv('../../data/fMRI_ToM/spcorr/spcorr_langlocSN_n149.csv')
lang.ROI = lang.ROI.map(lang_fROIs)
lang_mean = lang.groupby('ROI')['Fisher transformed correlation coefficients'].mean()
lang_sem = lang.groupby('ROI')['Fisher transformed correlation coefficients'].sem()

lang_mean = lang_mean[lang_order]
lang_sem = lang_sem[lang_order]

ToM_fROIs = {
    1: 'LTPJ',
    2: 'LDMPFC',
    3: 'LMMPFC',
    4: 'LVMPFC',
    5: 'LPC',
    6: 'RTPJ',
    7: 'RDMPFC',
    8: 'RMMPFC',
    9: 'RVMPFC',
    10: 'RPC',
}

ToM_order = [ToM_fROIs[x] for x in range(1, 11)]

ToM = pd.read_csv('../../data/fMRI_ToM/spcorr/spcorr_ToMshort_n149.csv')
ToM.ROI = ToM.ROI.map(ToM_fROIs)
ToM_mean = ToM.groupby('ROI')['Fisher transformed correlation coefficients'].mean()
ToM_sem = ToM.groupby('ROI')['Fisher transformed correlation coefficients'].sem()

ToM_mean = ToM_mean[ToM_order]
ToM_sem = ToM_sem[ToM_order]

print('Language network mean (core LH): %s' % lang_mean.iloc[:5].mean())
print('Language network mean (RH): %s' % lang_mean.iloc[5:10].mean())
print('Language network mean (AngG): %s' % lang_mean.iloc[10:].mean())
print('ToM network mean: %s' % ToM_mean.mean())

mean = pd.concat([lang_mean, ToM_mean], axis=0)
sem = pd.concat([lang_sem, ToM_sem], axis=0)

x_pos_lang = np.concatenate(
    [
        np.arange(len(lang_lh_core)),
        np.arange(len(lang_rh)) + len(lang_lh_core) + 1,
        np.arange(len(lang_angg)) + len(lang_lh_core) + len(lang_rh) + 2,
    ]
)
x_pos_ToM = np.arange(len(ToM_mean)) + len(lang_mean) + 3

ax = plt.gca()

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(top='off', bottom='on', left='off', right='off', labelleft='on', labelbottom='on')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, which='major', axis='y', ls='--', lw=.5, c='k', alpha=.3)
ax.axhline(y=0, lw=1, c='gray', alpha=1)

plt.bar(x_pos_lang, lang_mean, color='r')
plt.errorbar(x_pos_lang, y=lang_mean, yerr=lang_sem, ls='none', color='k')
plt.xticks(x_pos_lang, lang_mean.index, rotation=30, ha='right')

plt.bar(x_pos_ToM, ToM_mean, color='b')
plt.errorbar(x_pos_ToM, y=ToM_mean, yerr=ToM_sem, ls='none', color='k')
plt.xticks(np.concatenate([x_pos_lang, x_pos_ToM]), mean.index, rotation=30, ha='right')

plt.gcf().set_size_inches(8, 2)

plt.tight_layout()

plt.savefig('paper/img/spcorr_lang.png', dpi=300)



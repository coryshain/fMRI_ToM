import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

def get_cond(x):
    if x[:3] == 'bel':
        return 'False Belief'
    return 'False Photo'

def get_cond2(x):
    if x[:3] == 'TOM':
        return 'False Belief'
    return 'Phys Change'

df = pd.read_csv('feats/tom.wsj02to21-gcg15-nol-prtrm-3sm-synproc-+c_+u_+b5000_parsed.5-kenlm.unigram.dlt.lc.itemid-avg.itemmeasures', sep=' ')
df['condition'] = df.itemid.apply(get_cond)

cols = [
    'nwords',
    'nsent',
    'noF',
    'dltcvm',
    'unigramsurp',
    'fwprob5surp',
    'totsurp'
]

name_map = {
    'nwords': 'Num Words',
    'nsent': 'Num Sents',
    'noF': 'Constituent End',
    'wlen': 'Word Length',
    'unigramsurp': 'Unigram Surprisal',
    'fwprob5surp': '5-gram Surprisal',
    'totsurp': 'PCFG Surprisal',
    'dltcvm': 'Integration Cost',
    'dlts': 'Storage Cost',
    'embddepthMin': 'Embedding Depth'
}

df = df[cols + ['condition']]
df = pd.melt(df, id_vars='condition')

fig, axs = plt.subplots(2, 4)

for i, ax in enumerate(axs.flatten()):
    if i >= len(cols):
        ax.axis('off')
        ax.legend(legend.legendHandles, labels, loc='center right', bbox_to_anchor=(1, 0.5))
    else:
        col = cols[i]
        _df = df[df['variable'] == col]
        _df = _df.sort_values('condition', ascending=False)
        cond1, cond2 = 'False Photo', 'False Belief'
        print(col)
        print(cond1, cond2)
        print(ttest_ind(_df[_df['condition'] == cond1]['value'], _df[_df['condition'] == cond2]['value']))
        # ax = sns.violinplot(x='variable', y='value', hue='condition', inner='quartile', data=_df, palette='muted', split=True, ax=ax)
        ax = sns.boxplot(x='variable', y='value', hue='condition', data=_df, palette='binary', ax=ax)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.tick_params(top='off', bottom='off', left='on', right='off', labelleft='on', labelbottom='on')
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylabel(name_map[col])
        ax.get_xaxis().set_visible(False)
        ax.set_xlabel(None)
        if i < len(cols) - 1:
            ax.get_legend().remove()
        else:
            legend = ax.get_legend()
            labels = (x.get_text() for x in legend.get_texts())
            ax.get_legend().remove()
            # ax.legend(legend.legendHandles, labels, loc='center right', bbox_to_anchor=(2., 0.5))


fig.tight_layout()
plt.subplots_adjust(hspace=0.3)
fig.set_size_inches(6.5, 3)
fig.savefig('ling_feature_distribution.png', dpi=300)
plt.close('all')

print()

df = pd.read_csv('feats/tom2deen.wsj02to21-gcg15-nol-prtrm-3sm-synproc-+c_+u_+b5000_parsed.5-kenlm.unigram.dlt.lc.docid-avg.itemmeasures', sep=' ')
df['condition'] = df.docid.apply(get_cond2)
df['nsent'] = 3
df = df[cols + ['condition']]
df = pd.melt(df, id_vars='condition')

fig, axs = plt.subplots(2, 4)

for i, ax in enumerate(axs.flatten()):
    if i >= len(cols):
        ax.axis('off')
        ax.legend(legend.legendHandles, labels, loc='center right', bbox_to_anchor=(1, 0.5))
    else:
        col = cols[i]
        _df = df[df['variable'] == col]
        cond1, cond2 = 'Phys Change', 'False Belief'
        print(col)
        print(ttest_ind(_df[_df['condition'] == cond1]['value'], _df[_df['condition'] == cond2]['value']))
        # ax = sns.violinplot(x='variable', y='value', hue='condition', inner='quartile', data=_df, palette='muted', split=True, ax=ax)
        ax = sns.boxplot(x='variable', y='value', hue='condition', data=_df, palette='binary', ax=ax)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.tick_params(top='off', bottom='off', left='on', right='off', labelleft='on', labelbottom='on')
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylabel(name_map[col])
        ax.get_xaxis().set_visible(False)
        ax.set_xlabel(None)
        if i < len(cols) - 1:
            ax.get_legend().remove()
        else:
            legend = ax.get_legend()
            labels = (x.get_text() for x in legend.get_texts())
            ax.get_legend().remove()
            # ax.legend(legend.legendHandles, labels, loc='center right', bbox_to_anchor=(2., 0.5))

fig.tight_layout()
plt.subplots_adjust(hspace=0.3)
fig.set_size_inches(6.5, 3)
fig.savefig('ling_feature_distribution_tom2.png', dpi=300)
plt.close('all')

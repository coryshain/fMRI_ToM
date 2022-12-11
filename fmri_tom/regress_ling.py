import pickle
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import argparse

def z(x):
    return (x - x.mean()) / x.std()

def var(m, weights):
    w = np.array(weights)
    while len(w.shape) < 2:
        w = w[..., None]
    cov = m.cov_params()

    return np.dot(w.T, np.dot(cov, w))

ling_baselines = [
    ['nwords'],
    ['nsent'],
    ['noF'],
    ['unigramsurp'],
    ['fwprob5surp'],
    ['totsurp'],
    ['dltcvm'],
]
ling_baselines.append([x[0] for x in ling_baselines])

if __name__ == '__main__':
    argparser = argparse.ArgumentParser('''
    Regress linguistic features in ToM experiment
    ''')
    args = argparser.parse_args()

    print('Fitting GLM...')

    b_v_p = {}
    for baseline in ling_baselines:
        baseline_key = '.'.join(baseline)
        models = {}

        for network in ('lang', 'ToM'):
            df = pd.read_csv('../../data/fMRI_ToM/ToMshort_n149/%sfROIs_ToMshort_item_n149_data.csv' % network)
            ling_feats = pd.read_csv('feats/tom.wsj02to21-gcg15-nol-prtrm-3sm-synproc-+c_+u_+b5000_parsed.5-kenlm.unigram.dlt.lc.itemid-avg.itemmeasures', sep=' ')
            df = df.merge(ling_feats, left_on='Effect', right_on='itemid')
            df['condition'] = df.itemid.str[:3]

            for k, v in df.groupby(['Subject', 'ROI']):
                k = k + (network,)
                y = v['EffectSize']
                if len(y[~y.isna()]):
                    X = sm.add_constant(v[baseline])
                    X = X.rename(lambda x: 'Intercept' if x == 'const' else x, axis=1)
                    for pred in baseline:
                        X[pred] = z(X[pred])
                    m = sm.OLS(y, X)
                    result = m.fit()
                    models[k] = result
                    b = result.params.loc[baseline].values[None, ...]
                    delta = (b * X[baseline]).sum(axis=1)
                    y_new = y.values - delta.values
                    v['EffectSizeNew'] = y_new
                    if k not in b_v_p:
                        bel = v[v.condition == 'bel'].EffectSize.mean()
                        pho = v[v.condition == 'pho'].EffectSize.mean()
                        b_v_p[k] = pd.DataFrame({
                            'Subject': k[0],
                            'ROI': k[1],
                            'Network': k[2],
                            'Effect_bel': [bel],
                            'Effect_pho': [pho],
                            'BvP': [bel - pho],
                        })
                    bel_new = v[v.condition == 'bel'].EffectSizeNew.mean()
                    pho_new = v[v.condition == 'pho'].EffectSizeNew.mean()
                    bvp = bel_new - pho_new
                    delta = bvp - b_v_p[k]['BvP']
                    b_v_p[k]['Effect_%s' % baseline_key] = np.sum(b)
                    b_v_p[k]['Effect_bel.%s' % baseline_key] = bel_new
                    b_v_p[k]['Effect_pho.%s' % baseline_key] = pho_new
                    b_v_p[k]['BvP.%s' % baseline_key] = bvp
                    b_v_p[k]['BvPdelta.%s' % baseline_key] = delta
                else:
                    print('Model %s had all NaN response measures. Skipping...\n' % ', '.join(k))

        if not os.path.exists('ling_glm'):
            os.makedirs('ling_glm')

        with open('ling_glm/glm.%s.obj' % baseline_key, 'wb') as f:
            pickle.dump(models, f)

    b_v_p = pd.concat([b_v_p[k] for k in b_v_p], axis=0)
    b_v_p.to_csv('tom_effects_ling_controls.csv', index=False)
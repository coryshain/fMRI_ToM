import sys
import os
import re
import argparse
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

def get_stars(x):
    if x > 0.1:
        return ''
    if x > 0.05:
        return '.'
    if x > 0.01:
        return '*'
    if x > 0.001:
        return '**'
    return '***'

def get_network_fdr(x):
    if x.fROI != 'overall':
        return x.network + 'Ind'
    return x.network

def format_froi(x):
    if x.startswith('roi'):
        return 'roi' + '%02d' % int(x[3:])
    return x

def compute_row(path, stars=True):
    basename_parts = os.path.basename(path).split('.')[0].split('_')
    contrast = basename_parts[0]
    network = basename_parts[1]
    if len(basename_parts) == 3:
        roi = basename_parts[2]
    else:
        roi = 'overall'
    roi = format_froi(roi)

    conv = 't'
    sing = 'f'
    coefname = '---'
    beta = '---'
    se = '---'
    t = '---'
    p = '---'

    with open(path, 'r') as f:
        in_lrt = False
        in_full = False
        in_fixef = False
        for line in f:
            if 'isSingular' in line:
                sing = 't'
            if 'failed to converge' in line:
                conv = 'f'
            if not line.strip():
                if in_full and in_fixef:
                    in_full = False
                    in_fixef = False
            if line.startswith('Full model'):
                in_full = True
            elif line.startswith('Fixed effects:') or line.startswith('Coefficients:'):
                in_fixef = True
            elif line.startswith('Effect') or line.startswith('ismentalTRUE:isverbalTRUE') or (line.startswith('(Intercept)') and '-' in contrast):
                if in_full and in_fixef:
                    coefname, beta, se, t = line.strip().split()[:4]
                    beta = float(beta)
                    se = float(se)
                    t = float(t)
            elif line.strip().startswith('npar') or line.strip().startswith('Res.Df'):
                in_lrt = True
            elif line.startswith('m1') or line.startswith('m_full') or line.startswith('2'):
                if in_lrt:
                    line_parts = line.strip().split()
                    if '*' in line_parts[-1] or line_parts[-1] == '.':
                        line_parts = line_parts[:-1]
                    p = float(line_parts[-1])
                    in_lrt = False
            elif line.startswith('Ablated model'):
                break

    return {
        'network': network,
        'contrast': contrast,
        'fROI': roi,
        'conv': conv,
        'sing': sing,
        'coefname': coefname,
        'beta': beta,
        'se': se,
        't': t,
        'p': p
    }


def correct_p(x):
    p = fdrcorrection(x['p'].values, method='negcorr')[1]
    x['p_fdr'] = p
    return x


if __name__ == '__main__':
    argparser = argparse.ArgumentParser('''
    Get table of significance values for conlen tests.
    ''')
    args = argparser.parse_args()
    cols = ['network', 'contrast', 'fROI', 'conv', 'sing', 'coefname', 'beta', 'se', 't', 'p']
    rows = []
    directory = 'tests'
    paths = sorted([os.path.join(directory, x) for x in os.listdir(directory) if x.endswith('.lrt.summary.txt')])
    rows += [compute_row(path) for path in paths]

    rows = sorted(rows, key=lambda x: (x['network'], x['contrast'], x['fROI']))

    df = pd.DataFrame(rows)
    df['network_fdr'] = df[['network', 'fROI']].apply(get_network_fdr, axis=1)
    df_p = df.groupby(['network_fdr', 'contrast']).apply(correct_p)[['network_fdr', 'contrast', 'fROI', 'p_fdr']]
    df = pd.merge(df, df_p, on=['network_fdr', 'contrast', 'fROI'])
    df['signif'] = df['p_fdr'].apply(get_stars)
    df.to_csv('ToM_signif.csv', index=False)

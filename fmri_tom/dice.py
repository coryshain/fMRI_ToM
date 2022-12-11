import pandas as pd

# LANG fROI vs p < 0.001 unc
df = pd.read_csv('../../data/fMRI_ToM/overlap/langlocSN_mROI_unc0001_n149.csv')
df['0.001-func ratio'] = df.unc0001.astype(float) / df.mROI.astype(float)
df['% func in 0.001'] = df.overlap.astype(float) / df.mROI.astype(float) * 100

print('LANG fROI vs. p unc < 0.001:')
print('Mean size ratio (p < 0.001 to fROI): %s' % df['0.001-func ratio'].mean())
print('Mean DICE: %s' % df['mROI-unc0001'].mean())
print('Mean %% func covered by p < 0.001: %s' % df['% func in 0.001'].mean())
print()

# ToM fROI vs p < 0.001 unc
df = pd.read_csv('../../data/fMRI_ToM/overlap/ToMshort_mROI_unc0001_n149.csv')
df['0.001-func ratio'] = df.unc0001.astype(float) / df.mROI.astype(float)
df['% func in 0.001'] = df.overlap.astype(float) / df.mROI.astype(float) * 100

print('ToM fROI vs. p unc < 0.001:')
print('Mean size ratio (p < 0.001 to fROI): %s' % df['0.001-func ratio'].mean())
print('Mean DICE: %s' % df['mROI-unc0001'].mean())
print('Mean %% func covered by p < 0.001: %s' % df['% func in 0.001'].mean())
print()

# AngG vs TPJ

df = pd.read_csv('../../data/fMRI_ToM/overlap/TPJ-AngG_n149.csv')
print('AngG vs. TPJ LH:')
print('Mean DICE: %s' % df['LTPJ-AngG'].mean())
print('AngG vs. TPJ RH:')
print('Mean DICE: %s' % df['RTPJ-AngG'].mean())
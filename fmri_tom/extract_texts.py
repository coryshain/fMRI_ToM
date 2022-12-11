import pandas as pd

def get_id(row):
    num = row['itemnumall']
    cond = row['condition']
    if cond == 'fp':
        out = 'pho_%02d'
    else:
        out = 'bel_%02d'
    return out % num

item_codes = [
    21, 44, 4, 27, 1, 9, 39, 20, 19, 40, 14, 26, 24, 48, 2, 43, 3, 28, 31, 33
]

df = pd.read_csv('feats/itemmappings_manual.csv')
df = df[df.itemnumall.isin(item_codes)]
df['itemid'] = df[['itemnumall', 'condition']].apply(get_id, axis=1)
df = df[['itemid', 'story']]
df['word'] = df.story.str.split()
del df['story']
df = df.explode('word')
df['sentend'] = df.word.str.contains('.', regex=False)
df['sentid'] = df.sentend.cumsum().shift(1).fillna(0).astype(int)
del df['sentend']
df.to_csv('feats/tom.itemmeasures', index=False, sep=' ')


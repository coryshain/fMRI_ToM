import pandas as pd

if __name__ == '__main__':
    itemmeasures = []
    docid = None
    sentid = -1
    nwords = None
    with open('feats/tom2_deen.txt', 'r') as f:
        for line in f:
            if line.startswith('TOM') or line.startswith('PHYS'):
                if docid is not None:
                    print(docid)
                    print('%d words' % nwords)
                    print()
                docid = line.strip()
                nwords = 0
            elif line:
                sentid += 1
                for sentpos, word in enumerate(line.strip().split()):
                    itemmeasures.append((docid, sentid, sentpos + 1, word))
                    nwords += 1
    itemmeasures = pd.DataFrame(itemmeasures, columns=['docid', 'sentid', 'sentpos', 'word'])

    itemmeasures.to_csv('tom2deen.itemmeasures', sep=' ', index=False, na_rep='NaN')
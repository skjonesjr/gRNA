from pathlib import Path
import multiprocessing, difflib
import pandas


STEM1 = 'UCUAC'
STEM2 = 'GUAGA'


def main(i):
    df = pandas.read_csv(f'sequences_{i:03d}.csv')
    df['dist1'] = df['sequence1'].apply(lambda x: sum(a[0] == '+' for a in difflib.ndiff(STEM1, x.replace('T', 'U'))))
    df['dist2'] = df['sequence2'].apply(lambda x: sum(a[0] == '+' for a in difflib.ndiff(STEM2, x.replace('T', 'U'))))
    df.to_csv(f'sequences_{i:03d}_diffs.csv', index=False)


if __name__ == '__main__':
    n_processes = multiprocessing.cpu_count() - 1
    args = list(set([int(p.stem.split('_')[1]) for p in Path('.').glob('*.csv')]))

    with multiprocessing.Pool(n_processes) as pool:
        pool.map(main, args)

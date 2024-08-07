# %% [markdown]
# Setup
from pathlib import Path
import collections, itertools, subprocess, tempfile

import tqdm
import numpy as np
import pandas
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

from Bio import SeqFeature
from Bio.Blast import NCBIXML
import ViennaRNA

import sklearn.pipeline, sklearn.preprocessing, sklearn.ensemble, sklearn.model_selection, sklearn.cluster, sklearn.svm, sklearn.metrics, sklearn.linear_model


sklearn.set_config(transform_output='pandas')

with open('.env') as f:
    CONFIG = dict([line.strip().split('=', maxsplit=1) for line in f.readlines()])
RAW = Path(CONFIG['RAW'])
INTERIM = Path(CONFIG['INTERIM'])


direct_repeat = 'AACACCGTAATTTCTACTCTTGTAGAT'  # also forward primer binding site
constant_region = 'AGCTTGGCGTAACTAGATCTTG'  # reverse primer binding site


# %% [markdown]
# # Read data
# Sequence:
# direct_repeat + spacer + 6xT + barcode + GA + PAM + target + constant_region

df = (pandas
      .read_excel('https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4104/MediaObjects/41592_2017_BFnmeth4104_MOESM79_ESM.xlsx', header=1)
      .rename(columns={
          "Guide sequence (5' to 3')": 'spacer',
          "Synthetic target sequence  (5' to 3')": 'syn target',  # barcode + GA + PAM + target + constant_region
          "Target sequence (5' to 3')": 'target',
          'Indel frequency (%)': 'indel'
      })
      )
df


# %% [markdown]
# Spacers are equal to targets?

np.all(df['spacer'] == df['target'])


# %% [markdown]
# Syn targets all end with constant_region?

print(df.iloc[0]['syn target'])
print(df.iloc[1]['syn target'])

df['syn target'].apply(lambda x: x.endswith(constant_region)).all()


# %% [markdown]
# Plot secondary structure

def draw(
    structure: str,
    sequence: str | None = None,
    colors: list | None = None,
    node_size: int = 10,
    scale: float = (1, 1),
    shift: tuple[float, float] = (0, 0),
    ax: plt.Axes | None = None
) -> None:
    """
    Draw a structure
    """
    if sequence is None:
        sequence = [''] * len(structure)
    if colors is None:
        colors = ['skyblue'] * len(structure)
    hide_axis = False
    if ax is None:
        ax = plt.subplot(111)
        hide_axis = True

    # Place nucleotides in a graph with pairings obtained from the structure
    G = nx.Graph()
    stack = []
    for i, (nt, char, color) in enumerate(zip(sequence, structure, colors, strict=True)):
        G.add_node(i, label=nt, color=color)
        if char == '(':
            stack.append(i)
        elif char == ')':
            j = stack.pop()
            G.add_edge(j, i)

        if i > 0:
            G.add_edge(i - 1, i)

    # This layout seems to work well for gRNAs
    pos = nx.kamada_kawai_layout(G)
    # but we need to rotate the structures
    flipped_pos = {node: (1 - x, y) for node, (x, y) in pos.items()}
    rotated_pos = {node: (-y, x) for node, (x, y) in flipped_pos.items()}
    scaled_pos = {node: (x * scale[0], y * scale[1]) for node, (x, y) in rotated_pos.items()}
    shifted_pos = {node: (x + shift[0], y + shift[1]) for node, (x, y) in scaled_pos.items()}

    nx.draw_networkx(G, shifted_pos, with_labels=False, node_size=node_size,
                     node_color=[G.nodes[n]['color'] for n in G.nodes],
                     edge_color='lightgray',
                     hide_ticks=False,
                     ax=ax)
    nx.draw_networkx_labels(G, shifted_pos, labels={n: G.nodes[n]['label'] for n in G.nodes}, font_size=12, hide_ticks=False, ax=ax)
    if hide_axis:
        ax.set_axis_off()


idx = 149
seq = SeqFeature.Seq(direct_repeat + df.iloc[idx]['spacer']).transcribe()
structure, mfe = ViennaRNA.fold(str(seq))
colors = ['salmon'] * len(direct_repeat) + ['skyblue'] * len(df.iloc[idx]['spacer'])
draw(structure, sequence=seq, colors=colors, node_size=100)
plt.title(f'MFE = {mfe:.2f}')


# %% [markdown]
# # Plot editing efficiency

df['indel'].hist(bins=50, grid=False)
plt.title('Cas12a gRNA editing efficiency of 1251 targets')
plt.xlabel('Indel frequency (%)')
plt.ylabel('Count')
plt.savefig(ROOT / 'kim_2017_efficiency.png', dpi=300)


# %% [markdown]
# Analyze the impact of having several A's in a row in a spacer

aseqs = list(reversed(['A' * i for i in range(1, 6)]))


def find_a(spacer):
    for aseq in aseqs:
        if spacer.find(aseq) != -1:
            n_As = len(aseq)
            break
    else:
        n_As = 0

    return n_As


df['n_As'] = df['spacer'].apply(find_a)

sns.violinplot(data=df, x='n_As', y='indel', inner=None)
sns.stripplot(data=df, x='n_As', y='indel', color='black', size=3)
sns.pointplot(data=df, x='n_As', y='indel', color='orange', estimator='median')
plt.title("Number of consequent A's vs indel frequency")


# %% [markdown]
# Editing efficiency vs MFE of spacer

df['spacer MFE'] = df['spacer'].apply(lambda x: ViennaRNA.fold(x)[1])
sns.regplot(data=df, x='spacer MFE', y='indel', marker='.', color='black')


# %% [markdown]
# Editing efficiency vs MFE of the whole gRNA

df['gRNA+tail MFE'] = df['spacer'].apply(lambda x: ViennaRNA.fold(direct_repeat + x)[1])
sns.regplot(data=df, x='gRNA+tail MFE', y='indel', marker='.', color='black')


# %% [markdown]
# Editing efficiency vs MFE of the whole gRNA minus the first 8 bases that get cut off

df['gRNA MFE'] = df['spacer'].apply(lambda x: ViennaRNA.fold(direct_repeat[8:] + x)[1])
sns.regplot(data=df, x='gRNA MFE', y='indel', marker='.', color='black')


# %% [markdown]
# Editing efficiency vs GC content

def _get_gc_content(x):
    counts = collections.Counter(x)
    return (counts.get('C', 0) + counts.get('G', 0)) / counts.total()


def _gc_bin(x):
    k = int(100 * x) // 10 * 10
    if k < 40:
        return '<40'
    elif k >= 80:
        return '>80'
    else:
        return f'{k}-{k + 10}'


df['GC content'] = df['spacer'].apply(_get_gc_content)
df['GC bin'] = df['GC content'].apply(_gc_bin)
# plt.scatter(df['GC content'], df['indel'])
order = ['<40', '40-50', '50-60', '60-70', '70-80', '>80']
sns.swarmplot(data=df, x='GC bin', y='indel', size=1, order=order)
sns.pointplot(data=df, x='GC bin', y='indel', linestyles='', order=order)


# %% [markdown]
# Editing efficiency as  a function of number of paired bases

df['structure'] = df['spacer'].apply(lambda x: ViennaRNA.fold(direct_repeat + x)[0])
ct = df['structure'].apply(lambda x: x[len(direct_repeat):].count('(') + x[len(direct_repeat):].count(')'))
# Or count pairings in the entire gRNA
# ct = df['structure'].apply(lambda x: x.count('(') + x.count(')'))

plt.figure(figsize=(12, 6))
ax = sns.stripplot(data=df, x=ct, y='indel', dodge=True, alpha=.25, zorder=1, color='black')
ax = sns.pointplot(
    data=df, x=ct, y='indel', linestyle='none',
    errorbar=None, markers='d', color='red', estimator='median'
)

colors = ['salmon'] * len(direct_repeat) + ['skyblue'] * len(df.loc[0, 'spacer'])
for i, c in enumerate(np.unique(ct)):
    struct = df.loc[ct == c, 'structure'].iloc[0]
    draw(struct, node_size=.5, scale=(.5, 20), shift=(i, 100), colors=colors, ax=ax)

ax.set_xlabel('Number of paired bases in the spacer')
ax.set_ylabel('% indels')
ax.set_yticks(ticks=[0, 20, 40, 60, 80, 100])


# %% [markdown]
# Editing efficiency as  a function of number of paired bases

df['structure'] = df['spacer'].apply(lambda x: ViennaRNA.fold(direct_repeat + x)[0])
ct = df['structure'].apply(lambda x: x[len(direct_repeat):].count('(') + x[len(direct_repeat):].count(')'))
# Or count pairings in the entire gRNA
# ct = df['structure'].apply(lambda x: x.count('(') + x.count(')'))

plt.figure(figsize=(8, 6))
ax = sns.stripplot(data=df, x=ct, y='indel', dodge=True, alpha=.25, zorder=1, color='black')
ax = sns.pointplot(
    data=df, x=ct, y='indel', linestyle='none',
    errorbar=None, markers='d', color='red', estimator='median'
)

struct = df.loc[ct == 0, 'structure'].iloc[0]
draw(struct, node_size=3, scale=(3, 20), shift=(-5, 20), colors=colors, ax=ax)
struct = df.loc[ct == 19, 'structure'].iloc[0]
draw(struct, node_size=3, scale=(3, 20), shift=(20, 0), colors=colors, ax=ax)

ax.set_xlabel('Predicted paired bases in the spacer')
ax.set_ylabel('% indels')
ax.set_yticks(ticks=[0, 20, 40, 60, 80, 100])


# %% [markdown]
# # Predictions
#
# ## Make features

spacer_size = len(df['spacer'].iloc[0])
cols = []
X = []
for k in range(1, 3):
    # all k-mers
    for nts in itertools.product('ATGC', repeat=k):
        kmer = ''.join(nts)
        kmer_sum = []
        for pos in range(0, spacer_size - k + 1):
            out = df['spacer'].apply(lambda x: int(x[pos: pos + k] == kmer))
            kmer_sum.append(out)
            X.append(out)
            cols.append(f'{kmer}{pos:02d}')
        out = pandas.concat(kmer_sum, axis=1).sum(axis=1)
        X.append(out)
        cols.append(kmer)

out = df['spacer'].apply(_get_gc_content)
X.append(out)
cols.append('GC_content')

X.append((out <= 9).astype(int))
cols.append('GC_<=9')
X.append((out > 9).astype(int))
cols.append('GC_>9')

X.append(df['spacer MFE'])
cols.append('spacer MFE')

X = pandas.concat(X, axis=1)
X.columns = cols

X


# %% [markdown]
# ## Train the model

thresh = df['indel'].quantile(.8)  # equals to 50.0
y = df['indel'] >= thresh

cv = sklearn.model_selection.ShuffleSplit(n_splits=10, train_size=.75, random_state=0)
pipe = sklearn.pipeline.Pipeline([
    ('scaler', sklearn.preprocessing.StandardScaler()),
    ('reg', sklearn.linear_model.LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, max_iter=1000, random_state=0))

])

scores = []
for train, test in tqdm.tqdm(cv.split(X, y)):
    pipe.fit(X.loc[train], y[train])
    # feature selection
    sel = [X.columns[i] for i, coef in enumerate(pipe['reg'].coef_[0]) if coef != 0]
    pipe.fit(X.loc[train, sel], y[train])
    score = pipe.score(X.loc[test, sel], y[test])
    scores.append(score)

print(f'Mean score: {np.mean(scores):.2f}')


# %% [markdown]
# ## Try prediction instead of classification

y = df['indel']

cv = sklearn.model_selection.ShuffleSplit(n_splits=10, train_size=.75, random_state=0)
pipe = sklearn.pipeline.Pipeline([
    ('scaler', sklearn.preprocessing.StandardScaler()),
    ('reg', sklearn.linear_model.ElasticNetCV(cv=4))
])

scores = []
for train, test in tqdm.tqdm(cv.split(X, y)):
    pipe.fit(X.loc[train], y[train])
    sel = [X.columns[i] for i, coef in enumerate(pipe['reg'].coef_) if coef != 0]
    pipe.fit(X.loc[train, sel], y[train])
    y_pred = pipe.predict(X.loc[test, sel])
    score = pipe.score(X.loc[test, sel], y[test])
    scores.append(score)


print(f'Mean score: {np.mean(scores):.2f}')  # Mean score = .29

plt.scatter(y[test], y_pred)
plt.xlabel('measured indel')
plt.ylabel('predicted indel')
plt.title('Last fold results')


# %% [markdown]
# # Try RiNAlmo
# ## Get features

# Run on a GPU server and copy the resulting file to your computer:
# python get_rinalmo_features.py

feats = np.load(ROOT / 'rinalmo-giga-v1_features.npy')
feats.shape


# %% [markdown]
# ## Dimensionality reduction
# Check if all dimensions have enough variance

X = feats.mean(axis=1)
std = X.std(axis=0)
plt.plot(np.sort(std))  # Yeah, all have some info


# %% [markdown]
# How many dimensions are sufficient?

pca = sklearn.decomposition.PCA()
pca.fit(X)
plt.plot(range(1, X.shape[0] + 1), pca.explained_variance_ratio_)
plt.xlabel('n_components')
plt.ylabel('explained_variance_ratio')
plt.xlim([0, 20])  # 20 components is totally enough


# %% [markdown]
# Let's use ElasticNet

cv = sklearn.model_selection.ShuffleSplit(n_splits=10, train_size=.75, random_state=0)
pipe = sklearn.pipeline.Pipeline([
    ('scaler', sklearn.preprocessing.StandardScaler()),
    # ('pca', sklearn.decomposition.PCA(n_components=200)),
    ('reg', sklearn.linear_model.ElasticNet())
])

X = feats.mean(axis=1)
y = df['indel']
sklearn.model_selection.cross_val_score(pipe, X, y, cv=cv)


# %% [markdown]
# Let's try to predict as in Kim et al.

X = feats.mean(axis=1)
y = df['indel']

cv = sklearn.model_selection.ShuffleSplit(n_splits=10, train_size=.75, random_state=0)
pipe = sklearn.pipeline.Pipeline([
    ('scaler', sklearn.preprocessing.StandardScaler()),
    ('reg', sklearn.linear_model.ElasticNet())
])

scores = []
for train, test in tqdm.tqdm(cv.split(X, y)):
    pipe.fit(X[train], y[train])
    sel = [coef != 0 for coef in pipe['reg'].coef_]
    pipe.fit(X[train][:, sel], y[train])
    y_pred = pipe.predict(X[test][:, sel])
    score = pipe.score(X[test][:, sel], y[test])
    scores.append(score)


print(f'Mean score: {np.mean(scores):.2f}')  # Mean score = .19

plt.scatter(y[test], y_pred)
plt.xlabel('measured indel')
plt.ylabel('predicted indel')
plt.title('Last fold results')


# %% [markdown]
# # Find Kim et al. genes in the human genome

genome_path = 'human_genome_db' / 'GRCh38_latest_genomic.fna'


# %% [markdown]
# Create human genome database (need to run once only)

subprocess.call(['makeblastdb',
                '-in', RAW / genome_path,
                 '-out', INTERIM / genome_path
                 '-parse_seqids',
                 '-dbtype', 'nucl',
                 ])


# %% [markdown]
# Perform local alignment
# Takes a few minutes

results = []
for target in tqdm.tqdm(df['target'].unique()):
    target_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
    target_file.write(target)
    target_file.close()

    output_file = tempfile.NamedTemporaryFile(delete=False)
    output_file.close()

    cmd = ['blastn',
           '-db', INTERIM / genome_path,
           '-task', 'blastn-short',
           '-query', target_file.name,
           # '-perc_identity', '100',  # would be nice but doesn't do anything
           '-word_size', str(len(target)),
           '-out', output_file.name,
           '-outfmt', '5'
           ]

    # Run BLASTN search
    subprocess.call(cmd)

    records = NCBIXML.parse(open(output_file.name))

    for rec in records:
        for alignment in rec.alignments:
            header = {
                'hit_id': alignment.hit_id,
                'def': alignment.hit_def,
                'accession': alignment.accession,
                'length': alignment.length,
            }
            for hsp in alignment.hsps:
                hit_data = {}
                # Store all fields except the ones that start with _ or frame / strand
                for attr in dir(hsp):
                    if attr not in ['frame', 'strand'] and not attr.startswith('_'):
                        hit_data[attr] = getattr(hsp, attr)
                # Frame info stored here
                hit_data['query_frame'] = hsp.frame[0]
                hit_data['sbjct_frame'] = hsp.frame[1]
            results.append({**header, **hit_data})

df = pandas.DataFrame(results)
df.to_csv(ROOT / 'kim_targets_matches_to_human.csv', index=False)
df


# %% [markdown]
# # Any matches to oncogenes? (answer: none)
# Must run 2-get_oncogenes.py first

onco_targets = pandas.read_csv(ROOT / 'cas12a_targets_in_top200_oncogenes.csv')
for target in df['query']:
    if target[:21] in onco_targets['target']:
        print(target)

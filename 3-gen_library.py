# %%
from pathlib import Path
import random, importlib, difflib

import tqdm
import pandas
import matplotlib.pyplot as plt
import ViennaRNA

from ntgen import grna, generation, insertion, deletion, secondary, visualization
from ntgen.grna import gRNA

importlib.reload(grna)
importlib.reload(generation)
importlib.reload(deletion)
importlib.reload(secondary)
importlib.reload(visualization)

with open('.env') as f:
    CONFIG = dict([line.strip().split('=', maxsplit=1) for line in f.readlines()])
ROOT = Path(CONFIG['ROOT'])

COMPL = str.maketrans('ATUCG', 'UAAGC')

# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # Generate a library of gRNAs


# Wild-type variants
TAIL5 = 'AAUU'
STEM1 = 'UCUAC'
STEM2 = 'GUAGA'
LINKER = 'U'
LOOP = 'UCUU'  # As
LOOPS = {
    'Hk': 'UAUU',  # Adurb193 and Hk
    'Adurb336': 'UGUG',  # Adurb336
    'As': 'UCUU',  # As
    'Fn': 'UGUU',  # Fn, Fn3
    'Lb': 'UAAGU',  # Lb
    'Mb': 'UGUUU',  # Mb
    'Pd': 'UUCG',  # Pd
    'Pi': 'UUGU',  # Pi
}

# From Kim et al., 2017
KIM_TAIL5 = 'AACACCGUAAUU'
# KIM_SCAFFOLD = 'AACACCGTAATTTCTACTCTTGTAGAT'.replace('T', 'U')

# From Jedrzejczyk et al., 2022
HAIRPIN = 'GCCGAAAGGC'
HAIRPIN_STRUCTURE = '(((....)))'


# %% [markdown]
# # Generate all possible stem structures
stem_structures = generation.stem_structures(max_size=16)

# Plot
df_stem_structures = pandas.DataFrame(stem_structures, columns=['stem1', 'stem2'])
df_stem_structures['size'] = df_stem_structures.apply(lambda x: len(x['stem1'] + x['stem2']), axis=1)
agg = df_stem_structures.groupby('size').size()
agg.plot.bar()
plt.ylabel('count of distinct structures')
agg


# %% [markdown]
# Draw stem structural variation samples

uq_sizes = list(set(df_stem_structures['size']))
n_samples = 5
fig, axes = plt.subplots(nrows=n_samples + 1, ncols=len(uq_sizes), figsize=(12, 8))
rng = random.Random(0)
for col, size in enumerate(uq_sizes):
    samples = df_stem_structures[df_stem_structures['size'] == size]
    inds = samples.index.tolist()
    rng.shuffle(inds)

    axes[0][col].text(.5, .5, str(size), ha='center', va='center', fontsize=12)

    i = 0
    for row, idx in enumerate(inds[:n_samples]):
        stem1 = samples.loc[idx, 'stem1']
        stem2 = samples.loc[idx, 'stem2']
        structure = f'....{stem1}...{stem2}.'
        colors = ['skyblue'] * 4 + ['salmon'] * len(stem1) + ['skyblue'] * 3 + ['salmon'] * len(stem2) + ['skyblue']
        visualization.draw(structure, colors=colors, node_size=10, ax=axes[row + 1][col])

    for row in range(n_samples + 1):
        axes[row][col].axis('off')

plt.savefig('structures.png', bbox_inches='tight', pad_inches=0, transparent=False, dpi=300)


# %% [markdown]
# # Generate all possible loops

loops = []
for size in tqdm.trange(1, 12):
    loops += generation.loop_sequences(size)

df_loops = pandas.DataFrame(loops, columns=['loop'])
df_loops['size'] = df_loops['loop'].apply(lambda x: len(x))
agg = df_loops.groupby('size').size()
agg.plot.bar(y='n_norm')
agg


# %% [markdown]
# Generate all possible sequences that match these structures

# This is slow! Run it in parallel on a multicore machine
# python gen_stem_sequences.py


# %% [markdown]
# Unique structures of the wild-type gRNA scaffold deletions

# seqs = deletion.loop_to_loop(STEM1 + 'NNNN' + STEM2)
# structures = set()
# for seq in tqdm.tqdm(seqs):
#     structure, _ = ViennaRNA.fold(seq)
#     structures.add(structure)

# print('Number of unique structures of the wild-type gRNA scaffold deletion:', len(structures))
# structures


# %% [markdown]
# # Select targets
#
# ## Select Kim et al. targets:
# 40 with 0 editing efficiency
# 45 with in-between
# and top-15

df_kim = (pandas
          .read_excel('https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4104/MediaObjects/41592_2017_BFnmeth4104_MOESM79_ESM.xlsx', header=1)
          .rename(columns={
              "Guide sequence (5' to 3')": 'spacer',
              "Synthetic target sequence  (5' to 3')": 'syn target',  # barcode + GA + PAM + target + constant_region
              "Target sequence (5' to 3')": 'target',
              'Indel frequency (%)': 'indel'
          })
          )

df_kim['spacer'] = df_kim['spacer'].apply(lambda x: x.replace('T', 'U'))

kim_worst = df_kim[df_kim['indel'] == 0].head(40)
kim_best = df_kim.sort_values(by='indel').tail(15)

df_kim_medium = df_kim.loc[(df_kim['indel'] > 0) & (~df_kim.index.isin(kim_best.index))]
bin_edges = range(0, 100, 10)
bins = pandas.cut(
    df_kim_medium['indel'],
    bins=bin_edges,
    labels=False,
    include_lowest=True
)
kim_medium = (df_kim_medium
              .groupby(bins)
              .apply(lambda x: x.sample(5, replace=False,random_state=0))
              .reset_index(drop=True)
              )

kim_selection = pandas.concat([
    kim_worst,
    kim_medium,
    kim_best
], keys=['worst', 'medium', 'best'])

kim_selection.loc['medium', 'indel'].hist(bins=bin_edges)


# %% [markdown]
# ## Select oncogenes


# %% [markdown]
# ### Get predicted structures for oncogenes
# (only run once; next time load from a file)

if False:
    df_onco = pandas.read_csv(ROOT / 'cas12a_targets_in_top200_oncogenes.csv')
    df_onco['target'] = df_onco['target'].apply(lambda x: x.replace('T', 'U'))
    # mfe = df_onco['target'].apply(lambda x: pandas.Series(secondary.compute_mfe(x), index=['predicted_structure', 'mfe', 'delta']))

    out = []
    for target in tqdm.tqdm(df_onco['target']):
        out.append(secondary.compute_mfe(target))
    out = pandas.DataFrame(out, columns=['predicted_structure', 'mfe', 'delta'])
    df_onco = pandas.concat([df_onco, out], axis=1)

else:
    df_onco = pandas.read_csv(ROOT / 'cas12a_targets_in_top200_oncogenes_with_mfe.csv')

df_onco


# %% [markdown]
# Select distinct structures from the first exon only, the lowest predicted
# MFE, and the highest delta to the next predicted structure

onco_targets = (df_onco
                .loc[df_onco['exon_no'] == 0]
                .groupby('predicted_structure')
                .apply(lambda x: x.sort_values(by=['mfe', 'delta'], ascending=[True, False])
                       .iloc[0])
                .loc[:, 'target']
                .sample(5, random_state=0)
                .tolist()
                )
len(onco_targets)


# %% [markdown]
# ## Put Kim and onco targets together

kim_truncated = kim_selection.copy()
kim_truncated['spacer'] = kim_truncated['spacer'].apply(lambda x: x[:-2])
targets = {
    'kim_orig': kim_selection['spacer'].tolist(),
    'worst': kim_truncated.loc['worst', 'spacer'].tolist() + onco_targets,
    'best': kim_truncated.loc['best', 'spacer'].tolist() + onco_targets,
    'all': kim_truncated['spacer'].tolist() + onco_targets
}
for key, values in targets.items():
    print(key, len(values))


# %% [markdown]
# # Generate sequences and add them to a dict

gen_seqs = []


# %% [markdown]
# ## 1. Can we reproduce Kim et al. findings?
#
# Generate sequences with Kim et al. targets and their scaffolds

targets_name = 'kim_orig'
for spacer in targets[targets_name]:
    grna = gRNA(tail5=KIM_TAIL5, spacer=spacer, structure='.(.)..')
    gen_seqs.append({
        'kind': 'replication',
        'parts': '',
        'comment': 'kim_scaffold',
        'targets': targets_name,
        'n_targets': len(targets[targets_name]),
        'distance': 0,
        'total_stem_structures': 1,
        'n_stem_structures': 1,
        'total_sequences': 1,
        'n_sequences': 1,
        'grna': grna
    })

grna.draw(node_size=80)


# %% [markdown]
# Generate WT sequences with all targets

targets_name = 'all'
for spacer in targets[targets_name]:
    grna = gRNA(spacer=spacer, structure='.(.)..')
    gen_seqs.append({
        'kind': 'replication',
        'parts': '',
        'comment': 'wild-type_scaffold',
        'targets': targets_name,
        'n_targets': len(targets[targets_name]),
        'distance': 0,
        'total_stem_structures': 1,
        'n_stem_structures': 1,
        'total_sequences': 1,
        'n_sequences': 1,
        'grna': grna
    })

grna.draw(node_size=80)


# %% [markdown]
# ## 2. Can we disrupt good targets?


# %% [markdown]
# Hybridize with itself

combos = [
    (('tail5', ), (TAIL5, ), 0),
    (('tail5', 'stem1'), (TAIL5, STEM1), 0),
    (('tail5', 'stem1', 'loop'), (TAIL5, STEM1, LOOP), 0),
    (('loop', 'stem2'), (LOOP, STEM2), len(TAIL5 + STEM1)),
    (('stem2', ), (STEM2, ), len(TAIL5 + STEM1 + LOOP)),
]
targets_name = 'best'
for combo_name, parts, offset in combos:
    for spacer in targets[targets_name]:
        rev = spacer.translate(COMPL)[::-1]

        new_parts = {}
        new_parts_len = 0
        start = 0
        for part_name, part in zip(combo_name, parts):
            new_parts[part_name] = rev[start: start + len(part)]
            if part_name == 'stem1':
                new_parts['stem2'] = new_parts[part_name].translate(COMPL)[::-1]
            elif part_name == 'stem2':
                new_parts['stem1'] = new_parts[part_name].translate(COMPL)[::-1]
            start += len(part)
            new_parts_len += len(part)

        part1 = '(' * new_parts_len
        part2 = ')' * new_parts_len
        structure = '.' * offset + part1 + '.' * (len(TAIL5 + STEM1 + LOOP + STEM2 + LINKER + spacer) - len(part1 + part2) - offset) + part2
        grna = gRNA(spacer=spacer, structure=structure, **new_parts)
        gen_seqs.append({
            'kind': 'disruption',
            'parts': '+'.join(combo_name),
            'comment': '',
            'targets': targets_name,
            'n_targets': len(targets[targets_name]),
            'distance': -1,
            'total_stem_structures': 1,
            'n_stem_structures': 1,
            'total_sequences': 1,
            'n_sequences': 1,
            'grna': grna
        })

    grna.draw(node_size=80)


# %% [markdown]
# Add poly-A sequence in the scaffold

grna = gRNA(stem1='AAAAA', loop='AAAA', stem2='UUUUU', spacer=targets['best'][0], structure='.(.)...')
gen_seqs.append({
    'kind': 'disruption',
    'parts': 'stem1+loop+stem2',
    'comment': 'poly-A',
    'targets': 'custom',
    'n_targets': 1,
    'distance': -1,
    'total_stem_structures': -1,
    'n_stem_structures': -1,
    'total_sequences': 1,
    'n_sequences': 1,
    'grna': grna
})
grna.draw(node_size=80)


# %% [markdown]
# Add a hairpin in the stem

grna = gRNA(tail5=TAIL5 + STEM1, stem1='AGAUAGGCCGC', loop='UUCG', stem2='GCGGCCUAUCU', linker='', spacer=targets['best'][0], structure='.(.)...')
gen_seqs.append({
    'kind': 'disruption',
    'parts': 'tail5+stem1+loop+stem2+linker',
    'comment': 'terminator',
    'targets': 'custom',
    'n_targets': 1,
    'distance': -1,
    'total_stem_structures': -1,
    'n_stem_structures': -1,
    'total_sequences': 1,
    'n_sequences': 1,
    'grna': grna
})
grna.draw(node_size=80)


# %% [markdown]
# ## 3. Can we fix targets that do not work?
#
# Adding 5' and 3' hairpin

rev_tail5 = TAIL5.translate(COMPL)[::-1]
hairpin5 = rev_tail5[:3] + 'GCG' + TAIL5
hairpin_structure = '(' * 3 + '.' * 4 + ')' * 3
targets_name = 'worst'
for spacer in targets[targets_name]:
    rev_spacer = spacer.translate(COMPL)[::-1]
    hairpin3 = spacer[-2:][::-1] + rev_spacer[2:5]
    structure = hairpin_structure + '(' * len(STEM1) + '.' * len(LOOP) + ')' * len(STEM2) + '.' * (len(LINKER + spacer) - 5) + hairpin_structure
    grna = gRNA(tail5=hairpin5, tail3=hairpin3, spacer=spacer, structure=structure)
    gen_seqs.append({
        'kind': 'fix',
        'parts': 'tail5+tail3',
        'comment': 'hairpin',
        'targets': targets_name,
        'n_targets': len(targets[targets_name]),
        'distance': -1,
        'total_stem_structures': 1,
        'n_stem_structures': 1,
        'total_sequences': 1,
        'n_sequences': 1,
        'grna': grna
    })

grna.draw(node_size=80)


# %% [markdown]
# All stem structures with a random 5' tail and loop

targets_name = 'all'
rng = random.Random(0)
structures = rng.choices(list(stem_structures), k=3)
loops_sel = df_loops.loc[df_loops['size'] <= 8, 'loop']
for struct1, struct2 in tqdm.tqdm(structures):
    idx = sorted(list(stem_structures)).index((struct1, struct2))
    stem_seqs = pandas.read_csv(ROOT / 'grna' / f'sequences_{idx:03d}.csv')
    entry = stem_seqs.iloc[0]
    assert entry['structure1'] == struct1
    assert entry['structure2'] == struct2

    tail = rng.choice(loops_sel)
    loop = rng.choice(loops_sel)

    for spacer in targets['all']:
        structure = '.' * len(tail) + struct1 + '.' * len(loop) + struct2 + '.' * len(LINKER) + '.' * len(spacer)
        grna = gRNA(tail5=tail, stem1=entry['sequence1'], loop=loop, stem2=entry['sequence2'], spacer=spacer, structure=structure)
        gen_seqs.append({
            'kind': 'fix',
            'parts': 'stem',
            'comment': 'structure',
            'targets': targets_name,
            'n_targets': len(targets[targets_name]),
            'distance': -1,
            'total_stem_structures': len(stem_structures),
            'n_stem_structures': len(structures),
            'total_sequences': len(stem_seqs) * len(loops_sel) ** 2,
            'n_sequences': 1,
            'grna': grna
        })

    grna.draw(node_size=80)


# %% [markdown]
# Minimize hybridization with the scaffold

# Load sequences for the wild-type stems

stems_sorted = sorted(list(stem_structures))
idx = stems_sorted.index(('(((((', ')))))'))
stem_seqs = pandas.read_csv(ROOT / 'grna' / f'sequences_{idx:03d}.csv')


# %% [markdown]
# Generate sequences

scaffold = TAIL5 + STEM1 + LOOP + STEM2 + LINKER
scaffold_struct = '.' * len(TAIL5) + '(' * len(STEM1) + '.' * len(LOOP) + ')' * len(STEM2) + '.' * len(LINKER)
s1 = 1  # index of STEM1 in parts
s2 = 3  # index of STEM2 in parts

rng = random.Random(0)
out = []
for spacer in tqdm.tqdm(targets['all']):
    parts = [
        [TAIL5, '.'],
        [STEM1, '('],
        [LOOP, '.'],
        [STEM2, ')'],
        [LINKER, '.']
    ]
    grna = None
    for c in tqdm.trange(200, leave=False):
        scaffold = ''.join([p[0] for p in parts])
        structure, mfe = ViennaRNA.fold(scaffold + spacer)
        s = 0
        stem2_changed = False
        any_changes = False
        linker_only = False
        for i, (part, kind) in enumerate(parts):
            if i == s2 and stem2_changed:
                continue

            expected_structure = kind * len(part)
            pred_structure = structure[s: s + len(part)]
            if pred_structure != expected_structure:
                any_changes = True
                if kind == '.':
                    parts[i][0] = rng.choice(df_loops.loc[df_loops['size'] == len(part), 'loop'].values)
                    if i == 4:
                        linker_only = True
                else:
                    entry = rng.choice(stem_seqs.values)
                    if kind == '(':
                        parts[i][0] = entry[2]
                        parts[s2][0] = entry[3]
                    else:
                        parts[i][0] = entry[3]
                        parts[s1][0] = entry[2]

            s += len(part)

        grna = gRNA(tail5=parts[0][0], stem1=parts[1][0], loop=parts[2][0], stem2=parts[3][0], linker=parts[4][0], spacer=spacer, structure=structure)

        if not any_changes:
            if c == 0:
                grna = None  # don't want sequences that are already ok
            break
    else:
        if not linker_only or c == 0:
            raise ValueError

    if grna is None:
        continue

    out.append({
        'kind': 'fix',
        'parts': 'tail5+stem1+loop+stem2+linker',
        'comment': '',
        'targets': 'custom',
        'n_targets': -1,
        'distance': -1,
        'total_stem_structures': -1,
        'n_stem_structures': 1,
        'total_sequences': -1,
        'n_sequences': 1,
        'grna': grna
    })

for entry in out:
    entry['n_targets'] = len(out)
    gen_seqs.append(entry)

out[-1]['grna'].draw(node_size=80)

# grna = gRNA(spacer=out[-1].spacer)
# structure, mfe = ViennaRNA.fold(grna.sequence)
# grna = gRNA(spacer=out[-1].spacer, structure=structure)
# grna.draw(node_size=80)


# %% [markdown]
# # Understanding what matters in the wild-type gRNA
# ## 1. Substitutions in each component independently


# %% [markdown]
# Stem

n_samples = 2  # samples per condition
targets_name = 'best'
stems_sorted = sorted(list(stem_structures))
inds = [i for i, (s1, s2) in enumerate(stems_sorted) if len(s1 + s2) == len(STEM1 + STEM2)]
for idx in inds:
    stem_seqs = pandas.read_csv(ROOT / 'grna' / f'sequences_{idx:03d}.csv')
    stem_seqs['distance'] = stem_seqs.apply(lambda x: (
        sum(a != b for a, b in zip(STEM1, x['sequence1'].replace('T', 'U'))),
        sum(a != b for a, b in zip(STEM2, x['sequence2'].replace('T', 'U')))
    ), axis=1)

    if stem_seqs.loc[0, 'structure1'] == '(' * len(STEM1) and stem_seqs.loc[0, 'structure2'] == ')' * len(STEM2):
        suffix = 'structure-preserving'
    else:
        suffix = 'structure-independent'

    for dist1, dist2 in sorted(stem_seqs['distance'].unique()):
        if dist1 == 0 and dist2 == 0:
            continue

        sel = stem_seqs[stem_seqs['distance'] == (dist1, dist2)]
        n_samples_real = min(n_samples, len(sel))
        samples = sel.sample(n_samples_real, replace=False, random_state=0)

        for _, entry in samples.iterrows():
            scaffold = '.' * len(TAIL5) + entry['structure1'] + '.' * len(LOOP) + entry['structure2'] + '.' * len(LINKER)

            for spacer in targets[targets_name]:
                structure = scaffold + '.' * len(spacer)
                grna = gRNA(
                    stem1=entry['sequence1'],
                    stem2=entry['sequence2'],
                    spacer=spacer,
                    structure=structure)
                gen_seqs.append({
                    'kind': 'substitution',
                    'parts': 'stem',
                    'comment': f'{suffix} dist={dist1}+{dist2}',
                    'targets': targets_name,
                    'n_targets': len(targets[targets_name]),
                    'distance': dist1 + dist2,
                    'total_stem_structures': len(inds),
                    'n_stem_structures': len(inds),
                    'total_sequences': len(sel),
                    'n_sequences': n_samples_real,
                    'grna': grna
                })

grna.draw(node_size=80)


# %% [markdown]
# Tail, loop, and linker

n_samples = 5
parts = [
    ('tail5', TAIL5),
    ('loop', LOOP),
    ('linker', LINKER)
]
targets_name = 'best'
for part_name, wt_part in parts:
    df_loops_sel = df_loops[df_loops['size'] == len(wt_part)]
    distances = df_loops_sel['loop'].apply(lambda x: sum(a != b for a, b in zip(wt_part, x.replace('T', 'U'))))

    for dist in sorted(distances.unique()):
        if dist == 0:
            continue
        sel = df_loops_sel[distances == dist]
        n_samples_real = min(n_samples, len(sel))
        samples = sel.sample(n_samples_real, replace=False, random_state=0)

        for _, entry in samples.iterrows():
            part = {part_name: entry['loop']}
            for spacer in targets[targets_name]:
                grna = gRNA(spacer=spacer, structure='.(.)..', **part)
                gen_seqs.append({
                    'kind': 'substitution',
                    'parts': part_name,
                    'comment': '',
                    'targets': targets_name,
                    'n_targets': len(targets[targets_name]),
                    'distance': dist,
                    'total_stem_structures': 1,
                    'n_stem_structures': 1,
                    'total_sequences': len(sel),
                    'n_sequences': n_samples_real,
                    'grna': grna
                })

    grna.draw(node_size=80)


# %% [markdown]
# ## 2. Deletions


# %% [markdown]
# Tail, loop, linker

parts = [
    ('tail5', TAIL5),
    ('loop', LOOP),
    ('linker', LINKER)
]
targets_name = 'best'
min_loop_size = 3
for part_name, wt_part in parts:
    dels = deletion.loop_to_loop(wt_part)
    if part_name == 'loop':
        dels_sel = [seq for seq in dels if len(seq) >= min_loop_size]
    else:
        dels_sel = dels

    for seq in sorted(dels_sel):
        part = {part_name: seq}
        for spacer in targets[targets_name]:
            grna = gRNA(spacer=spacer, structure='.(.)..', **part)
            gen_seqs.append({
                'kind': 'deletion',
                'parts': part_name,
                'comment': '',
                'targets': targets_name,
                'n_targets': len(targets[targets_name]),
                'distance': -1,
                'total_stem_structures': 1,
                'n_stem_structures': 1,
                'total_sequences': len(dels_sel),
                'n_sequences': len(dels_sel),
                'grna': grna
            })

    grna.draw(node_size=80)


# %% [markdown]
# Stem (retaining some stem structure)

dels = deletion.stem_to_stem((STEM1, STEM2))
dels_df = pandas.DataFrame([a + b for a, b in dels], columns=['stem1', 'stem2', 'struct1', 'struct2'])
counts = dels_df.groupby(['struct1', 'struct2']).size()
targets_name = 'best'
for stem, struct in list(dels):
    stem1, stem2 = stem
    struct1, struct2 = struct
    for spacer in targets[targets_name]:
        structure = '.' * len(TAIL5) + struct1 + '.' * len(LOOP) + struct2 + '.' * len(LINKER + spacer)
        grna = gRNA(stem1=stem1, stem2=stem2, spacer=spacer, structure=structure)
        gen_seqs.append({
            'kind': 'deletion',
            'parts': 'stem',
            'comment': 'paired',
            'targets': targets_name,
            'n_targets': len(targets[targets_name]),
            'distance': -1,
            'total_stem_structures': len(counts),
            'n_stem_structures': len(counts),
            'total_sequences': counts[(struct1, struct2)],
            'n_sequences': counts[(struct1, struct2)],
            'grna': grna
        })

grna.draw(node_size=80)


# %% [markdown]
# ## 3. Insertions


# %% [markdown]
# Tail, loop, and linker

n_samples = 5
max_insertions = 4
parts = [
    ('tail5', TAIL5),
    ('loop', LOOP),
    ('linker', LINKER)
]
targets_name = 'best'
for part_name, wt_part in tqdm.tqdm(parts):
    for dist in range(1, max_insertions + 1):
        df_loops_sel = df_loops[df_loops['size'] == len(wt_part) + dist]
        distances = df_loops_sel['loop'].apply(lambda x: sum(a[0] == '+' for a in difflib.ndiff(wt_part, x.replace('T', 'U'))))
        sel = df_loops_sel[distances == dist]
        n_samples_real = min(n_samples, len(sel))
        samples = sel.sample(n_samples_real, replace=False, random_state=0)

        for _, entry in samples.iterrows():
            part = {part_name: entry['loop']}
            for spacer in targets[targets_name]:
                grna = gRNA(spacer=spacer, structure='.(.)..',**part)
                gen_seqs.append({
                    'kind': 'insertion',
                    'parts': part_name,
                    'comment': '',
                    'targets': targets_name,
                    'n_targets': len(targets[targets_name]),
                    'distance': dist,
                    'total_stem_structures': 1,
                    'n_stem_structures': 1,
                    'total_sequences': len(sel),
                    'n_sequences': n_samples_real,
                    'grna': grna
                })

grna.draw(node_size=80)


# %% [markdown]
# Stem
# Generate stems

ins = insertion.stem(STEM1, STEM2, max_dist=5)
df_stem_insertions = pandas.DataFrame(ins, columns=['stem1', 'stem2', 'structure1', 'structure2'])
df_stem_insertions['distance'] = df_stem_insertions.apply(lambda x: len(x['stem1'] + x['stem2']), axis=1)
df_stem_insertions


# %% [markdown]
# Add a selection to the final set

structures = df_stem_insertions.apply(lambda x: f'{x["structure1"]}{x["structure2"]}', axis=1)
n_samples = 2
targets_name = 'best'
for dist in tqdm.tqdm(df_stem_insertions['distance'].unique()):
    dist_sel = df_stem_insertions[df_stem_insertions['distance'] == dist]
    structure_sel = structures[df_stem_insertions['distance'] == dist]
    uq_structures = structure_sel.unique()
    for structure in uq_structures:
        sel = dist_sel[structure_sel == structure]
        n_samples_real = min(n_samples, len(sel))
        samples = sel.sample(n_samples_real, replace=False, random_state=0)
        for _, entry in samples.iterrows():
            for spacer in targets[targets_name]:
                structure = '.' * len(TAIL5) + entry['structure1'] + '.' * len(LOOP) + entry['structure2'] + '.' * len(LINKER + spacer)
                grna = gRNA(stem1=entry['stem1'], stem2=entry['stem2'], spacer=spacer, structure=structure)
                gen_seqs.append({
                    'kind': 'insertion',
                    'parts': 'stem1+stem2',
                    'comment': '',
                    'targets': targets_name,
                    'n_targets': len(targets[targets_name]),
                    'distance': dist,
                    'total_stem_structures': len(uq_structures),
                    'n_stem_structures': len(uq_structures),
                    'total_sequences': len(sel),
                    'n_sequences': n_samples_real,
                    'grna': grna
                })

grna.draw(node_size=80)


# %% [markdown]
# # Make the final DataFrame of all constructs

constructs = []
for entry in gen_seqs:
    grna = entry['grna']
    construct = {k: v for k, v in entry.items() if k not in ['distance', 'total_sequences', 'n_sequences', 'grna']}
    construct['stem_structure'] = grna.get_part_structure('stem1') + grna.get_part_structure('stem2')
    construct['total_sequences'] = entry['total_sequences']
    construct['n_sequences'] = entry['n_sequences']
    construct['distance'] = entry['distance']
    construct['structure'] = grna.structure
    for part_name, part in zip(gRNA.part_names, grna.parts):
        construct[part_name] = part
    constructs.append(construct)

constructs = pandas.DataFrame(constructs)
constructs


# %% [markdown]
# Export

counts = constructs.groupby(constructs.loc[:, :'distance'].columns.tolist()).size()
counts.name = 'n_constructs'
counts.to_excel('counts.xlsx')

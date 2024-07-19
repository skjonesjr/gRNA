# %%
from datetime import datetime
from pathlib import Path
import random, importlib, difflib, mmap, io
import urllib.request

import tqdm
import pandas
import matplotlib.pyplot as plt
import ViennaRNA
import fuzzysearch
import fpdf

from Bio.Seq import Seq
import Bio.Restriction, Bio.SeqIO

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
FINAL = Path(CONFIG['FINAL'])

COMPL = str.maketrans('ATUCG', 'UAAGC')
DNA_COMPL = str.maketrans('ATUCG', 'TAAGC')
BACKBONE = str(list(Bio.SeqIO.parse('../data/2.08 - pUC-GW-Kan-HDV-BsaI.fa', 'fasta'))[0].seq).upper()
BACKBONE1 = BACKBONE[:419]
BACKBONE2 = BACKBONE[532:]
SEED_REGION_SIZE = 6
RESTRICTION_ENZYME = Bio.Restriction.BsaI
MAX_LIBRARY_SIZE = 12472 - 50  # 50 sequences left for GenScript's quality checks
PAM = 'TTTA'
PROMOTER_EXT = 'GG'


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
# # Define construct
# ## Define parts


class Construct:

    def __init__(self,
                 backbone1: str,
                 primer5: str,
                 barcode: str,
                 constant: str,
                 pam: str,
                 promoter: str,
                 promoter_ext: str,
                 grna: gRNA,
                 primer3: str,
                 backbone2: str,
                 restriction_enzyme: str | None = None,
                 primer_tol=1,
                 min_off_target_size=6,
                 primer_off_tol=1,
                 ):

        target_rev = grna.spacer.translate(DNA_COMPL)[::-1]
        pam_rev = pam.translate(DNA_COMPL)[::-1]
        self.backbone1 = backbone1
        self.primer5 = primer5
        self.pam = pam
        self.grna = grna
        self.primer3_rev_compl = primer3.translate(DNA_COMPL)[::-1]
        self.backbone2 = backbone2

        parts = [backbone1, self.primer5, barcode, constant, target_rev, pam_rev, promoter, promoter_ext] + grna.parts + [self.primer3_rev_compl, backbone2]
        self.parts = [p.replace('U', 'T') for p in parts]
        self.part_names = ['backbone1', 'primer5', 'barcode', 'constant', 'target_rev_compl', 'pam_rev_compl', 'promoter', 'promoter_ext'] + grna.part_names + ['primer3_rev_compl', 'backbone2']

        self.restriction_enzyme = restriction_enzyme
        if isinstance(restriction_enzyme, str):
            self.restriction_enzyme = getattr(Bio.Restriction, restriction_enzyme)
        self.primer_tol = primer_tol
        self.min_off_target_size = min_off_target_size
        self.primer_off_tol = primer_off_tol

    @property
    def sequence(self):
        return ''.join(self.parts)

    @property
    def is_valid(self):
        return not self.has_re_sites and not self.has_off_targets and not self.has_off_primers

    @property
    def re_sites(self):
        linear = self.backbone1 == '' and self.backbone2 == ''
        return detect_re_sites(self.sequence, linear=linear)

    @property
    def has_re_sites(self):
        return len(self.re_sites) != 2

    @property
    def off_targets(self):
        seed = self.pam + self.grna.spacer[:self.min_off_target_size].replace('U', 'T')
        return self.sequence.count(seed) + self.sequence.translate(DNA_COMPL)[::-1].count(seed)

    @property
    def has_off_targets(self):
        return self.off_targets != 1

    @property
    def off_primers(self):
        rev_compl = self.sequence.translate(DNA_COMPL)
        return [
            self.n_fuzzy_matches(self.primer5, self.sequence) != 1,
            self.n_fuzzy_matches(self.primer5, rev_compl) != 0,
            self.n_fuzzy_matches(self.primer3_rev_compl, self.sequence) != 1,
            self.n_fuzzy_matches(self.primer3_rev_compl, rev_compl) != 0,
        ]

    @property
    def has_off_primers(self):
        return any(self.off_primers)

    def n_fuzzy_matches(self, short_seq, long_seq):
        return len(fuzzysearch.find_near_matches(short_seq, long_seq, max_l_dist=self.primer_off_tol))

    def to_dict(self):
        out = {}
        for part_name, part in zip(self.part_names, self.parts):
            if part_name.startswith('backbone'):
                continue
            out[part_name] = part
        return out


def get_barcodes(size=17, n_errors=2):
    freebarcodes_url = 'https://raw.githubusercontent.com/finkelsteinlab/freebarcodes/master/barcodes/barcodes{}-{}.txt'
    url = freebarcodes_url.format(size, n_errors)
    with urllib.request.urlopen(url) as response:
        freebarcodes = response.read().decode().split()
    return freebarcodes


def detect_re_sites(seq, linear=True):
    seq = Seq('N' * 10 + seq + 'N' * 10)
    return RESTRICTION_ENZYME.search(seq, linear=linear)


def detect_off_targets(seq, spacer):
    seed = PAM + spacer[:SEED_REGION_SIZE].replace('U', 'T')
    return seq.count(seed) + seq.translate(DNA_COMPL)[::-1].count(seed)


def construct_factory(grna, barcode, constant_buffer='', skip_ext=False):
    primer5 = 'TCAATCTGTGGTCTCTCAGG'
    primer3 = 'CGTGTTCTAGGTCTCAGGCC'  # 5' to 3'
    return Construct(
        backbone1=BACKBONE1,
        primer5=primer5,
        barcode=barcode,
        constant='CAGAT' + constant_buffer,
        pam=PAM,
        promoter='TAATACGACTCACTATA',
        promoter_ext=PROMOTER_EXT if not skip_ext else '',
        grna=grna,
        primer3=primer3,
        backbone2=BACKBONE2,
        restriction_enzyme=RESTRICTION_ENZYME,
        min_off_target_size=SEED_REGION_SIZE,
    )


# %% [markdown]
# # Select targets


# %% [markdown]
# ## Invalid seeds
#
# Backbone has some PAM sequences, so all spacers whose seed region
# matches those sequences are invalid

region_size = len(PAM) + SEED_REGION_SIZE
bad_seeds = set()
for fragment in [BACKBONE1, BACKBONE2, BACKBONE2[-region_size:] + BACKBONE1[:region_size]]:
    for frag in [fragment, fragment.translate(DNA_COMPL)[::-1]]:
        idx = -1
        while True:
            idx = frag.find(PAM, idx + 1)
            if idx == -1:
                break
            seed = frag[idx + len(PAM): idx + region_size]
            bad_seeds.add(seed)

print('Bad seeds:', bad_seeds)


# %% [markdown]
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

df_kim = df_kim[
    df_kim['spacer'].apply(lambda x: (len(detect_re_sites(x)) == 0) & (x[:SEED_REGION_SIZE] not in bad_seeds))
]
df_kim['spacer'] = df_kim['spacer'].apply(lambda x: x.replace('T', 'U'))
df_kim


# %% [markdown]
# ## Select oncogenes
#
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
                .loc[df_onco['target'].apply(lambda x: len(detect_re_sites(x.replace('U', 'T'))) == 0)]
                .loc[~df_onco['target'].apply(lambda x: x[:SEED_REGION_SIZE].replace('U', 'T') in bad_seeds)]
                .groupby('predicted_structure', group_keys=True)
                .apply(lambda x: x.sort_values(by=['mfe', 'delta'], ascending=[True, False])
                       .iloc[0])
                .loc[:, 'target']
                )
len(onco_targets)


# %% [markdown]
# Combine targets together

targets_kim_orig = df_kim['spacer']
kim_truncated = df_kim.copy()
kim_truncated['spacer'] = kim_truncated['spacer'].apply(lambda x: x[:-2])

targets_kim = kim_truncated['spacer'].tolist()
targets_onco = onco_targets.tolist()
targets = targets_kim + targets_onco

targets_best20 = kim_truncated.sort_values(by='indel').tail(20)['spacer']

worst_sel = kim_truncated['indel'] == 0
targets_worst10 = kim_truncated[worst_sel].sample(10, random_state=0)['spacer'].tolist()

kim_medium = kim_truncated[(~worst_sel) & (~kim_truncated.index.isin(targets_best20.index))]  # .sample(10, random_state=0)['spacer'].tolist()
bin_edges = range(0, 100, 10)
bins = pandas.cut(
    kim_medium['indel'],
    bins=bin_edges,
    labels=False,
    include_lowest=True
)
targets_medium10 = (kim_medium
                    .groupby(bins)
                    .apply(lambda x: x.sample(10, replace=False,random_state=0))
                    .reset_index(drop=True)
                    .loc[:, 'spacer']
                    .tolist()
                    )
targets_best20 = targets_best20.tolist()
targets_best5 = targets_best20[-5:]

targets_mix = targets_best20 + targets_medium10 + targets_worst10 + onco_targets.sample(10, random_state=0).tolist()

print('Number of targets:', len(targets))
print('Number of targets_mix:', len(targets_mix))


# %% [markdown]
# # Select barcodes
# Barcodes must not overlap with targets

test_barcodes = get_barcodes(size=17, n_errors=2)
barcodes = []
for barcode in tqdm.tqdm(test_barcodes):
    no_off_target = [t[:SEED_REGION_SIZE] not in barcode for t in targets]
    re_sites = detect_re_sites(barcode)
    if all(no_off_target) and len(re_sites) == 0:
        barcodes.append(barcode)
print('Number of all barcodes:', len(test_barcodes))
print('Number of valid barcodes:', len(barcodes))


# %% [markdown]
# # Generate sequences and add them to a dict

gen_seqs = {}


# %% [markdown]
# ## 1. Can we reproduce Kim et al. findings?
#
# Generate sequences with Kim et al. targets and their scaffolds

gen_seqs[0] = []
# for spacer in targets[targets_name]:
for spacer in targets_kim_orig:
    grna = gRNA(tail5=KIM_TAIL5, spacer=spacer, structure='.(.)..')
    gen_seqs[0].append({
        'kind': 'replication',
        'parts': '',
        'comment': 'kim_scaffold',
        'targets': 'kim_orig',
        'distance': 0,
        'total_stem_structures': 1,
        'n_stem_structures': 1,
        'total_sequences': 1,
        'n_sequences': 1,
        'grna': grna
    })

grna.draw(node_size=80)


# %% [markdown]
# Generate WT sequences with all targets, including all oncogenes

gen_seqs[1] = []
for spacer in targets:
    grna = gRNA(spacer=spacer, structure='.(.)..')
    comment = 'oncogene' if spacer in onco_targets else 'kim'
    gen_seqs[1].append({
        'kind': 'replication',
        'parts': '',
        'comment': 'wild-type_scaffold',
        'targets': 'kim_mix100' if spacer in targets_kim else 'onco',
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
gen_seqs[2] = []
for combo_name, parts, offset in combos:
    for spacer in targets:
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
        gen_seqs[2].append({
            'kind': 'disruption',
            'parts': '+'.join(combo_name),
            'comment': '',
            'targets': 'kim_best' if spacer in targets_kim else 'onco5',
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

gen_seqs[3] = []
for spacer in targets_best20:
    grna = gRNA(stem1='AAAAA', loop='AAAA', stem2='UUUUU', spacer=spacer, structure='.(.)...')
    gen_seqs[3].append({
        'kind': 'disruption',
        'parts': 'stem1+loop+stem2',
        'comment': 'poly-A',
        'targets': 'kim_best5',
        'distance': -1,
        'total_stem_structures': -1,
        'n_stem_structures': -1,
        'total_sequences': 1,
        'n_sequences': 1,
        'grna': grna
    })

grna.draw(node_size=80)


# %% [markdown]
# Add poly-U after the stem, replacing part of spacer

poly_u_size = 10
gen_seqs[4] = []
for spacer in targets_best5:
    spacer_with_u = poly_u_size * 'U' + spacer[poly_u_size:]
    grna = gRNA(spacer=spacer_with_u, structure='.(.)...')
    gen_seqs[4].append({
        'kind': 'disruption',
        'parts': 'spacer',
        'comment': 'poly-U spacer',
        'targets': 'best5_custom',
        'distance': -1,
        'total_stem_structures': 1,
        'n_stem_structures': 1,
        'total_sequences': 1,
        'n_sequences': 1,
        'grna': grna
    })

grna.draw(node_size=80)


# %% [markdown]
# Omit GG after promoter.
# Since this change in not withing gRNA, we'll have to introduce this change
# when producing the whole construct

gen_seqs[5] = []
for spacer in targets_best20:
    grna = gRNA(structure='.(.)...', spacer=spacer)
    gen_seqs[5].append({
        'kind': 'disruption',
        'parts': 'promoter_ext',
        'comment': 'no GG',
        'targets': 'kim_best5' if spacer in targets_kim else 'onco5',
        'distance': -1,
        'total_stem_structures': 1,
        'n_stem_structures': 1,
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
gen_seqs[6] = []
for spacer in targets:
    rev_spacer = spacer.translate(COMPL)[::-1]
    hairpin3 = spacer[-2:][::-1] + rev_spacer[2:5]
    structure = hairpin_structure + '(' * len(STEM1) + '.' * len(LOOP) + ')' * len(STEM2) + '.' * (len(LINKER + spacer) - 5) + hairpin_structure
    grna = gRNA(tail5=hairpin5, tail3=hairpin3, spacer=spacer, structure=structure)
    gen_seqs[6].append({
        'kind': 'fix',
        'parts': 'tail5+tail3',
        'comment': 'hairpin',
        'targets': 'kim_worst' if spacer in targets_kim else 'onco5',
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

rng = random.Random(0)
structures = stem_structures  # rng.sample(list(stem_structures), 3)

# If already computed before, skip to save time
if 'n_seqs' not in globals():
    n_seqs = {}
linker = 'C'
resite1 = RESTRICTION_ENZYME.site
resite2 = RESTRICTION_ENZYME.site.translate(DNA_COMPL)[::-1]
loops_sel = df_loops.loc[df_loops['size'] <= 8, 'loop']
gen_seqs[7] = []
for struct1, struct2 in tqdm.tqdm(structures):
    if struct1 == '' and struct2 == '':
        continue
    idx = sorted(list(stem_structures)).index((struct1, struct2))
    parts = []
    with open(ROOT / 'interim' / f'sequences_{idx:03d}.csv') as f:
        header = f.readline().strip().split(',')
        if idx not in n_seqs:
            buf = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            n_seqs[idx] = buf.read().count(b'\n') - 1  # last line is empty
        # Loop over entries to find the one with no RE sites
        while True:
            data = f.readline().strip().split(',')
            entry = {k:v for k,v in zip(header, data)}
            assert entry['structure1'] == struct1
            assert entry['structure2'] == struct2
            tail = rng.choice(loops_sel)
            loop = rng.choice(loops_sel)

            grna = gRNA(tail5=tail, stem1=entry['sequence1'], loop=loop, stem2=entry['sequence2'], linker=linker)
            seq = grna.sequence.replace('U', 'T')
            re_sites = detect_re_sites('GG' + seq)
            re_overlap1 = any(seq.endswith(resite1[:i]) for i in range(4, len(resite1)))
            re_overlap2 = any(seq.endswith(resite2[:i]) for i in range(4, len(resite2)))
            if len(re_sites) == 0 and not re_overlap1 and not re_overlap2:
                parts.append((tail, loop, entry))
            if len(parts) == 2:
                break

    for tail, loop, entry in parts:
        for spacer in targets_best20:
            structure = '.' * len(tail) + struct1 + '.' * len(loop) + struct2 + '.' * len(linker) + '.' * len(spacer)
            grna = gRNA(tail5=tail, stem1=entry['sequence1'], loop=loop,
                        stem2=entry['sequence2'], linker=linker,
                        spacer=spacer, structure=structure)

            gen_seqs[7].append({
                'kind': 'fix',
                'parts': 'stem',
                'comment': 'structure',
                'targets': 'kim_best3',
                'distance': -1,
                'total_stem_structures': len(stem_structures),
                'n_stem_structures': len(structures),
                'total_sequences': n_seqs[idx] * len(loops_sel) ** 2,
                'n_sequences': 1,
                'grna': grna
            })

grna.draw(node_size=80)


# %% [markdown]
# Minimize hybridization with the scaffold

# Load sequences for the wild-type stems

stems_sorted = sorted(list(stem_structures))
idx = stems_sorted.index(('(((((', ')))))'))
stem_seqs = pandas.read_csv(ROOT / 'interim' / f'sequences_{idx:03d}.csv')


# %% [markdown]
# Generate sequences

scaffold = TAIL5 + STEM1 + LOOP + STEM2 + LINKER
scaffold_struct = '.' * len(TAIL5) + '(' * len(STEM1) + '.' * len(LOOP) + ')' * len(STEM2) + '.' * len(LINKER)
s1 = 1  # index of STEM1 in parts
s2 = 3  # index of STEM2 in parts
max_iters = 1000

rng = random.Random(0)
gen_seqs[8] = []
for spacer in tqdm.tqdm(targets):
    parts = [
        [TAIL5, '.'],
        [STEM1, '('],
        [LOOP, '.'],
        [STEM2, ')'],
        [LINKER, '.']
    ]
    grna = None
    for c in range(max_iters):
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

    gen_seqs[8].append({
        'kind': 'fix',
        'parts': 'tail5+stem1+loop+stem2+linker',
        'comment': '',
        'targets': 'kim_mix100_sel' if spacer in targets_kim else 'onco5',
        'distance': -1,
        'total_stem_structures': -1,
        'n_stem_structures': 1,
        'total_sequences': -1,
        'n_sequences': 1,
        'grna': grna
    })

# Visualize the "bad" structure before fixing
grna = gRNA(spacer=out[-1]['grna'].spacer)
structure, mfe = ViennaRNA.fold(grna.sequence)
grna = gRNA(spacer=out[-1]['grna'].spacer, structure=structure)
grna.draw(node_size=80)

# Visualize the fixed structure
out[-1]['grna'].draw(node_size=80)


# %% [markdown]
# # Understanding what matters in the wild-type gRNA
# ## 1. Substitutions in each component independently


# %% [markdown]
# Stem

n_samples = 2  # samples per condition
stems_sorted = sorted(list(stem_structures))
inds = [i for i, (s1, s2) in enumerate(stems_sorted) if len(s1 + s2) == len(STEM1 + STEM2)]

gen_seqs[9] = []
for idx in inds:
    stem_seqs = pandas.read_csv(ROOT / 'interim' / f'sequences_{idx:03d}.csv')
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
        samples = sel.sample(frac=1, random_state=0)
        # samples = sel.sample(n_samples_real, replace=False, random_state=0)
        samples_good = []
        for _, entry in samples.iterrows():
            for spacer in targets_mix:
                grna = gRNA(
                    stem1=entry['sequence1'],
                    stem2=entry['sequence2'],
                    spacer=spacer
                )
                re_sites = detect_re_sites(grna.sequence.replace('U', 'T'))
                n_off_targets = detect_off_targets(grna.sequence.replace('U', 'T'), grna.spacer.replace('U', 'U'))
                if len(re_sites) != 0 or n_off_targets != 0:
                    break
            else:
                samples_good.append(entry)
                if len(samples_good) == n_samples_real:
                    break
        else:  # all are bad, so take the first two
            samples_good = samples.iloc[:n_samples_real].to_dict(orient='records')

        for entry in samples_good:
            scaffold = '.' * len(TAIL5) + entry['structure1'] + '.' * len(LOOP) + entry['structure2'] + '.' * len(LINKER)

            for spacer in targets_mix:
                structure = scaffold + '.' * len(spacer)
                grna = gRNA(
                    stem1=entry['sequence1'],
                    stem2=entry['sequence2'],
                    spacer=spacer,
                    structure=structure)
                gen_seqs[9].append({
                    'kind': 'substitution',
                    'parts': 'stem',
                    'comment': f'{suffix} dist={dist1}+{dist2}',
                    'targets': 'kim_mix23' if spacer in targets_kim else 'onco5',
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
gen_seqs[10] = []
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
            for spacer in targets_mix:
                grna = gRNA(spacer=spacer, structure='.(.)..', **part)
                gen_seqs[10].append({
                    'kind': 'substitution',
                    'parts': part_name,
                    'comment': '',
                    'targets': 'kim_mix23' if spacer in targets_kim else 'onco5',
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
min_loop_size = 3
gen_seqs[11] = []
for part_name, wt_part in parts:
    dels = deletion.loop_to_loop(wt_part)
    if part_name == 'loop':
        dels = [seq for seq in dels if len(seq) >= min_loop_size]

    if len(dels) == 8:
        dels_sel = rng.sample(list(dels), 6)
    else:
        dels_sel = dels

    for seq in sorted(dels_sel):
        part = {part_name: seq}
        for spacer in targets:
            grna = gRNA(spacer=spacer, structure='.(.)..', **part)
            gen_seqs[11].append({
                'kind': 'deletion',
                'parts': part_name,
                'comment': '',
                'targets': 'kim_best' if spacer in targets_kim else 'onco5',
                'distance': -1,
                'total_stem_structures': 1,
                'n_stem_structures': 1,
                'total_sequences': len(dels),
                'n_sequences': len(dels_sel),
                'grna': grna
            })

    grna.draw(node_size=80)


# %% [markdown]
# Stem (retaining some stem structure)

dels = deletion.stem_to_stem((STEM1, STEM2))
dels_df = pandas.DataFrame([a + b for a, b in dels], columns=['stem1', 'stem2', 'struct1', 'struct2'])
counts = dels_df.groupby(['struct1', 'struct2']).size()
# Drop entries if there are more that 5 for that condition
dels_df = dels_df.groupby(['struct1', 'struct2']).apply(lambda x: x.sample(n=min(len(x), 6), random_state=0)).reset_index(drop=True)

gen_seqs[12] = []
for _, (stem1, stem2, struct1, struct2) in tqdm.tqdm(dels_df.iterrows()):
    for spacer in targets:
        structure = '.' * len(TAIL5) + struct1 + '.' * len(LOOP) + struct2 + '.' * len(LINKER + spacer)
        grna = gRNA(stem1=stem1, stem2=stem2, spacer=spacer, structure=structure)
        gen_seqs[12].append({
            'kind': 'deletion',
            'parts': 'stem',
            'comment': 'paired',
            'targets': 'kim_best' if spacer in targets_kim else 'onco5',
            'distance': -1,
            'total_stem_structures': len(counts),
            'n_stem_structures': len(counts),
            'total_sequences': counts[(struct1, struct2)],
            'n_sequences': len(dels_df[(dels_df['struct1'] == struct1) & (dels_df['struct2'] == struct2)]),
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
gen_seqs[13] = []
for part_name, wt_part in tqdm.tqdm(parts):
    for dist in range(1, max_insertions + 1):
        df_loops_sel = df_loops[df_loops['size'] == len(wt_part) + dist]
        distances = df_loops_sel['loop'].apply(lambda x: sum(a[0] == '+' for a in difflib.ndiff(wt_part, x.replace('T', 'U'))))
        sel = df_loops_sel[distances == dist]
        n_samples_real = min(n_samples, len(sel))
        samples = sel.sample(n_samples_real, replace=False, random_state=0)

        for _, entry in samples.iterrows():
            part = {part_name: entry['loop']}
            for spacer in targets_mix:
                grna = gRNA(spacer=spacer, structure='.(.)..', **part)
                gen_seqs[13].append({
                    'kind': 'insertion',
                    'parts': part_name,
                    'comment': '',
                    'targets': 'kim_mix23' if spacer in targets_kim else 'onco5',
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

if False:
    ins = insertion.stem(STEM1, STEM2, max_dist=5)
    df_stem_insertions = pandas.DataFrame(ins, columns=['stem1', 'stem2', 'structure1', 'structure2'])
    df_stem_insertions['distance'] = df_stem_insertions.apply(lambda x: len(x['stem1'] + x['stem2']), axis=1)
    df_stem_insertions.to_csv(ROOT / 'stem_insertions.csv', index=False)
else:
    df_stem_insertions = pandas.read_csv(ROOT / 'stem_insertions.csv')

df_stem_insertions


# %% [markdown]
# Add a selection to the final set

structures = df_stem_insertions.apply(lambda x: f'{x["structure1"]}{x["structure2"]}', axis=1)
n_samples = 2
gen_seqs[14] = []
for dist in tqdm.tqdm(df_stem_insertions['distance'].unique()):
    dist_sel = df_stem_insertions[df_stem_insertions['distance'] == dist]
    structure_sel = structures[df_stem_insertions['distance'] == dist]
    uq_structures = structure_sel.unique()
    for structure in uq_structures:
        sel = dist_sel[structure_sel == structure]
        n_samples_real = min(n_samples, len(sel))
        samples = sel.sample(n_samples_real, replace=False, random_state=0)
        for _, entry in samples.iterrows():
            for spacer in targets_mix:
                structure = '.' * len(TAIL5) + entry['structure1'] + '.' * len(LOOP) + entry['structure2'] + '.' * len(LINKER + spacer)
                grna = gRNA(stem1=entry['stem1'], stem2=entry['stem2'], spacer=spacer, structure=structure)
                gen_seqs[14].append({
                    'kind': 'insertion',
                    'parts': 'stem1+stem2',
                    'comment': '',
                    'targets': 'kim_mix23' if spacer in targets_kim else 'onco5',
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
#
# ## Identify bad constructs

bad_targets = set()
max_grna_size = 0
bad_constructs = []
for entries in tqdm.tqdm(list(gen_seqs.values())):
    for entry in entries:
        grna = entry['grna']
        max_grna_size = max(max_grna_size, len(grna.sequence))
        construct = construct_factory(
            grna,
            barcode='N' * 20,
            skip_ext=entry['comment'] == 'no GG',
        )
        if not construct.is_valid:
            bad_targets.add(grna.spacer)
            bad_constructs.append(construct)

print('Number of bad targets:', len(bad_targets))


# %% [markdown]
# ## Create the final list of targets

kim_sel = df_kim[(~df_kim['spacer'].isin(bad_targets)) & (~kim_truncated['spacer'].isin(bad_targets))]

kim_zero = kim_sel[kim_sel['indel'] == 0]
sel_worst = kim_sel['spacer'].apply(lambda x: x[:-2] in targets_worst10)
kim_worst = kim_sel[sel_worst]
kim_worst = pandas.concat([kim_worst, kim_zero[~sel_worst].sample(30, random_state=0)])

kim_best = kim_sel.sort_values(by='indel').tail(15)
assert kim_best['spacer'].apply(lambda x: x[:-2] in targets_best20).all()

kim_sel_medium = kim_sel[kim_sel['spacer'].apply(lambda x: x[:-2] in targets_medium10)]

bin_edges = range(0, 100, 10)
bins = pandas.cut(
    kim_sel_medium['indel'],
    bins=bin_edges,
    labels=False,
    include_lowest=True
)
kim_medium = (kim_sel_medium
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
# ## Put Kim and onco targets together

onco_sel = onco_targets.sample(5, random_state=0).tolist()
kim_sel_trunc = kim_selection.copy()
kim_sel_trunc['spacer'] = kim_sel_trunc['spacer'].apply(lambda x: x[:-2])
targets_final = {
    'kim_orig': kim_selection['spacer'].tolist(),
    'kim_worst': kim_sel_trunc.loc['worst', 'spacer'].tolist() + onco_sel,
    'kim_best': kim_sel_trunc.loc['best', 'spacer'].tolist() + onco_sel,
    'kim_best5': kim_sel_trunc.loc['best', 'spacer'].tolist()[-5:],  # the very best
    'kim_best3': kim_sel_trunc.loc['best', 'spacer'].tolist()[-3:],  # the very best
    'kim_mix23': kim_sel_trunc.loc['best', 'spacer'].tolist() + kim_sel_trunc.loc['medium', 'spacer'].tolist()[:1] + kim_sel_trunc.loc['worst', 'spacer'].tolist()[:2] + onco_sel,
    'kim_mix100': kim_sel_trunc['spacer'].tolist(),
    # 'kim_mix100_sel': kim_sel_trunc['spacer'].tolist(),
    'onco5': onco_sel,
    'onco': onco_targets.tolist()
}
for key, values in targets_final.items():
    print(key, len(values))


# %% [markdown]
# ## Compile the final list as a DataFrame

available_barcodes = list(barcodes)  # copy
constructs_list = []
bad_targets = set()
for entries in tqdm.tqdm(gen_seqs.values()):
    for entry in tqdm.tqdm(entries, leave=False):
        if entry['kind'] == 'replication' and entry['comment'] == 'wild-type_scaffold' and entry['targets'] == 'kim_mix100':
            repeats = 2  # two repeats in case one barcode fails
        else:
            repeats = 1

        grna = entry['grna']
        if entry['targets'] in targets_final:
            targets_sel = targets_final[entry['targets']]
            if grna.spacer not in targets_sel:
                continue
        elif entry['targets'] == 'best5_custom':
            targets_sel = targets_final['kim_best5']
            if grna.spacer in bad_targets:
                continue
        elif entry['targets'] == 'kim_mix100_sel':
            targets_sel = targets_final['kim_mix100']
            if grna.spacer not in targets_sel:
                continue
        else:
            raise ValueError('No such target:', entry['targets'])

        if entry['targets'].startswith('kim'):
            spacer = grna.spacer[:len(targets[0])]
            indel = kim_sel_trunc.loc[kim_sel_trunc['spacer'] == spacer, 'indel']
            assert len(indel) == 1
            indel = indel.iloc[0]
        else:
            indel = ''

        # Search for a barcode that keeps construct valid
        for _ in range(repeats):
            for i, barcode in enumerate(available_barcodes):
                construct = construct_factory(
                    grna,
                    barcode=barcode,
                    skip_ext=entry['comment'] == 'no GG',
                    constant_buffer=''.join(rng.choices('ATCG', k=max_grna_size - len(grna.sequence)))
                )
                if construct.is_valid:
                    available_barcodes.pop(i)
                    break
            else:
                raise ValueError

            # targets_prefix = '+onco5' if entry['targets'] not in ['kim_orig', 'best5', 'best3', 'all+oncogenes'] else ''

            construct_dict = {
                'kind': entry['kind'],
                'parts': entry['parts'],
                'comment': entry['comment'],
                'targets': entry['targets'],
                'n_targets': len(targets_sel),
                'kim_indel': indel,
                'total_stem_structures': entry['total_stem_structures'],
                'n_stem_structures': entry['n_stem_structures'],
                'stem_structure': grna.get_part_structure('stem1') + grna.get_part_structure('stem2'),
                'total_sequences': entry['total_sequences'],
                'n_sequences': entry['n_sequences'],
                'distance': entry['distance'],
                'structure': grna.structure,
            }
            construct_dict.update(construct.to_dict())
            constructs_list.append(construct_dict)

constructs = pandas.DataFrame(constructs_list)
constructs.to_csv(FINAL / f"grnas_{datetime.now().strftime('%Y%m%d')}.csv", index=False)
constructs


# %% [markdown]
# ## Export overview of the library

order = ['replication', 'disruption', 'fix', 'substitution', 'deletion', 'insertion']
row_size = 3
count = 0
data = []
for kind in order:
    if kind == 'substitution':
        columns = ['parts']
    else:
        columns = ['parts', 'comment']

    groups = constructs[constructs['kind'] == kind].groupby(columns).groups

    for cols, inds in groups.items():
        if kind == 'substitution' and parts == 'stem1+stem2':
            sel = constructs[(constructs['kind'] == 'substitution') & (constructs['parts'] == 'stem1+stem2')]
            n_stem_structures = sel['n_stem_structures'].unique().sum()
        else:
            n_stem_structures = constructs.loc[inds[0], 'n_stem_structures']

        if len(cols) == 2:
            parts, comment = cols
        else:
            parts = cols
            comment = ''

        construct = constructs.loc[inds[0]]
        grna = gRNA(tail5=construct['tail5'],
                    stem1=construct['stem1'],
                    loop=construct['loop'],
                    stem2=construct['stem2'],
                    linker=construct['linker'],
                    spacer=construct['spacer'],
                    tail3=construct['tail3'],
                    structure=construct['structure']
                    )

        plt.ioff()  # do not display the plot
        plt.figure()
        _ = grna.draw(node_size=80)
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight', dpi=200)
        buf.seek(0)

        if count == 0:
            tmp = []

        tmp.append(((kind, parts, comment, n_stem_structures, len(inds)), buf))

        count += 1

        if count == row_size:
            data.append(tmp)
            count = 0


# %% [markdown]
# Export to pdf

header = '\n'.join(['Kind', 'Parts', 'Comment', 'Structures', 'gRNAs'])
col_widths = [20, 30, 50, 30, 50, 30, 50]

pdf = fpdf.FPDF(orientation="landscape")
pdf.add_page()
with pdf.table(first_row_as_headings=False, col_widths=col_widths) as table:
    for data_row in data:
        row = table.row()

        pdf.set_font('helvetica', 'B', size=10)
        row.cell(header)

        pdf.set_font('helvetica', size=10)
        for annot, img in data_row:
            text = list(annot)  # copy
            if len(text[1]) > 16:
                text[1] = text[1][:11] + '...'
            row.cell('\n'.join([str(t) if t != -1 else '' for t in text]))
            row.cell(img=img)

idx = 5000
construct = constructs.iloc[idx]
pdf.add_page()

pdf.set_font('helvetica', 'B', size=30)
pdf.write(text=f'Random library member (#{idx})')
pdf.ln(h=15)

pdf.set_font('helvetica', size=20)
descr = ', '.join([c for c in construct.loc['kind': 'comment'] if c != ''])
pdf.cell(text=f'Description: {descr}')
pdf.ln(h=10)

pdf.set_x(80)
pdf.cell(text="Sequence 5' to 3'")
pdf.ln(h=10)

resite1 = RESTRICTION_ENZYME.site
resite2 = RESTRICTION_ENZYME.site.translate(DNA_COMPL)[::-1]
for name, part in construct.loc['primer5':'primer3_rev_compl'].items():
    pdf.set_font('helvetica', 'B', size=20)
    pdf.cell(w=70, text=name)
    if name == 'primer5':
        i = part.find(resite1)
        pdf.set_font('courier', size=20)
        pdf.write(text=part[:i])
        pdf.set_font('courier', 'B', size=20)
        pdf.write(text=part[i: i + len(resite1)])
        pdf.set_font('courier', size=20)
        pdf.write(text=part[i + len(resite1):])
    elif name == 'primer3_rev_compl':
        i = part.find(resite2)
        pdf.set_font('courier', size=20)
        pdf.write(text=part[:i])
        pdf.set_font('courier', 'B', size=20)
        pdf.write(text=part[i: i + len(resite2)])
        pdf.set_font('courier', size=20)
        pdf.write(text=part[i + len(resite2):])
    else:
        pdf.set_font('courier', size=20)
        pdf.write(text=part)
    pdf.ln()

pdf.ln(h=10)

pdf.set_font('helvetica', size=20)
pdf.cell(text=f"{RESTRICTION_ENZYME.__name__} recognition site highlighted")

pdf.output(FINAL / f"overview_{datetime.now().strftime('%Y%m%d')}.pdf")


# %% [markdown]
# ## Export counts

counts = constructs.groupby(constructs.loc[:, :'distance'].columns.tolist()).size()
counts.name = 'n_constructs'
counts.to_excel('counts.xlsx')



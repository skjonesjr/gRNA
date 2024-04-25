# %% [markdown]
# # Generate a library of gRNAs

import typing, itertools, random, collections

import tqdm
import matplotlib.pyplot as plt
import pandas
import networkx as nx

import ViennaRNA


# Wild-type variants
TAIL5 = 'AAUU'
STEM1 = 'UCUAC'
STEM2 = 'GUAGA'
LINKER = 'U'
LOOPS = [
    'UCUU',  # AsCas12a
    'UGUU',  # FnCas12a
    'UAAGU',  # LbCas12a
    'UGUUU'  # MbCas12a
]


# %% [markdown]
# Main functions


def generate_stems(
    max_size: int = 16,  # maximum total number nucleotides in the stem
    min_stem_pairs: int = 2  # min number of paired bases in a row to say that a generated loop forms a self-dimer; such sequences are discarded
) -> set[tuple[str]]:
    """
    Generates all possible stem structures up to max_size
    """
    # Generate stems with up to max_size // 2 pairs (no loose bases yet!)
    stem_fragments = ['P' * i for i in range(min_stem_pairs, 2 * min_stem_pairs)]

    all_pairs = set()
    for n_pairs in range(max_size // 2 // min_stem_pairs + 1):
        base_stem = '|' * n_pairs
        max_unpaired = max_size - len(base_stem) * min_stem_pairs

        out = set([base_stem])
        stems = set([base_stem])

        # Take stem and add one loose base at all possible locations,
        # then add another one to the resulting stems at all possible locations, and so on.
        # If current stem is of size, then we can at most add max_size - size loose pairs
        while len(stems) > 0:
            stem = stems.pop()
            # Insert a single base at all possible locations
            for i in range(1, len(stem)):
                new_stem = stem[:i] + '.' + stem[i:]
                if len(new_stem) - len(base_stem) <= max_unpaired:
                    stems.add(new_stem)
                    out.add(new_stem)  # add to the final list

        # Make all combinations of stem1 and stem2 that are <= max_size
        for o1, o2 in itertools.product(out, repeat=2):
            if len(o1 + o2) - 2 * n_pairs > max_size:
                continue

            for repls in itertools.product(stem_fragments, repeat=n_pairs):
                new_stem1 = o1
                new_stem2 = o2
                for repl in repls:
                    new_stem1 = new_stem1.replace('|', repl, 1).replace('P', '(')
                    new_stem2 = new_stem2.replace('|', repl, 1).replace('P', ')')
                if len(new_stem1 + new_stem2) > max_size:
                    continue
                all_pairs.add((new_stem1, new_stem2[::-1]))

    return all_pairs


def generate_loops(
    size,  # loop length
    min_stem_pairs=2,  # min number of paired bases in a row to say that a generated loop forms a self-dimer; such sequences are discarded
    min_loop_size=3,  # min number of bases between the paired bases required to form a self-dimer; such sequences are discarded
    max_homopolymer_size=5,  # max number of identical bases in a row
    frac_window=10,  # sliding window size to check fractions
    min_frac=.2,  # min fraction of each base required within the window
    max_frac=.3,  # max fraction of each base required within the window
) -> set[str]:
    """
    Generates all possible loops of a given size

    Verifies that:
    - There are no homopolymers longer than the given limit
    - Within each sliding window, each base occurs no less and no more that the given limits
    """
    def get_allowed_bases(
        seq: str,
        allowed='ATGC'
    ) -> str:
        """
        Returns all possible bases that, when added to seq, would not cause self-dimerization
        """
        # We'll go over all k-mers from 5' to 3' and take note which bases would cause self-dimerization
        # Find the last position of k-mers to check
        end = min(
            # Self-dimer can only start forming at
            # min_loop_size + 2 * min_stem_pairs from the end
            # E.g., when filling in the 7th position in this:
            # 01234567
            # GATCCCA_
            # we have to check dimerization with 0 and 1 only
            len(seq) - min_loop_size - 2 * min_stem_pairs + 2,
            # This is the last position from which we still get a full k-mer
            len(seq) - min_stem_pairs + 1
        )
        # go over all k-mer positions
        for i in range(end):
            # for this k-mer, its reverse complement k-mer
            # is not allowed (otherwise a self-dimer would form)...
            disallowed_kmer = rev_compl[seq[i: i + min_stem_pairs]]
            # ...so we cannot complete seq with the last base of this disallowed k-mer
            allowed = allowed.replace(disallowed_kmer[-1], '')
            if len(allowed) == 0:
                break
        return allowed

    def check_frac(
            seq: str,  # sequence to check
            size: int | None = None  # size of the total sequence
    ) -> bool:
        """
        Checks if the fraction of each base is between min_frac and max_frac
        """
        if size is None:
            size = len(seq)

        # min number of bases that must be present
        min_count = int(min_frac * (size - 1)) + 1
        # max number of bases that must be present
        max_count = int(max_frac * size)

        if size >= frac_window:
            for base in 'ATGC':
                # Slide the window
                for i in range(len(seq) - frac_window + 1):
                    counts = collections.Counter(seq[i: i + frac_window])
                    if counts[base] > max_count:
                        return False
                    # If we have a certain count of bases at the moment,
                    # in the best case we can have size - len(seq) more
                    # of them at the end of generation. So that number
                    # must be at least min_count or more.
                    if counts[base] + size - len(seq) < min_count:
                        return False
            # If we're here, it means that the sequence passed all sliding window tests
            return True
        else:  # No check for sequences shorter than the window
            return True

    bases = 'ATGC'
    # Precompute all reverse complements to our k-mer (k=min_stem_pairs)
    compl_table = str.maketrans(bases, 'TACG')
    rev_compl = {''.join(kmer): ''.join(kmer).translate(compl_table)[::-1] for kmer in itertools.product(bases, repeat=min_stem_pairs)}

    # Start from a single base and grow until the required loop length is reached
    seqs = set(list(bases))
    for _ in range(1, size):
        for seq in seqs.copy():
            seqs.remove(seq)
            allowed_bases = get_allowed_bases(seq)
            for base in allowed_bases:
                new_seq = seq + base
                # Any too long homopolymers?
                if new_seq[-max_homopolymer_size:] == base * max_homopolymer_size:
                    continue
                # Check if fractions are alright
                # Helps prune the options and speed up
                if check_frac(new_seq, size=size):
                    seqs.add(new_seq)

    # Now check if base fractions of the final sequences are acceptable using a sliding window
    if size >= frac_window:
        final = set()
        for seq in seqs:
            if check_frac(seq):
                final.add(seq)
    else:
        final = seqs

    return final


def loop_deletions(sequence: str) -> set[str]:
    """
    Generates all possible deletions of a given sequence
    """
    out = set()
    seqs = set([sequence])
    for _ in range(1, len(sequence)):
        for seq in seqs.copy():
            seqs.remove(seq)
            for i in range(len(seq)):
                new_seq = seq[:i] + seq[i + 1:]
                seqs.add(new_seq)
                out.add(new_seq)

    return out


def draw(
    structure: str,
    colors: list | None = None,
    node_size: int = 10,
    ax: plt.Axes | None = None
) -> None:
    """
    Draw a structure
    """
    if colors is None:
        colors = ['skyblue'] * len(structure)
    if ax is None:
        ax = plt.subplot(111)

    # Place nucleotides in a graph with pairings obtained from the structure
    G = nx.Graph()
    stack = []
    for i, (char, color) in enumerate(zip(structure, colors, strict=True)):
        G.add_node(i, color=color)
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

    _ = nx.draw(G, rotated_pos, with_labels=False, node_size=node_size,
                node_color=[G.nodes[n]['color'] for n in G.nodes],
                edge_color='lightgray',
                ax=ax)


# %% [markdown]
# # Generate all possible stem structures

stems = generate_stems(max_size=16)

# Plot
df_stems = pandas.DataFrame(stems, columns=['stem1', 'stem2'])
df_stems['size'] = df_stems.apply(lambda x: len(x['stem1'] + x['stem2']), axis=1)
agg = df_stems.groupby('size').size()
agg.plot.bar()
agg


# %% [markdown]
# Draw stem structural variation samples

uq_sizes = list(set(df_stems['size']))
n_samples = 5
fig, axes = plt.subplots(nrows=n_samples + 1, ncols=len(uq_sizes), figsize=(12, 8))
rng = random.Random(0)
for col, size in enumerate(uq_sizes):
    samples = df_stems[df_stems['size'] == size]
    inds = samples.index.tolist()
    rng.shuffle(inds)

    axes[0][col].text(.5, .5, str(size), ha='center', va='center', fontsize=12)

    i = 0
    for row, idx in enumerate(inds[:n_samples]):
        stem1 = samples.loc[idx, 'stem1']
        stem2 = samples.loc[idx, 'stem2']
        structure = f'....{stem1}...{stem2}.'
        colors = ['skyblue'] * 4 + ['salmon'] * len(stem1) + ['skyblue'] * 3 + ['salmon'] * len(stem2) + ['skyblue']
        draw(structure, colors=colors, node_size=10, ax=axes[row + 1][col])

    for row in range(n_samples + 1):
        axes[row][col].axis('off')

plt.savefig('structures.png', bbox_inches='tight', pad_inches=0, transparent=False, dpi=300)


# %% [markdown]
# # Generate all possible loops

loops = []
for size in tqdm.trange(1, 12):
    loops += generate_loops(size)

df_loops = pandas.DataFrame(loops, columns=['loop'])
df_loops['size'] = df_loops['loop'].apply(lambda x: len(x))
agg = df_loops.groupby('size').size()
agg.plot.bar(y='n_norm')
agg


# %% [markdown]
# # Wild-type deletions
# All possible loop deletions of the wild-type Cas12a

wt_loop_dels = set.union(*[loop_deletions(loop) for loop in LOOPS])
print('Number of WT loop deletions:', len(wt_loop_dels))


# %% [markdown]
# All possible tail deletions of the wild-type Cas12a

wt_tail_dels = loop_deletions('AAUU')
print('Number of WT tail deletions:', len(wt_tail_dels))



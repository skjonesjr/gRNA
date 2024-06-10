import collections, itertools
import tqdm


def stem_structures(
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
        # Given current stem's size, we can have at most add max_size - size loose pairs
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


def stem_sequences(
    structure: str,
    loops_dict: dict[int, list[str]]
) -> set[str]:

    pairs = ['AU', 'CG', 'UA', 'GC']

    parts = []
    indices = []
    stack = []
    dot_stack = []
    for i, char in enumerate(structure):
        if char == '.':
            dot_stack.append(i)
        else:
            if len(dot_stack) > 0:
                indices.extend(dot_stack)
                parts.append(loops_dict[len(dot_stack)])
                dot_stack = []

            if char == '(':
                stack.append(i)
            elif char == ')':
                j = stack.pop()
                indices.extend([j, i])
                parts.append(pairs)

    if len(dot_stack) > 0:
        indices.extend(dot_stack)
        parts.append(loops_dict[len(dot_stack)])

    out = set()
    for combo in tqdm.tqdm(itertools.product(*parts), desc='Stem sequences'):
        flat_combo = ''.join(combo)
        seq = ''.join([nt for _, nt in sorted(zip(indices, flat_combo))])
        out.add(seq)

    return out


def loop_sequences(
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
        allowed='AUGC'
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
            # GAUCCCA_
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
            for base in 'AUGC':
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

    bases = 'AUGC'
    # Precompute all reverse complements to our k-mer (k=min_stem_pairs)
    compl_table = str.maketrans(bases, 'UACG')
    rev_compl = {''.join(kmer): ''.join(kmer).translate(compl_table)[::-1] for kmer in itertools.product(bases, repeat=min_stem_pairs)}

    if size == 0:
        return set('')

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

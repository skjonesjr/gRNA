from ntgen import generation


def loop_to_loop(sequence: str) -> set[str]:
    """
    Generates all possible deletions of a given sequence
    """
    out = set([''])
    seqs = set([sequence])
    for _ in range(1, len(sequence)):
        for seq in seqs.copy():
            seqs.remove(seq)
            for i in range(len(seq)):
                new_seq = seq[:i] + seq[i + 1:]
                seqs.add(new_seq)
                out.add(new_seq)

    return out


def stem_to_loop(
    sequence: tuple[str, str] | list[str, str]
) -> set[str]:
    """
    Deletes stem in all possible ways to result in a loop (i.e., a structure that does not form a self-dime)
    """
    stem = sequence[0] + sequence[1]
    loops = [generation.loop_sequences(s, min_stem_pairs=1, min_frac=0, max_frac=1) for s in range(len(stem) + 1)]
    loops = set.union(*loops)
    out = loop_to_loop(stem)
    for loop in out.copy():
        if loop not in loops:
            out.remove(loop)
    return out


def stem_to_stem(
    sequence: tuple[str, str] | list[str, str],
) -> set[str]:
    """
    Generates all possible sequences of all possible deletions of a given structure
    """
    translate = {n1: n2 for n1, n2 in zip('AUGC', 'UACG')}

    # Generate possible stem structures of this size
    stem_size = len(sequence[0] + sequence[1])
    structures = generation.stem_structures(max_size=stem_size)

    # Generate possible loop sequences that we'll use for checking self-dimers
    loops = [generation.loop_sequences(s, min_stem_pairs=1, min_frac=0, max_frac=1) for s in range(stem_size + 1)]

    # Generate possible deletions for each strand of the stem
    seqs1 = set([sequence[0]]) | loop_to_loop(sequence[0])
    seqs2 = set([sequence[1]]) | loop_to_loop(sequence[1])

    out = set()
    # Go over all structures and check matches to the deletions
    for stem1, stem2 in structures:
        for seq1 in seqs1:
            # Check 1: Do lengths match for stem1?
            if len(seq1) != len(stem1):
                continue

            for seq2 in seqs2:
                # Check 2: Do lengths match for stem2?
                if len(seq2) != len(stem2):
                    continue

                stack = []
                nodes = list(seq1 + seq2)
                # Now get base pairings for this structure
                for i, (char, nt) in enumerate(zip(stem1 + stem2, seq1 + seq2)):
                    if char == '(':
                        stack.append(i)
                    elif char == ')':
                        j = stack.pop()
                        nodes[i] = (nt, nodes[j])
                        nodes[j] = (nodes[j], nt)

                # Check that loops don't form self-dimers
                # First, let's find all groups of unpaired bases, e.g.,
                # in the structure ((12(())345))
                # where unpaired bases were replaced with numbers
                # we want to find the group consisting of 12345
                # and check that the do not form a self-dimer

                # We'll go from each side of the stem
                # looking for these loops.
                # Thus, we add dummy nodes | to make stems equal length
                size = max(len(stem1), len(stem2))
                stem1_ = stem1 + '|' * (size - len(stem1))
                seq1_ = seq1 + '|' * (size - len(seq1))
                stem2_ = '|' * (size - len(stem2)) + stem2
                seq2_ = '|' * (size - len(seq2)) + seq2

                dot_groups = []
                stack1 = []
                stack2 = []
                for s1, n1, s2, n2 in zip(stem1_, seq1_, stem2_[::-1], seq2_[::-1]):
                    if s1 == '.':
                        stack1.append(n1)
                    if s2 == '.':
                        stack2.append(n2)
                    if s1 == '(' and s2 == ')':  # end of dot group
                        if len(stack1) > 0 or len(stack2) > 0:
                            dot_groups.append((stack1, stack2[::-1]))
                            stack1 = []
                            stack2 = []

                # Now go over all groups and check that they are within the possible loop structures
                for stack1, stack2 in dot_groups:
                    group = ''.join(stack1 + stack2)
                    if group not in loops:
                        continue

                # nodes is a list of tuples (base pairs) and single chars (unpaired bases)
                for ns in nodes:
                    if isinstance(ns, tuple):
                        # A pair must be composed of a nucleotide and its complement
                        if ns[0] != translate[ns[1]]:
                            break
                else:  # Everything matches, store this sequence
                    out.add(((seq1, seq2), (stem1, stem2)))

    out.discard((('', ''), ('', '')))
    return out

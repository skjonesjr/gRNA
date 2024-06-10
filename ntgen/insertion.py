import itertools
import tqdm
import Bio.Align

from ntgen import generation


def loop(seq: str, size=1) -> set:
    """
    All possible sequences where `size` of bases are inserted.
    """
    out = []
    inds_combos = itertools.combinations_with_replacement(range(len(seq) + 1), size)
    nts_combos = itertools.combinations_with_replacement('AUGC', size)
    for inds, nts in itertools.product(inds_combos, nts_combos):
        new_seq = seq
        for idx, nt in zip(sorted(inds, reverse=True), nts):
            new_seq = new_seq[:idx] + nt + new_seq[idx:]
        out.append(new_seq)
    return set(out)


def stem(stem1: str,
         stem2: str,
         max_dist: int = 5,
         stem_structures: set | None = None,
         min_stem_pairs: int = 2
         ):
    """
    Given two stem sequences, generate all possible insertions that result in one of the accepted stem structures.
    """

    compl = str.maketrans('ATUCG', 'UAAGC')
    aligner = Bio.Align.PairwiseAligner(gap_score=-.5)
    if stem_structures is None:
        stem_structures = generation.stem_structures(max_size=len(stem1 + stem2) + max_dist)

    insertions1 = [(s, loop(stem1, size=s)) for s in range(max_dist + 1)]
    insertions2 = [(s, loop(stem2, size=s)) for s in range(max_dist + 1)]
    out = set()
    for s1, stems1 in insertions1:
        for s2, stems2 in insertions2:
            if s1 + s2 > max_dist:
                continue

            for stem1, stem2 in tqdm.tqdm(set(itertools.product(stems1, stems2))):
                rev_stem2 = stem2.translate(compl)[::-1]
                for alignment in sorted(aligner.align(stem1, rev_stem2)):
                    struct1 = ''
                    struct2 = ''
                    for i, j in alignment.indices.T:
                        if i == -1:
                            struct2 += '.'
                        elif j == -1:
                            struct1 += '.'
                        elif stem1[i] == rev_stem2[j]:
                            struct1 += '('
                            struct2 += ')'
                        else:
                            struct1 += '.'
                            struct2 += '.'
                        if i == -1 and j == -1:
                            raise ValueError

                    if (struct1, struct2) not in stem_structures:
                        continue

                    out.add((stem1, stem2, struct1, struct2[::-1]))

    return out

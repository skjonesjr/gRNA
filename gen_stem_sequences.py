import multiprocessing, csv
import tqdm
from ntgen import generation, secondary


max_size = 16
min_stem_pairs = 2


def main(i, struct1, struct2):

    loops_dict = {}
    for size in range(1, max_size - 4 * min_stem_pairs + 1):
        loops_dict[size] = generation.loop_sequences(size)

    sequences = generation.stem_sequences(struct1 + struct2, loops_dict)

    f = open(f'sequences_{i:03d}.csv', 'w', newline='')
    writer = csv.writer(f, lineterminator='\n')
    writer.writerow(['structure1', 'structure2', 'sequence1', 'sequence2'])
    # writer.writerow(['structure', 'sequence', 'predicted_structure', 'mfe', 'delta'])

    desc = f'MFE for {i}'
    for seq in tqdm.tqdm(sequences, desc=desc):
        # pred_structure, mfe, delta = secondary.compute_mfe(seq)
        # writer.writerow([stem1, seq, pred_structure, mfe, delta])
        writer.writerow([struct1, struct2,
                         seq[:len(struct1)], seq[len(struct1):]
                         ])

    f.close()


if __name__ == '__main__':
    n_processes = multiprocessing.cpu_count() - 1
    stem_structures = generation.stem_structures(max_size=max_size, min_stem_pairs=min_stem_pairs)
    args = [(i, struct1, struct2) for i, (struct1, struct2) in enumerate(sorted(list(stem_structures)))]

    with multiprocessing.Pool(n_processes) as pool:
        pool.starmap(main, args)


# df = []
# for result in results:
#     df.extend(result)

# df = pandas.DataFrame(df, columns=['structure', 'sequence', 'mfe', 'delta'])
# df.to_csv('stem_sequences.csv', index=False)

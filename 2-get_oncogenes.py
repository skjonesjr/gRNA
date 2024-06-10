# %% [markdown]
# # Get Cas12a targets in oncogenes

from pathlib import Path
import re

import tqdm
import matplotlib.pyplot as plt
import pandas
import mygene
from Bio import Entrez, SeqIO


with open('.env') as f:
    CONFIG = dict([line.strip().split('=', maxsplit=1) for line in f.readlines()])
ROOT = Path(CONFIG['ROOT'])

# Must define email address in case NCBI wants to contact you
Entrez.email = 'joneslabeu@gmail.com'


# %% [markdown]
# # Find oncogenes
#
# Download the list of oncogenes from CancerMine

onco_df = pandas.read_csv('https://zenodo.org/records/7689627/files/cancermine_collated.tsv', sep='\t')
onco_df


# %% [markdown]
# Sort most cited unique genes

top_genes = (onco_df
             .groupby('gene_entrez_id')
             .first()
             .reset_index()
             .sort_values(by='citation_count', ascending=False)
             )
plt.plot(onco_df['citation_count'].iloc[:1000])
top_genes


# %% [markdown]
# # Get sequences of all human chromosomes

name = 'Homo_sapiens.GRCh38.dna.chromosome.{}.fa'

chromosomes = {}
for i in tqdm.tqdm(list(range(1, 23)) + list('XY')):
    chromosome = list(SeqIO.parse(ROOT / 'human_genome' / name.format(i), 'fasta'))
    assert len(chromosome) == 1
    chromosomes[str(i)] = chromosome[0]


# %% [markdown]
# # Get a list of genes with exon start and stop positions for all genes

mg = mygene.MyGeneInfo()
top_n_genes = top_genes.iloc[:200]  # go over that many top oncogenes only
gene_list = mg.querymany(
    top_n_genes['gene_entrez_id'],
    scopes='entrezgene',
    fields='exons',
    species='human'
)

assert len(gene_list) == len(top_n_genes)


# %% [markdown]
# # Extract target sequences

PAMs = ['TTTA', 'TTTG', 'TTTC']
target_size = 21  # target size in our construct

targets = []
missing_chromosomes = set()
for gene, name in zip(gene_list, top_n_genes['gene_normalized']):
    gene_id = gene['query']

    if 'exons' not in gene:
        raise ValueError(f'Gene {gene_id} has no annotated exons')

    # There are many alternative splicings – go over them all
    for exon_no, transcript in enumerate(gene['exons']):
        chrom = transcript['chr']
        # Find exons that are within CDS (cause not all exons are translated),
        # match PAM, and have the correct length
        for start, end in transcript['position']:
            if start >= transcript['cdsstart'] and end <= transcript['cdsend']:
                # Get exon's sequence
                if chrom not in chromosomes:
                    missing_chromosomes.add(chrom)
                    continue
                exon_seq = str(chromosomes[chrom][start: end + 1].seq)
                for pam in PAMs:
                    # Find the PAM in the sequence
                    for match in re.finditer(pam, exon_seq):
                        # Target starts right after PAM
                        target_start = match.end()
                        target = exon_seq[target_start: target_start + target_size]
                        # Check that the target is not outside the CDS
                        if len(target) != target_size:
                            continue
                        targets.append({
                            'chromosome': chrom,
                            'gene_id': gene_id,
                            'name': name,
                            'transcript_id': transcript['transcript'],
                            'exon_no': exon_no,
                            'PAM': pam,
                            'start': start + target_start,
                            'end': start + target_start + target_size,
                            'target': target,
                        })
                        exon_no += 1

targets = pandas.DataFrame(targets)
targets.to_csv(ROOT / f'cas12a_targets_in_top{len(top_n_genes)}_oncogenes.csv', index=False)
print('Missing chromosomes:', missing_chromosomes)
targets


# %% [markdown]
# Alternative method to obtain target sequences that does not require to have the human genome downloaded but is MUCH slower

PAMs = ['TTTA', 'TTTG', 'TTTC']
target_size = 21
mg = mygene.MyGeneInfo()

targets = []
for gene, name in zip(gene_list, top_n_genes['gene_normalized']):
    gene_id = gene['query']

    if 'exons' not in gene:
        raise ValueError(f'Gene {gene_id} has no annotated exons')

    # There are many alternative splicings – go over them all
    for transcript in gene['exons']:
        transcript_id = transcript['transcript']

        # Get the transcript from NCBI
        # Automatically limited to the required 3 requests/sec
        handle = Entrez.efetch(
            db='nucleotide',
            id=transcript_id,
            rettype='gb',
            retmode='text'
        )
        record = SeqIO.read(handle, 'genbank')
        handle.close()

        # Find where the coding sequence starts
        # Not all variants will have a coding sequence!
        for feat in record.features:
            if feat.type == 'CDS':
                cds_start = feat.location.start
                cds_end = feat.location.end
                break

        # Find exons that are within CDS (cause not all exons are translated),
        # match PAM, and have the correct length
        exon_no = 0
        for feat in record.features:
            if feat.type == 'exon' and feat.location.start >= cds_start and feat.location.end <= cds_end:
                # Get exon's sequence
                seq = feat.extract(record).seq
                for pam in PAMs:
                    # Find the PAM in the sequence
                    for match in re.finditer(pam, seq):
                        start = match.end()  # start of the target sequence
                        spacer = seq[start: start + target_size]
                        # Check that the spacer is not outside the CDS
                        if len(spacer) != target_size:
                            continue
                        targets.append({
                            'gene_id': gene_id,
                            'name': name,
                            'transcript_id': transcript_id,
                            'exon_no': exon_no,
                            'PAM': pam,
                            'start': start,
                            'end': start + target_size,
                            'spacer': str(spacer),
                        })
                        exon_no += 1

targets = pandas.DataFrame(targets)
targets.to_csv(ROOT / f'cas12a_targets_in_top{len(gene_list)}_oncogenes.csv', index=False)
targets

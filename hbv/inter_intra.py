from .preset import DB
from preset.fasta import load_fasta
from .genotype_distance import calc_intra_distance
from .genotype_distance import calc_inter_distance


def get_inter_intra(
        folder=DB / 'genotype', exclude_genotype=['RF', 'overall']):

    seq_files = []
    for i in folder.iterdir():
        if i.name.endswith('_seq.fasta'):
            genotype = i.stem.replace('_seq', '')
            if genotype in exclude_genotype:
                continue

            seq_files.append((genotype, i))

    intra_dist = {}

    all_aligned_RT = []
    for genotype, i in seq_files:
        aligned_RT = load_fasta(i, with_name=False)
        all_aligned_RT.extend(aligned_RT)

        print(genotype, len(aligned_RT))

        dist = calc_intra_distance(aligned_RT)
        intra_dist[genotype] = dist

    # calc_covariation(all_aligned_RT)
    calc_inter_distance(DB / 'aligned_consensus.fasta', intra_dist)

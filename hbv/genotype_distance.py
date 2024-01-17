from .preset import DB
from collections import defaultdict
from preset.file_format import dump_csv
from preset.fasta import load_fasta
from itertools import combinations
from statistics import mean


def calc_inter_distance(consensus_file, intra_dist):

    consensus_list = load_fasta(consensus_file)
    consensus_list.pop('overall')
    consensus_list = [
        (g, s)
        for g, s in consensus_list.items()
    ]

    inter_distance = defaultdict(dict)
    for (g1, s1), (g2, s2) in combinations(consensus_list, 2):
        # print(f'inter distance {g1}, {g2},', calc_seq_distance(s1, s2))
        inter_distance[g1][g2] = calc_seq_distance(s1, s2)
        inter_distance[g2][g1] = calc_seq_distance(s1, s2)

    for g, dist in intra_dist.items():
        if g not in inter_distance:
            continue
        inter_distance[g][g] = dist

    report = []
    for g1, value in inter_distance.items():
        rec = {
            'genotype': g1
        }
        for g2 in sorted(value.keys()):
            rec[g2] = f"{round(value[g2] * 100, 2)}%"

        report.append(rec)

    dump_csv(DB / 'inter_distance.csv', report)


def calc_intra_distance(aligned_seq):
    distance_list = []
    for i, j in combinations(aligned_seq, 2):
        distance_list.append(calc_seq_distance(i, j))

    return mean(distance_list)


def calc_seq_distance(seq1, seq2):
    assert len(seq1) == len(seq2)
    length = len(seq1)
    diff = 0
    for i, j in zip(seq1, seq2):
        if i != j:
            diff += 1

    return diff / length

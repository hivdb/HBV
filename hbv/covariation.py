from .preset import DB
from collections import defaultdict
from preset.file_format import dump_csv
from itertools import combinations
from scipy.stats import chi2_contingency


def calc_covariation(aligned_RT, round_number=0):
    conserve_mut = get_conserve_mut(aligned_RT)
    pos_mut = calc_mut_count(aligned_RT)
    paired_pos_mut = calc_paired_mut_count(aligned_RT)

    pcnt_cut = 1 / (10 ** (round_number + 2))

    report = []
    total = len(aligned_RT)
    for (p1, p2), cnt in paired_pos_mut.items():
        if p1 in conserve_mut:
            continue
        if p2 in conserve_mut:
            continue
        if (pos_mut[p1] / total) < pcnt_cut:
            continue
        if (pos_mut[p2] / total) < pcnt_cut:
            continue

        if (pos_mut[p1] / total) > 0.5:
            continue
        if (pos_mut[p2] / total) > 0.5:
            continue

        chisq = chi2_contingency([
                [cnt, pos_mut[p1] - cnt],
                [
                    pos_mut[p2] - cnt,
                    total - cnt - (pos_mut[p1] - cnt) - (pos_mut[p2] - cnt)]
            ])

        report.append({
            'mut1': p1,
            'mut2': p2,
            'TT': cnt,
            'TF': pos_mut[p1] - cnt,
            'FT': pos_mut[p2] - cnt,
            'FF': total - cnt - (pos_mut[p1] - cnt) - (pos_mut[p2] - cnt),
            'chisq': chisq.statistic,
            'p': chisq.pvalue,
        })

    dump_csv(DB / 'covariation.csv', report)


def get_conserve_mut(aligned_RT):
    pos_mut = defaultdict(set)
    for seq in aligned_RT:
        for ofst, mut in enumerate(list(seq)):
            pos = ofst + 1
            pos_mut[pos].add(mut)

    conserve_mut = [
        f"{pos}{list(mut_list)[0]}"
        for pos, mut_list in pos_mut.items()
        if len(mut_list) == 1
    ]

    return conserve_mut


def calc_mut_count(aligned_RT):
    pos_mut = defaultdict(int)
    for seq in aligned_RT:
        for ofst, mut in enumerate(list(seq)):
            if len(mut) > 1:
                continue
            if mut == 'X':
                continue

            pos = ofst + 1
            pos_mut[f"{pos}{mut}"] += 1

    return pos_mut


def calc_paired_mut_count(aligned_RT):

    paired_mut = defaultdict(int)
    for seq in aligned_RT:
        pos_mut = []
        for ofst, mut in enumerate(list(seq)):
            if len(mut) > 1:
                continue
            if mut == 'X':
                continue

            pos = ofst + 1
            pos_mut.append(f"{pos}{mut}")

        for p1, p2 in combinations(pos_mut, 2):
            paired_mut[(p1, p2)] += 1

    return paired_mut

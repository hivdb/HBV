from Bio import AlignIO
from .preset import DB
from collections import defaultdict
from preset.file_format import dump_csv
from preset.fasta import dump_fasta
from preset.fasta import load_fasta
from operator import itemgetter
from scipy.stats import entropy
from itertools import combinations
from statistics import mean


HBVDB_GENOTEYPE_START_STOP = {
    'A': (351, 696),
    'B': (352, 696),
    'C': (363, 714),
    'D': (351, 700),
    'E': (348, 692),
    'F': (349, 697),
    'G': (346, 690),
    'H': (346, 690),
    'RF': (359, 703),
}


def get_prevalence():

    folder = DB / 'hbvdb'
    genotype_files = []
    for i in folder.iterdir():
        if i.suffix != '.clu':
            continue

        genotype = i.name.split('_', 1)[0]
        genotype_files.append((genotype, i))

    genotype_files.sort(key=lambda x: x[0])

    # genotype_detect_RT(genotype_files)

    prevalence = []
    intra_dist = {}
    for genotype, i in genotype_files:
        prev, dist = get_genotype_prevalence(genotype, i)
        prevalence.extend(prev)
        intra_dist[genotype] = dist

    dump_csv(DB / 'hbvdb_prevalence.csv', prevalence)

    dump_pos_mut_by_genotype(prevalence)
    dump_pos_mut_by_mutation(prevalence)

    prevalence = get_overall_prevalance(prevalence)

    dump_csv(DB / 'hbvdb' / 'overall_prev.csv', prevalence)

    cons_seq = get_cons_seq(prevalence)

    dump_fasta(DB / 'hbvdb' / 'overall_cons.fasta', {'overall': cons_seq})

    align_consensus(DB / 'hbvdb')

    calc_inter_distance(DB / 'aligned_consensus.fasta', intra_dist)


def dump_pos_mut_by_genotype(prevalence, exclude_genotype=['RF']):

    prevalence = [
        i
        for i in prevalence
        if i['genotype'] not in exclude_genotype
    ]

    genotypes = sorted(list(set(
        i['genotype']
        for i in prevalence
    )))

    pos_mut_list = set([
        (i['pos'], i['mut'])
        for i in prevalence
    ])

    pos_mut_list = {
        (pos, mut): {
            g: 0
            for g in genotypes
        }
        for pos, mut in pos_mut_list
    }

    for i in prevalence:
        pos_mut_list[(i['pos'], i['mut'])][i['genotype']] = i['pcnt']

    report = []

    for (pos, mut), value in pos_mut_list.items():
        record = {
            'pos': pos,
            'mut': mut
        }
        for g, pcnt in value.items():
            record[g] = pcnt

        report.append(record)

    report.sort(key=itemgetter('pos', 'mut'))

    dump_csv(DB / 'genotype_compare.csv', report)


def dump_pos_mut_by_mutation(prevalence, exclude_genotype=['RF']):

    prevalence = [
        i
        for i in prevalence
        if i['genotype'] not in exclude_genotype
    ]

    muts = sorted(list(set(
        i['mut']
        for i in prevalence
        if i['mut'] not in ['X', 'del', 'ins']
    )))

    pos_genotype = set([
        (i['pos'], i['genotype'])
        for i in prevalence
    ])

    pos_genotype = {
        (pos, genotype): {
            m: 0
            for m in muts
        }
        for pos, genotype in pos_genotype
    }

    for i in prevalence:
        if i['mut'] in ['X', 'del', 'ins']:
            continue
        pos_genotype[(i['pos'], i['genotype'])][i['mut']] = i['pcnt']

    report = []

    for (pos, genotype), value in pos_genotype.items():
        record = {
            'pos': pos,
            'genotype': genotype
        }
        for m, pcnt in value.items():
            record[m] = pcnt

        e = entropy([
            pcnt
            for _, pcnt in value.items()
        ], base=2)

        record['entropy'] = e

        report.append(record)

    report.sort(key=itemgetter('pos', 'genotype'))

    dump_csv(DB / 'genotype_compare_by_mut.csv', report)


def align_consensus(folder, exclude_genotype=['RF']):

    fasta_files = []

    for i in folder.iterdir():
        if i.suffix != '.fasta':
            continue
        fasta_files.append(i)

    fasta_files.sort()

    consensus_list = {}
    for i in fasta_files:
        consensus_list.update(load_fasta(i))

    consensus_list = {
        k: v
        for k, v in consensus_list.items()
        if k not in exclude_genotype
    }

    dump_fasta(DB / 'aligned_consensus.fasta', consensus_list, 'overall')


def get_overall_prevalance(prevalence, exclude_genotype=['RF']):

    prevalence = [
        i
        for i in prevalence
        if i['genotype'] not in exclude_genotype
    ]

    genotype_num_seq = {}
    for i in prevalence:
        genotype = i['genotype']
        if genotype in exclude_genotype:
            continue

        total = i['total']
        genotype_num_seq[genotype] = total

    total = sum([
        t
        for g, t in genotype_num_seq.items()
    ])

    profile = defaultdict(dict)
    for i in prevalence:
        pos = i['pos']
        mut = i['mut']
        num = i['num']

        if pos not in profile:
            profile[pos] = defaultdict(int)

        profile[pos][mut] += num

    pos_cons = collect_consensus(profile)

    return get_prevalence_by_profile('all', profile, pos_cons, total)


def genotype_detect_RT(genotype_files):

    for genotype, i in genotype_files:

        aligned_seqs = get_aligned_sequence(i)

        print(genotype)

        if genotype in HBVDB_GENOTEYPE_START_STOP:
            start, stop = HBVDB_GENOTEYPE_START_STOP[genotype]
            print(genotype, 'length', stop - start)
            detect_RT(
                aligned_seqs,
                (start, 0, 10),
                (stop, -10, 0)
            )
        else:
            detect_RT(aligned_seqs)


def get_aligned_sequence(file_path):
    alignment = AlignIO.read(file_path, "clustal")

    seqs = [
        str(i.seq)
        for i in alignment
    ]

    return seqs


def detect_RT(
        aligned_seqs,
        start_window=(346, -5, 5),
        stop_window=(346+344, -5, 5)):

    start_probe = set()
    stop_probe = set()

    for i in aligned_seqs:
        pos, left, right = start_window
        start_probe.add(i[pos + left: pos + right])

        pos, left, right = stop_window
        stop_probe.add(i[pos + left: pos + right])

    pos, left, right = start_window
    print('start', pos, pos + left, pos + right)
    print(start_probe)

    pos, left, right = stop_window
    print('stop', pos, pos + left, pos + right)
    print(stop_probe)


def get_RT(aligned_seqs, start, stop):

    return [
        i[start: stop]
        for i in aligned_seqs
    ]


def get_genotype_prevalence(genotype, genotype_file):

    aligned_seqs = get_aligned_sequence(genotype_file)

    aligned_RT = get_RT(aligned_seqs, *HBVDB_GENOTEYPE_START_STOP[genotype])

    prevalence, aligned_RT = collect_prevalence(genotype, aligned_RT)

    # print(genotype)
    intra_distance = calc_intra_distance(aligned_RT)

    dump_csv(DB / 'hbvdb' / f'{genotype}_prev.csv', prevalence)

    cons_seq = get_cons_seq(prevalence)

    dump_fasta(DB / 'hbvdb' / f'{genotype}_cons.fasta', {genotype: cons_seq})
    # print(genotype, cons_seq, len(cons_seq), cons_seq.count('-'))

    return prevalence, intra_distance


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


def collect_prevalence(genotype, aligned_RT):

    num_total = len(aligned_RT)

    pos_mut = build_pos_mut(aligned_RT)

    pos_cons = collect_consensus(pos_mut)

    merge_pos_map = get_merge_pos_map(pos_cons)
    # print(merge_pos_map)

    pos_mut = build_pos_mut(aligned_RT, merge_pos_map)
    pos_cons = collect_consensus(pos_mut)

    return (
        get_prevalence_by_profile(genotype, pos_mut, pos_cons, num_total),
        get_aligned_seq_after_merge(aligned_RT, merge_pos_map)
    )


def get_aligned_seq_after_merge(aligned_seq, merge_pos_map):
    return merge_pos(aligned_seq, merge_pos_map)


def get_prevalence_by_profile(
        genotype, pos_mut, pos_cons, num_total, round_number=0):

    prevalence = []
    for pos, mut_list in pos_mut.items():
        for mut, num in mut_list.items():

            if mut.replace('-', '') != '':
                mut = mut.replace('-', '')
            else:
                mut = 'del'

            if len(mut) > 1:
                mut = 'ins'

            cons = pos_cons[pos].replace('-', '')

            pcnt = num / num_total

            if pcnt < (1 / (10 ** (round_number + 2))):
                continue

            # if mut == 'X':
            #     continue

            prevalence.append({
                'genotype': genotype,
                'cons': cons,
                'pos': pos,
                'mut': mut,
                'total': num_total,
                'num': num,
                'pcnt': round(pcnt * 100, round_number),
                'is_cons': 'Y' if mut == cons else ''
            })

    return prevalence


def build_pos_mut(aligned_seq, merge_pos_map={}):

    pos_mut = defaultdict(dict)

    if merge_pos_map:
        aligned_seq = merge_pos(aligned_seq, merge_pos_map)

    seq_length = len(aligned_seq[0])

    for seq in aligned_seq:

        for ofst in range(seq_length):
            mut = seq[ofst]
            pos = ofst + 1

            if pos not in pos_mut:
                pos_mut[pos] = defaultdict(int)

            pos_mut[pos][mut] += 1

    return pos_mut


def merge_pos(aligned_seq, merge_pos_map):

    merged_seq = []

    for seq in aligned_seq:
        seq = list(seq)

        for cons_pos, pos_list in merge_pos_map.items():
            cons_ofst = cons_pos - 1
            pos_list = sorted(pos_list)

            for pos in pos_list:
                ofst = pos - 1
                seq[cons_ofst] += seq[ofst]

        merge_pos_list = [
            p
            for _, pos_list in merge_pos_map.items()
            for p in pos_list
        ]

        seq = [
            mut
            for ofst, mut in enumerate(seq)
            if (ofst + 1) not in merge_pos_list
        ]

        merged_seq.append(seq)

    return merged_seq


def collect_consensus(pos_mut):
    pos_cons = {}
    for pos, mut_list in pos_mut.items():
        mut_list = [
            (mut, num)
            for mut, num in mut_list.items()
        ]
        mut_list.sort(key=lambda x: x[-1])
        cons_mut = mut_list[-1][0]
        pos_cons[pos] = cons_mut

    return pos_cons


def get_cons_seq(prevalence):

    cons_seq = {}
    for i in prevalence:
        cons_seq[i['pos']] = i['cons']

    cons_seq = [
        (pos, cons)
        for pos, cons in cons_seq.items()
    ]

    cons_seq.sort(key=lambda x: x[0])

    cons_seq = [
        i[-1]
        for i in cons_seq
    ]

    return ''.join(cons_seq)


def get_merge_pos_map(pos_cons):

    pos_that_cons_is_del = [
        pos
        for pos, cons in pos_cons.items()
        if cons == '-'
    ]

    if not pos_that_cons_is_del:
        return {}

    merge_pos_map = defaultdict(list)
    for pos in pos_that_cons_is_del:
        prev_pos = pos - 1

        while prev_pos > 1:
            mut = pos_cons[prev_pos]
            if mut != '-':
                break
            prev_pos -= 1

        merge_pos_map[prev_pos].append(pos)

    return dict(merge_pos_map)

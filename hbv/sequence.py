from Bio import AlignIO
from .preset import DB
from collections import defaultdict
from preset.fasta import dump_fasta


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


def get_sequences():

    genotype_files = load_alignment_file()

    # genotype_detect_RT(genotype_files)

    all_aligned_RT = []
    for genotype, i in genotype_files:
        aligned_RT = get_aligned_RT(genotype, i)

        all_aligned_RT.extend(aligned_RT)

    dump_fasta(DB / 'genotype' / 'overall_seq.fasta', {
        (idx + 1): ''.join(i)
        for idx, i in enumerate(all_aligned_RT)
    })


def load_alignment_file(folder=DB / 'hbvdb', exclude_genotype=['RF']):

    genotype_files = []
    for i in folder.iterdir():
        if i.suffix != '.clu':
            continue

        genotype = i.name.split('_', 1)[0]
        if genotype in exclude_genotype:
            continue
        genotype_files.append((genotype, i))

    genotype_files.sort(key=lambda x: x[0])

    return genotype_files


def get_aligned_RT(genotype, genotype_file):

    aligned_seqs = get_aligned_sequence(genotype_file)

    aligned_RT = get_RT(aligned_seqs, *HBVDB_GENOTEYPE_START_STOP[genotype])

    aligned_RT = fix_aligned_RT(genotype, aligned_RT)

    dump_fasta(DB / 'genotype' / f'{genotype}_seq.fasta', {
        (idx + 1): ''.join(i)
        for idx, i in enumerate(aligned_RT)
    })

    return aligned_RT


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


def get_aligned_sequence(file_path):
    alignment = AlignIO.read(file_path, "clustal")

    seqs = [
        str(i.seq)
        for i in alignment
    ]

    return seqs


def get_RT(aligned_seqs, start, stop):

    return [
        i[start: stop]
        for i in aligned_seqs
    ]


def fix_aligned_RT(genotype, aligned_RT):

    pos_mut = build_pos_mut(aligned_RT)

    pos_cons = collect_consensus(pos_mut)

    merge_pos_map = get_merge_pos_map(pos_cons)
    # print(genotype, merge_pos_map)

    fixed_aligned_RT = get_aligned_seq_after_merge(aligned_RT, merge_pos_map)

    return fixed_aligned_RT


def build_pos_mut(aligned_seq, merge_pos_map={}, exclude_mut=[]):

    pos_mut = defaultdict(dict)

    if merge_pos_map:
        aligned_seq = merge_pos(aligned_seq, merge_pos_map)

    seq_length = len(aligned_seq[0])

    for seq in aligned_seq:

        for ofst in range(seq_length):
            mut = seq[ofst]
            if mut in exclude_mut:
                continue

            pos = ofst + 1

            if pos not in pos_mut:
                pos_mut[pos] = defaultdict(int)

            pos_mut[pos][mut] += 1

    return pos_mut


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


def get_aligned_seq_after_merge(aligned_seq, merge_pos_map):
    return merge_pos(aligned_seq, merge_pos_map)


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
            process_indel(mut)
            for ofst, mut in enumerate(seq)
            if (ofst + 1) not in merge_pos_list
        ]

        merged_seq.append(seq)

    return merged_seq


def process_indel(mut):

    if mut.replace('-', '') != '':
        mut = mut.replace('-', '')
    else:
        mut = '-'

    if len(mut) > 1:
        mut = '_'

    return mut

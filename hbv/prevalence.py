from Bio import AlignIO
from .preset import DB
from preset.file_format import dump_csv
from preset.fasta import dump_fasta
from preset.fasta import load_fasta
from .genotype import dump_pos_mut_by_genotype
# from .genotype import dump_pos_mut_by_mutation
from preset.table import group_records_by
from collections import Counter
from .preset import AA
import re
from .preset import overall_pcnt_round
from .preset import genotype_pcnt_round
from .sequence import build_pos_mut
from .sequence import collect_consensus
from .usual_muts import calc_num_rear
from preset.file_format import load_csv


def get_prevalence(
        folder=DB / 'genotype',
        exclude_genotype=['RF', 'overall'],
        save_folder=DB / 'prevalence_1',
        filter=None):

    seq_files = []
    for i in folder.iterdir():
        if i.name.endswith('_seq.fasta'):
            genotype = i.stem.replace('_seq', '')
            if genotype in exclude_genotype:
                continue

            seq_files.append((genotype, i))

    prevalence = []

    all_aligned_RT = []
    for genotype, i in seq_files:
        aligned_RT = load_fasta(i, with_name=False)

        if filter:
            aligned_RT = filter(aligned_RT)

        all_aligned_RT.extend(aligned_RT)

        prev = get_prevalence_by_aligned_seq(genotype, aligned_RT)
        prevalence.extend(prev)

    overall_prev = get_overall_prevalance(all_aligned_RT)

    dump_csv(DB / 'genotype' / 'overall_prev.csv', overall_prev)

    # TODO a switch between two different ways of calc over all cons
    # save two version of overall cons to a single folder
    overall_cons_seq = get_cons_seq(overall_prev)
    # overall_cons_seq = get_overall_cons_by_genotype(DB / 'genotype')
    dump_fasta(
        DB / 'genotype' / 'overall_cons.fasta', {'overall': overall_cons_seq})

    prevalence += overall_prev

    for i in prevalence:
        i['overall_cons'] = overall_cons_seq[i['pos'] - 1]

    dump_csv(DB / save_folder / 'hbvdb_prevalence.csv', prevalence)

    dump_csv(DB / save_folder / 'hbvdb_prevalence_view.csv',
             get_prevalence_view(prevalence))

    dump_pos_mut_by_genotype(prevalence, save_folder)

    # dump_pos_mut_by_mutation(prevalence, save_folder)

    align_consensus(DB / 'genotype')


def get_prevalence_view(prevalence):

    non_normal = [
        i
        for i in prevalence
        if not re.match(AA, i['mut'])
    ]

    print('non normal muts', set([
        i['mut']
        for i in non_normal
        if i['mut'] != 'X'
    ]))

    prevalence = [
        i
        for i in prevalence
        if re.match(AA, i['mut'])
    ]

    return prevalence


def get_overall_cons_by_genotype(folder, exclude_genotype=['overall', 'RF']):
    fasta_files = []

    for i in folder.iterdir():
        if i.suffix != '.fasta':
            continue
        fasta_files.append(i)

    aligned_seq = []
    for i in fasta_files:
        for i, j in load_fasta(i).items():
            aligned_seq.append(j)

    num_total = len(aligned_seq)

    pos_mut = build_pos_mut(aligned_seq)

    pos_cons = collect_consensus(pos_mut)

    prevalence = get_prevalence_by_profile(
        'overall', pos_mut, pos_cons, num_total)

    cons_seq = get_cons_seq(prevalence)

    return cons_seq


def align_consensus(folder, exclude_genotype=['RF']):

    fasta_files = []

    for i in folder.iterdir():
        if i.suffix != '.fasta':
            continue
        if 'cons' not in i.name:
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


def get_overall_prevalance(aligned_RT):

    pos_mut = build_pos_mut(aligned_RT)

    pos_cons = collect_consensus(pos_mut)

    num_total = len(aligned_RT)

    return get_prevalence_by_profile(
        'overall', pos_mut, pos_cons, num_total,
        round_number=overall_pcnt_round)


def get_aligned_sequence(file_path):
    alignment = AlignIO.read(file_path, "clustal")

    seqs = [
        str(i.seq)
        for i in alignment
    ]

    return seqs


def get_prevalence_by_profile(
        genotype, pos_mut, pos_cons, num_total,
        round_number=genotype_pcnt_round):

    prevalence = []
    for pos, mut_list in pos_mut.items():
        # num_total = sum([
        #     num
        #     for _, num in mut_list.items()
        # ])
        for mut, num in mut_list.items():

            cons = pos_cons[pos]

            pcnt = num / num_total

            if pcnt < (1 / (10 ** (round_number + 2))):
                continue

            prevalence.append({
                'genotype': genotype,
                'cons': cons,
                'pos': pos,
                'mut': mut,
                'total': num_total,
                'num': num,
                'pcnt': (
                    round(pcnt * 100, round_number)
                    if round_number else round(pcnt * 100)
                ),
                # 'is_cons': 'Y' if mut == cons else ''
            })

    assert_prevalance_report(prevalence)

    return prevalence


def assert_prevalance_report(prevalence):

    for key, rows in group_records_by(prevalence, ['genotype', 'pos']).items():
        mutations = [
            i['mut']
            for i in rows
        ]

        key = dict(key)
        genotype = key['genotype']
        pos = key['pos']
        for mut, cnt in Counter(mutations).items():
            assert cnt == 1, f"{genotype}, {pos}, {mut}"


def get_prevalence_by_aligned_seq(genotype, aligned_RT):

    pos_mut = build_pos_mut(aligned_RT)
    pos_cons = collect_consensus(pos_mut)

    num_total = len(aligned_RT)

    prevalence = get_prevalence_by_profile(
        genotype, pos_mut, pos_cons, num_total)

    dump_csv(DB / 'genotype' / f'{genotype}_prev.csv', prevalence)

    cons_seq = get_cons_seq(prevalence)

    dump_fasta(DB / 'genotype' / f'{genotype}_cons.fasta', {genotype: cons_seq})

    # print(genotype, cons_seq, len(cons_seq), cons_seq.count('-'))

    return prevalence


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


# TODO the num rear highest for each genotype
def seq_filter(seqs, num_rear_highest=6):
    usual_muts = DB / 'prevalence_1' / 'genotype_compare.csv'
    usual_muts = load_csv(usual_muts)
    usual_muts = [
        (int(i['pos']), i['mut'])
        for i in usual_muts
        if i['is_usual'].lower() == 'yes'
    ]

    filtered = []
    for seq in seqs:
        rear_num = len(calc_num_rear(seq, usual_muts))
        if rear_num > num_rear_highest:
            continue
        else:
            filtered.append(seq)

    return filtered

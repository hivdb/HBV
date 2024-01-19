from preset.fasta import load_fasta
from .preset import DB
from preset.file_format import load_csv
from preset.file_format import dump_csv
from preset.plot import plot_hist
from collections import Counter
from operator import itemgetter


def view_usual_muts_by_genotype(
        folder=DB / 'genotype',
        usual_muts=DB / 'prevalence_1' / 'genotype_compare.csv',
        save_folder='prevalence_1',
        exclude_genotype=['RF']):

    fasta_files = {}

    for i in folder.iterdir():
        if i.suffix != '.fasta':
            continue
        if '_seq' not in i.name:
            continue

        genotype = i.stem.replace('_seq', '')
        if i.stem.replace('_seq', '') in exclude_genotype:
            continue

        fasta_files[genotype] = i

    usual_muts = load_csv(usual_muts)
    usual_muts = [
        (int(i['pos']), i['mut'])
        for i in usual_muts
        if i['is_usual'].lower() == 'yes'
    ]
    # print(len(usual_muts))

    for geno, f in fasta_files.items():
        seqs = load_fasta(f, with_name=False)
        rear_num = calc_num_rear_per_seq(seqs, usual_muts)

        hist = get_hist_summary(rear_num)

        x = [
            i[0]
            for i in hist
        ]

        y = [
            i[1]
            for i in hist
        ]

        dump_csv(DB / save_folder / 'rear_muts' / f'{geno}.csv', sorted([
            {
                'num_unusual_mut': i[0],
                'num_seq': i[1]
            }
            for i in hist
        ], key=itemgetter('num_unusual_mut')))

        plot_hist(
            DB / save_folder / 'rear_muts' / f'{geno}.png',
            x, y, geno,
            '# Unusual muts', '# Sequences')


def calc_num_rear_per_seq(seqs, usual_muts):

    rear_num = [
        len(calc_num_rear(seq, usual_muts))
        for seq in seqs
    ]

    return rear_num


def calc_num_rear(seq, usual_muts):

    rear_mut = []

    for ofst, mut in enumerate(seq):
        pos = ofst + 1

        if mut == 'X':
            continue

        if (pos, mut) not in usual_muts:
            rear_mut.append((pos, mut))

    return rear_mut


def get_hist_summary(numbers):

    return [
        (num_mut, num_seq)
        for num_mut, num_seq in Counter(numbers).items()
    ]

from .preset import DB
from preset.file_format import load_csv
from preset.file_format import dump_csv


def prepare_sierra_files(
        folder=DB / 'genotype',
        exclude_genotype=['RF'],
        save_folder=DB / 'Sierra program'):

    prev_files = []
    for i in folder.iterdir():
        if i.name.endswith('_prev.csv'):
            genotype = i.stem.replace('_prev', '')
            if genotype in exclude_genotype:
                continue

            prev_files.append((genotype, i))

    usual_muts = load_usual_muts()

    for g, i in prev_files:
        prevs = [
            {
                'gene': 'RT',
                'position': i['pos'],
                'aa': i['mut'],
                'percent': int(i['num']) / int(i['total']),
                'count': i['num'],
                'total': i['total'],
                'reason': 'PCNT',
                'isUnusual': not ((int(i['pos']), i['mut']) in usual_muts)
            }
            for i in load_csv(i)
        ]
        g = g.replace('overall', 'all')
        dump_csv(save_folder / f'rx-all_genotype-{g}.csv', prevs)


def load_usual_muts():
    usual_muts = DB / 'prevalence_1' / 'genotype_compare.csv'
    usual_muts = load_csv(usual_muts)
    usual_muts = [
        (int(i['pos']), i['mut'])
        for i in usual_muts
        if i['is_usual'].lower() == 'yes'
    ]

    return usual_muts

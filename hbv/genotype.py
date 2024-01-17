from .preset import DB
from preset.file_format import dump_csv
from operator import itemgetter
from scipy.stats import entropy


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

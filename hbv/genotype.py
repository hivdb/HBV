from .preset import DB
from preset.file_format import dump_csv
from preset.table import group_records_by
from operator import itemgetter
from scipy.stats import entropy
from math import ceil


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

    genotype_total = {
        i['genotype']: i['total']
        for i in prevalence
    }

    pos_mut_list = set([
        (i['pos'], i['mut'])
        for i in prevalence
    ])

    pos_mut_list = {
        (pos, mut): {
            g: (0, 0)
            for g in genotypes
        }
        for pos, mut in pos_mut_list
    }

    for i in prevalence:
        pos_mut_list[
            (i['pos'], i['mut'])][i['genotype']] = (
                i['pcnt'], i['num'])

    usual_mut = get_usual_mutation(prevalence)

    report = []

    for (pos, mut), value in pos_mut_list.items():
        record = {
            'pos': pos,
            'mut': mut
        }
        for g, (pcnt, num) in value.items():
            record[f"#{g} ({genotype_total[g]})"] = num
            record[f"%{g}"] = f'{pcnt}%'

        if (pos, mut) in usual_mut:
            record['is_usual'] = 'yes'

        report.append(record)

    report.sort(key=itemgetter('pos', 'mut'))

    dump_csv(DB / 'genotype_compare.csv', report)


def get_usual_mutation(prevalence):

    prevalence = [
        i
        for i in prevalence
        if i['genotype'] != 'overall'
    ]

    cutoff_config = prepare_cutoff_config(prevalence)

    usual_mut = []

    for key, value in group_records_by(prevalence, ['pos', 'mut']).items():
        key = dict(key)
        pos = key['pos']
        mut = key['mut']

        if apply_cutoff_config(cutoff_config, value):
            usual_mut.append((pos, mut))

    return usual_mut


def prepare_cutoff_config(prevalence):
    cutoff_config = {
    }

    # num_genotype = len(set(
    #     i['genotype']
    #     for i in prevalence
    # ))

    # num_cut = 2

    # for i in range(num_genotype):
    #     num_genotype_cut = num_genotype - i
    #     if num_genotype_cut == 0:
    #         break

    #     cutoff_config[num_cut] = num_genotype_cut
    #     num_cut *= 2

    # print(cutoff_config)

    genotype_limit = {
        i['genotype']: max(ceil(i['total'] * 0.005), 2)
        for i in prevalence
    }

    return genotype_limit


def apply_cutoff_config(cut_off_config, value):
    # for num_cut, num_genotype_cut in cut_off_config.items():

    #     num_at_least = len([
    #         i
    #         for i in value
    #         if i['num'] >= num_cut
    #     ])

    #     if num_at_least >= num_genotype_cut:
    #         return True

    for i in range(len(cut_off_config)):
        at_least = len(cut_off_config) - i

        num_at_least = len([
            r
            for r in value
            if r['num'] >= (cut_off_config[r['genotype']] * (i + 1))
        ])

        if num_at_least >= at_least:
            return True

    return False


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

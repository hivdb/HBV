from .preset import DB
from preset.file_format import dump_csv
from preset.table import group_records_by
from operator import itemgetter
from scipy.stats import entropy
from math import ceil
from .preset import usual_cutoff_func
from .preset import get_usual_mutation_func
import re
from .preset import AA
from Bio.Align import substitution_matrices


def dump_pos_mut_by_genotype(
        prevalence, save_folder, exclude_genotype=['RF']):

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
        (i['pos'], i['mut'], i['overall_cons'])
        for i in prevalence
    ])

    pos_mut_list = {
        (pos, mut, cons): {
            g: (0, 0)
            for g in genotypes
        }
        for pos, mut, cons in pos_mut_list
    }

    for i in prevalence:
        pos_mut_list[
            (
                i['pos'], i['mut'], i['overall_cons']
            )][i['genotype']] = (
                i['pcnt'], i['num'])

    get_usual_mutation = globals()[get_usual_mutation_func]

    usual_mut = get_usual_mutation(prevalence)

    report = []

    for (pos, mut, cons), value in pos_mut_list.items():

        record = {
            'pos': pos,
            'overall_cons': cons,
            'mut': mut,
        }

        # TODO assert function
        # asserter = []
        # overall_asserter = 0
        for g, (pcnt, num) in value.items():
            record[f"#{g} ({genotype_total[g]})"] = num
            record[f"%{g}"] = f'{pcnt}%'

            # if g != 'overall':
            #     asserter.append(num)
            # else:
            #     overall_asserter = num

        # assert overall_asserter >= sum(asserter), f"{overall_asserter}{asserter}"

        if (pos, mut) in usual_mut:
            record['is_usual'] = 'yes'

        record['blosum62'] = calc_amino_acid_substitution(
            mut, cons
        )
        # record['maybe_APOBEC'] = 'yes' if maybe_apobec(cons, mut) else ''

        report.append(record)

    report.sort(key=itemgetter('pos', 'mut'))

    dump_csv(DB / save_folder / 'genotype_compare.csv', report)


def get_usual_mutation_by_genotype_num(prevalence):

    prevalence = [
        i
        for i in prevalence
        if i['genotype'] != 'overall'
    ]

    prepare_cutoff_config, apply_cutoff_config = usual_cutoff_func
    prepare_cutoff_config = globals()[prepare_cutoff_config]
    apply_cutoff_config = globals()[apply_cutoff_config]

    cutoff_config = prepare_cutoff_config(prevalence)

    usual_mut = []

    for key, value in group_records_by(prevalence, ['pos', 'mut']).items():
        key = dict(key)
        pos = key['pos']
        mut = key['mut']

        if apply_cutoff_config(cutoff_config, value):
            usual_mut.append((pos, mut))

    return usual_mut


def gradually_by_genotype_and_num(prevalence):

    genotype_limit = {
        i['genotype']: max(ceil(i['total'] * 0.005), 2)
        for i in prevalence
    }

    return genotype_limit


def apply_gradually_by_genotype_and_num(cut_off_config, value):

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


def get_usual_mutation_by_overall_0001(prevalence):
    return [
        (i['pos'], i['mut'])
        for i in prevalence
        if (
            i['genotype'] == 'overall'
            and i['pcnt'] >= 0.1
            and re.match(AA, i['mut'])
        )
    ]


def dump_pos_mut_by_mutation(prevalence, save_folder, exclude_genotype=['RF']):

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

    dump_csv(DB / save_folder / 'genotype_compare_by_mut.csv', report)


blosum62 = substitution_matrices.load("BLOSUM62")
blosum80 = substitution_matrices.load('BLOSUM80')


def calc_amino_acid_substitution(a1, a2, special=['-', '_', 'X']):
    if (a1 in special) or (a2 in special):
        return ''
    return blosum62[a1][a2]

APOBEC = [
    ('A', 'T'),
    ('C', 'Y'),
    ('D', 'N'),
    ('E', 'K'),
    ('G', 'R'),
    ('G', 'K'),
    ('G', 'S'),
    ('G', 'E'),
    ('G', 'D'),
    ('G', 'N'),
    ('M', 'I'),
    ('R', 'K'),
    ('R', 'H'),
    ('R', 'Q'),
    ('S', 'N'),
    ('V', 'M'),
    ('V', 'I'),
    ('W', '*'),
]


def maybe_apobec(cons, mut):
    if (cons, mut) in APOBEC:
        return True

from pathlib import Path

WS = Path(__file__).resolve().parent.parent
DB = WS / 'database'

AA = r'[AC-IK-WY]'

overall_pcnt_round = 1
genotype_pcnt_round = 0

usual_cutoff_func = (
    'gradually_by_genotype_and_num', 'apply_gradually_by_genotype_and_num')

get_usual_mutation_func = 'get_usual_mutation_by_overall_0001'
# get_usual_mutation_func = 'get_usual_mutation_by_genotype_num'

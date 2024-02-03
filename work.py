from hbv.sequence import get_sequences
from hbv.prevalence import get_prevalence
from hbv.prevalence import seq_filter
from hbv.usual_muts import view_usual_muts_by_genotype
from hbv.inter_intra import get_inter_intra
from hbv.preset import DB
from hbv.sierra import prepare_sierra_files

get_sequences()
get_prevalence()
view_usual_muts_by_genotype()
# get_inter_intra()

get_prevalence(save_folder=DB/'prevalence_2', filter=seq_filter)
view_usual_muts_by_genotype(save_folder=DB/'prevalence_2')
prepare_sierra_files()

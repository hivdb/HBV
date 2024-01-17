from pathlib import Path
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def dump_fasta(file_path, seq_list, ref_name=None):
    # TODO: ref can be a name or a pair of name:seq

    file_path = Path(file_path)
    file_path.parent.mkdir(exist_ok=True, parents=True)

    reference = {}
    if ref_name:
        reference = {
            k: v
            for k, v in seq_list.items()
            if k == ref_name
        }
        seq_list = {
            k: v
            for k, v in seq_list.items()
            if k != ref_name
        }

    dump_list = [
        SeqRecord(
            id=str(k),
            seq=Seq(v),
            description=''
        )
        for k, v in reference.items()
    ]

    dump_list += [
        SeqRecord(
            id=str(k),
            seq=Seq(v),
            description=''
        )
        for k, v in seq_list.items()
    ]

    SeqIO.write(
        dump_list,
        str(file_path),
        'fasta')


def load_fasta(file_path, with_name=True):
    # TODO
    # return types
    # 1. simple seq
    # 2. seq list
    # 3. seq with name in dict
    # 4. list of {'name': a, 'seq': b}
    seq_list = {}
    for i in SeqIO.parse(str(file_path), 'fasta'):
        # TODO, name, id, desc
        seq_list[i.name] = str(i.seq)

    if with_name:
        return seq_list

    seq_list = list(seq_list.values())

    if len(seq_list) == 1:
        return seq_list[0]

    return seq_list

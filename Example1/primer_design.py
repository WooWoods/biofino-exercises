import sys

import primer3
from Bio import SeqIO


def load_seq(fseq):
    """载入上一步获取的基因序列"""
    handle = SeqIO.parse(fseq, 'fasta')
    return handle

def designer(seq):
    """传入一条序列，针对该序列设计引物"""
    # 引物设计的默认参数，如不符合需求可以做相应更改
    default_params = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_NUM_RETURN': 5,
        'PRIMER_OPT_SIZE': 25,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 28,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 5,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [[200, 400]],
        }
    # 参数更改示例
    # default_params['PRIMER_OPT_SIZE'] = 20
    try:
        result = primer3.bindings.designPrimers(
                {
                    'SEQUENCE_ID': seq.id,
                    'SEQUENCE_TEMPLATE': str(seq.seq)
                    },
                default_params
                )
    except Exception:
        result = None
    return result

def prepare_output(seq, result, fhandle):
    """将设计好的引物按一定的格式保存"""
    if result is None:
        fhandle.write(seq.id + '\n')
    else:
        tmp = [
                seq.id,
                result.get('PRIMER_LEFT_0_SEQUENCE'),
                result.get('PRIMER_RIGHT_0_SEQUENCE'),
                str(result.get('PRIMER_LEFT_0_TM')),
                str(result.get('PRIMER_RIGHT_0_TM')),
                str(result.get('PRIMER_LEFT_0')[0]),
                str(result.get('PRIMER_RIGHT_0')[0])
                ]
        fhandle.write('\t'.join(tmp) + '\n')

def main():
    handle = load_seq(fseq)
    fhandle = open('primers.txt', 'w')
    fhandle.write('seqID\tForward\tReverse\tTM_Forward\tTM_Reverse\tPosition_Forward\tPosition_Reverse\n')
    for record in handle:
        result = designer(record)
        prepare_output(record, result, fhandle)
    fhandle.close()



if __name__ == '__main__':
    fseq = sys.argv[1]
    main()
    



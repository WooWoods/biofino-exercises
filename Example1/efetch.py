import sys
import re

from Bio import Entrez, SeqIO


def load_HSP(fdb):
    """从refGene数据库中选出HSP家族的基因"""
    data = {}
    with open(fdb) as fh:
        for line in fh:
            if re.search(r'HSP', line):
                larr = line.split()
                gname = larr[-4]
                nm = larr[1]
                if gname not in data:
                    data[gname] = nm
    return data

def search_entrez(accession):
    """传入一个基因的accession，从NCBI上获取其序列"""
    # 设置为从API获取数据，以免高频率的访问对网站造成过重负担而被封禁ip
    Entrez.api_key = "MyAPIkey"
    Entrez.email = 'wujunjames@126.com'
    try:
        handle = Entrez.efetch(db='nucleotide', id=accession, rettype='gb', retmode='text')
        record = SeqIO.read(handle, 'genbank')
    except Exception:
        record = None
    return record

def main():
    data = load_HSP(fdb)
    with open('HSPfam.fa', 'w') as fh:
        for gname, acc in data.items():
            record = search_entrez(acc)
            # 逐个获取序列并保存为fasta格式
            if record is not None:
                SeqIO.write(record, fh, 'fasta')


if __name__ == '__main__':
    fdb = sys.argv[1]
    main()

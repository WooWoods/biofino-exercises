import sys
import re

from Bio import Entrez, SeqIO


def load_HSP(fdb):
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
            if record is not None:
                SeqIO.write(record, fh, 'fasta')


if __name__ == '__main__':
    fdb = sys.argv[1]
    main()

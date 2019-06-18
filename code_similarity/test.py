import os
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from sklearn.metrics import roc_auc_score, matthews_corrcoef
import pandas as pd

path = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\bindPredict_performance.txt'
data = pd.read_csv(path, sep=' ', header=0, engine='python')


prod = ((data['tp'] + data['fp']) * (data['tp'] + data['fn']) * (data['tn'] + data['fp']) * (data['tn'] + data['fn']))
mcc = (data['tp'] * data['tn'] - data['fp'] * data['fn']) / prod ** (0.5) #if prod != 0 else 0
data['mcc'] = mcc*100

print(data.mean())


def multiple_alignment(sequences, path, group_id, clustalw_command):
    # seqs = [Seq(x) for x in sequences if type(x) != Seq]
    records = [SeqRecord(sequences[x], id=str(x)) for x in range(0, len(sequences))]

    SeqIO.write(records, path + "\\" + group_id + ".fasta", "fasta")
    #	clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
    clustalw_cline = ClustalwCommandline(clustalw_command, infile=path + "\\" + group_id + ".fasta",
                                         outfile=path + "\\" + group_id + ".align")

    clustalw_cline()
    try:
        align = AlignIO.read(path + "\\" + group_id + ".align", "clustal")
    except ValueError as e:
        print(e)
        print(path + "\\" + group_id + ".align")
    return (align)


# seqs = [Seq('TRMRMRPWLEMQI'), Seq('TRMRMRPWLEMQ'), Seq('TRMRMRPWMQI'), Seq('MRMRPWLEMQI')]
# path = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\funfam_project\data'
# group_id = '1.1.1.1'
# clustalw_command = r'C:\Program Files (x86)\ClustalW2\clustalw2.exe'
#
# #align = multiple_alignment(seqs, path, group_id, clustalw_command)
# #clustalw_cline = ClustalwCommandline(clustalw_command, infile=path + '\\' + group_id + '.fasta')
#
# #clustalw_cline()
#
# name = 'abcdef'
# path = 'abcd/efg'
#
# print(os.path.join(path, name + '.aln'))
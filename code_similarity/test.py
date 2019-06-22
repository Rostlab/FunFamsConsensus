import os
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from sklearn.metrics import precision_score, accuracy_score, recall_score, f1_score, roc_auc_score, matthews_corrcoef
import pandas as pd
#
path = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\bindPredict_performance.txt'
data = pd.read_csv(path, sep=' ', header=0, engine='python')


prod = ((data['tp'] + data['fp']) * (data['tp'] + data['fn']) * (data['tn'] + data['fp']) * (data['tn'] + data['fn']))
mcc = (data['tp'] * data['tn'] - data['fp'] * data['fn']) / prod ** (0.5) #if prod != 0 else 0
data['mcc'] = mcc*100

print(data.mean())


def eval_cm(tp, fp, fn, tn, type):
    prec = tp / (tp + fp) if (tp + fp) != 0 else 1
    cov = tp / (tp + fn) if (tp + fn) != 0 else 1
    F1 = 2 * (cov * prec) / (cov + prec) if (cov + prec) != 0 else 0
    acc = (tp + tn) / (tp + tn + fp + fn)
    prod = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if prod == 0:
        mcc = np.NaN  # no annotated and/or no predicted binding sites
    else:
        mcc = (tp * tn - fp * fn) / prod ** (0.5)

    if type == 'prec':
        return prec
    elif type == 'cov':
        return cov
    elif type == 'F1':
        return F1
    elif type == 'acc':
        return acc
    elif type == 'mcc':
        return mcc


def compute_eval(predictions, annotation):
    trues = sum(predictions)
    falses = sum((predictions == False))
    tp = sum(predictions & annotation)
    fp = trues - tp
    fn = sum((predictions == False) & annotation)
    tn = falses - fn

    prec = tp / trues if trues != 0 else 1
    cov = tp / (tp + fn) if (tp + fn) != 0 else 1
    F1 = 2 * (cov * prec) / (cov + prec) if (cov + prec) != 0 else 0
    acc = (tp + tn) / (tp + tn + fp + fn)

    prod = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if prod == 0:
        mcc = np.NaN  # no annotated and/or no predicted binding sites
    else:
        mcc = (tp * tn - fp * fn) / prod ** (0.5)  # if prod != 0 else 0
        # mcc = matthews_corrcoef(annotation, predictions)

    return [prec, cov, F1, acc, mcc]

path = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\confusion_matrices_001.csv'
data = pd.read_csv(path, sep=',', header=0, engine='python')
print(data.shape)

for measure in ['prec','cov','F1','acc','mcc']:
    data[measure+'_ccs'] = np.vectorize(eval_cm)(data['tp_ccs'], data['fp_ccs'], data['fn_ccs'], data['tn_ccs'], measure)
for measure in ['prec', 'cov', 'F1', 'acc', 'mcc']:
    data[measure + '_cc'] = np.vectorize(eval_cm)(data['tp_cc'], data['fp_cc'], data['fn_cc'], data['tn_cc'],
                                                   measure)
for measure in ['prec','cov','F1','acc','mcc']:
    data[measure+'_ccs_cons'] = np.vectorize(eval_cm)(data['tp_ccs_cons'], data['fp_ccs_cons'], data['fn_ccs_cons'], data['tn_ccs_cons'], measure)
for measure in ['prec', 'cov', 'F1', 'acc', 'mcc']:
    data[measure + '_cc_cons'] = np.vectorize(eval_cm)(data['tp_cc_cons'], data['fp_cc_cons'], data['fn_cc_cons'], data['tn_cc_cons'],
                                                   measure)

#print(data.shape)
#counts = data['funfam'].value_counts()
#to_drop = data[data['funfam'].value_counts() == 1]
##data = data.drop(to_drop)
#group_size = data.groupby('funfam').size()
#print('groups < 2',group_size[group_size == 2])
#print(len(data['funfam']), len(set(data['funfam'])))
#data = data.groupby(['funfam']).filter(lambda x: group_size[x.funfam] >= 1)
data = data.groupby('funfam').mean()
print(data.shape)
#print(data.head())
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
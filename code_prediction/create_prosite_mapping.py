from Bio import SeqIO
from collections import defaultdict
import os


def read_prosite_files(path):
    uniprot_prosite_map = defaultdict(list)
    msa_files = os.listdir(
        r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\prosite_alignments')
    for file in msa_files:
        # print("reading", file)
        prosite_id = file[:-4]
        for record in SeqIO.parse(path + file, "fasta"):
            uniprot_id = record.id.split('|')[1].split('/')[0]
            start, end = record.id.split('/')[1][:-1].split('-')
            uniprot_prosite_map[uniprot_id].append((prosite_id, start, end))
    return uniprot_prosite_map


def read_prosite_mapping(path):
    uniprot_prosite_map = defaultdict(list)
    with open(path, 'r') as f:
        for line in f:
            line_split = line.split(',')
            prosite_id = line_split[0]
            uniprot_id = line_split[1]
            start, end = line_split[2:]
            uniprot_prosite_map[uniprot_id].append((prosite_id, start, end))

    return uniprot_prosite_map


uniprot_prosite_map = read_prosite_mapping(
    r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\prosite_alignments\prosite_mapping.txt')
print(len(uniprot_prosite_map))

# print("reading prosite files...")
# uniprot_prosite_map = read_prosite_files(
#     r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\prosite_alignments\\')
#
# print("writing mapping...")
# with open(
#         r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\prosite_alignments\prosite_mapping.txt',
#         'w') as f:
#     for key, values in uniprot_prosite_map.items():
#         for value in values:
#             f.write(value[0] + ',' + key + ',' + str(value[1]) + ',' + str(value[2]) + '\n')
# print("done")

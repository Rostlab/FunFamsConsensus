from collections import defaultdict
from itertools import chain

def read_used_entries(file):
    entries = set()
    with open(file, 'r') as f:
        for line in f:
            line = line[:-1] #remove line-break
            line_split = line.split(',')
            entries.add(tuple(line_split))
    print("number of entries read from file:",len(entries), "example entry:\n",list(entries)[0])
    return entries

def read_uniprot_pfam_mapping(file):
    uniprot_pfam_mapping = defaultdict(set)

    with open(file, "r") as f:
        for line in f:
            line = line.strip()
            line_split = line.split("\t")
            uniprot_id = line_split[0]
            start, end = map(int, line_split[1].split(','))
            pfam_id = line_split[2]

            uniprot_pfam_mapping[uniprot_id].add((start, end, pfam_id))

    return uniprot_pfam_mapping

def read_uniprot_prosite_mapping(path):
    uniprot_prosite_map = defaultdict(list)

    with open(path, 'r') as f:
        for line in f:
            line_split = line.split(',')
            prosite_id = line_split[0]
            uniprot_id = line_split[1]
            start, end = line_split[2:]
            uniprot_prosite_map[uniprot_id].append((prosite_id, start, end))

    return uniprot_prosite_map

def check_sizes(mapping):
    to_remove  = set()
    for k,v in mapping.items():
        if len(v) < 2:
            to_remove.update(v)
    return to_remove

path_funfam_entries = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\used_entries_funfam.txt'
path_pfam_entries = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\used_entries_pfam.txt'
path_prosite_entries = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\used_entries_prosite.txt'
path_ec_entries = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\used_entries_ec.txt'
path_uniprot_pfam = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\uniprot_pfam_map'
path_uniprot_prosite = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\prosite_mapping.txt'
path_out_pfam = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\pfam_subset.txt'

funfam_entries = read_used_entries(path_funfam_entries)
pfam_entries = read_used_entries(path_pfam_entries)
uniprot_pfam = read_uniprot_pfam_mapping(path_uniprot_pfam)
uniprot_prosite = read_uniprot_prosite_mapping(path_uniprot_prosite)



# remove_from_funfam_entries = set()
# for entry in funfam_entries:
#     if (entry[0], entry[1], entry[2]) not in pfam_entries:
#         remove_from_funfam_entries.add(entry)
# funfam_entries.difference_update(remove_from_funfam_entries)
#
# remove_from_pfam_entries = set()
# for entry in pfam_entries:
#     remove = True
#     for entry_f in funfam_entries:
#         if (entry_f[0], entry_f[1], entry_f[2]) == entry:
#             remove_from_pfam_entries.add(entry_f)
#             break
# pfam_entries.difference_update(remove_from_pfam_entries)

# for entry in pfam_entries.difference(funfam_entries):
#     funfam_entries.remove(entry)
#     pfam_entries.remove(entry)

print("constructing mappings...")
seq_to_funfam = defaultdict(list)
(seq_to_funfam[x[2]].append((x[0],[1])) for x in funfam_entries)

seq_to_uniprot = defaultdict(list)
uniprot_to_seq = defaultdict(list)
for superfamily, funfam, seq_id, uniprot_id in funfam_entries:
    seq_to_uniprot[seq_id].append(uniprot_id)
    uniprot_to_seq[uniprot_id].append(seq_id)
#print(len(seq_to_uniprot), len(set(chain(*seq_to_uniprot.values()))))
#print(len(uniprot_to_seq), len(set(chain(*uniprot_to_seq.values()))))


new_pfam_entries = set()
for entry in pfam_entries:
    uniprot_ids = seq_to_uniprot[entry[2]]
    #print(entry[2],uniprot_ids)
    for uniprot_id in uniprot_ids:
        new_pfam_entries.add(entry+(uniprot_id,))

print(len(new_pfam_entries))
entries = funfam_entries.symmetric_difference(new_pfam_entries)
print("number of entries in both groups:", len(entries))

funfam_to_seq = defaultdict(list)
pfam_to_seq = defaultdict(list)
seq_to_pfam = defaultdict(list)

for entry in entries:
    superfamily, funfam, e_id, uniprot_id = entry
    pfam_ids = uniprot_pfam[uniprot_id]
    for start, end, pfam_id in pfam_ids:
        seq_to_pfam[e_id].append(pfam_id)
        pfam_to_seq[pfam_id].append(e_id)
    funfam_to_seq[(superfamily, funfam)].append(e_id)

print("number of groups funfam:", len(funfam_to_seq), len(set(chain(*funfam_to_seq.values()))))
print("number of groups pfam:", len(pfam_to_seq), len(set(chain(*pfam_to_seq.values()))))

i = 0
print("start pruning...", len(funfam_to_seq),len(set(chain(*funfam_to_seq.values()))),len(pfam_to_seq), len(set(chain(*pfam_to_seq.values()))))
while True:
    i += 1
    remove_due_to_funfam = check_sizes(funfam_to_seq)
    remove_due_to_pfam = check_sizes(pfam_to_seq)
    #print(remove_due_to_pfam)
    #print(remove_due_to_pfam)
    if False:#not remove_due_to_funfam and not remove_due_to_pfam:
        break
    else:
        for entry in remove_due_to_funfam:
            pfam_ids = seq_to_pfam[entry]
            for pfam_id in pfam_ids:
                try:
                    pfam_to_seq[pfam_id].remove(entry)
                    if not pfam_to_seq[pfam_id]:
                        del pfam_to_seq[pfam_id]
                except ValueError:
                    print(entry)
                    continue
            funfams = seq_to_funfam[entry]
            for superfamily, funfam in funfams:
                try:
                    funfam_to_seq[(superfamily, funfam)].remove(entry)
                    if not funfam_to_seq[(superfamily, funfam)]:
                        del funfam_to_seq[(superfamily, funfam)]
                except ValueError:
                    print(entry)
                    continue
        for entry in remove_due_to_pfam:
            pfam_ids = seq_to_pfam[entry]
            for pfam_id in pfam_ids:
                try:
                    pfam_to_seq[pfam_id].remove(entry)
                    if not pfam_to_seq[pfam_id]:
                        del pfam_to_seq[pfam_id]
                except ValueError:
                    print(entry)
                    continue
            funfams = seq_to_funfam[entry]
            for superfamily, funfam in funfams:
                try:
                    funfam_to_seq[(superfamily, funfam)].remove(entry)
                    if not funfam_to_seq[(superfamily, funfam)]:
                        del funfam_to_seq[(superfamily, funfam)]
                except ValueError:
                    print(entry)
                    continue

    if i == 500:
        break
   # print("iteration: ", i, len(funfam_to_seq), len(pfam_to_seq))

print("done")
print(len(funfam_to_seq),len(set(chain(*funfam_to_seq.values()))),len(pfam_to_seq), len(set(chain(*pfam_to_seq.values()))))
with open(path_out_pfam, 'w') as f:
    for k,v in funfam_to_seq.items():
        for e_id in v:
            f.write(k[0]+','+k[1]+','+e_id+'\n')



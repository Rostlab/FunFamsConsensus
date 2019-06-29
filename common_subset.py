from collections import defaultdict

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
#pfam_entries = read_used_entries(path_pfam_entries)
uniprot_pfam = read_uniprot_pfam_mapping(path_uniprot_pfam)
uniprot_prosite = read_uniprot_prosite_mapping(path_uniprot_prosite)

print("constructing mappings...")
seq_to_funfam = defaultdict(list)
(seq_to_funfam[x[2]].append((x[0],[1])) for x in funfam_entries)

seq_to_uniprot = defaultdict(list)
(seq_to_uniprot[x[2]].append(x[3]) for x in funfam_entries)
uniprot_to_seq = defaultdict(list)
(seq_to_uniprot[x[3]].append(x[2]) for x in funfam_entries)

funfam_to_seq = defaultdict(list)
(funfam_to_seq[x[1]].append(x[2]) for x in funfam_entries)

pfam_to_seq = defaultdict(list)
#(pfam_to_seq[v[2]].extend(uniprot_to_seq(k)) for k,v in uniprot_pfam.values())

seq_to_pfam = defaultdict(list)
# for x in pfam_entries:
#     uniprot_ids = seq_to_uniprot[x[2]]
#     for uniprot_id in uniprot_ids:
#         pfam_ids = uniprot_pfam[uniprot_id]
#         seq_to_pfam[x[2]].extend(pfam_ids)

for entry in funfam_entries:
    superfamily, funfam, e_id, uniprot_id = entry
    pfam_id = uniprot_pfam[uniprot_id][2]
    seq_to_pfam[e_id].append(pfam_id)
    pfam_to_seq[pfam_id].append(e_id)

i = 0
print("start pruning...")
while True:
    i += 1
    remove_due_to_funfam = check_sizes(funfam_to_seq)
    remove_due_to_pfam = check_sizes(pfam_to_seq)
    if not remove_due_to_funfam and not remove_due_to_pfam:
        break
    else:
        for entry in remove_due_to_funfam:
            pfam_ids = seq_to_pfam[entry]
            for pfam_id in pfam_ids:
                pfam_to_seq[pfam_id].remove(entry)
                if not pfam_to_seq[pfam_id]:
                    del pfam_to_seq[pfam_id]
        for entry in remove_due_to_pfam:
            funfams = seq_to_funfam[entry]
            for superfamily, funfam in funfams:
                funfam_to_seq[(superfamily, funfam)].remove(entry)
                if not funfam_to_seq[(superfamily, funfam)]:
                    del funfam_to_seq[(superfamily, funfam)]

    print("finished iteration: ", i)

print("done")
with open(path_out_pfam, 'w') as f:
    for k,v in funfam_to_seq:
        f.write(k[0]+','+k[1]+','+v+'\n')



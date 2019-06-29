"""this script provides utility functions for:
- processing FunFamEntries w.r.t. their binding site annotation
- computing similarity scores
- creating a sequence alignment

"""
import os
import numpy as np
import math
from itertools import chain
from random import shuffle
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline

def standard_error(x):
    sd = np.std(x, ddof=1)
    se = sd / math.sqrt(len(x))
    return (se)

def process_funfam_entries(funfam_entries, uniprot_binding_site_mapping):
    relevant_entries = list()
    for entry in funfam_entries:
        entry.add_binding_sites(uniprot_binding_site_mapping)
        has_sites = entry.process_sites()
        if has_sites and entry.sites is not None and entry.sites:
            relevant_entries.append(entry)
    return relevant_entries


def multiple_similarity(sites):
    """
    returns an aggregate score for all binding sites in list of lists sites
    """
    scores = []
    for i in range(0, len(sites)):
        for j in range(i + 1, len(sites)):
            scores.append(pairwise_similarity(sites[i], sites[j]))
    if len(scores) == 0:
        return 0.0
        # raise ValueError("cannot compute similarity of zero sites")
    return (sum(scores) / len(scores))


def pairwise_similarity(sitesA, sitesB):
    """
    returns a simple similarity score for two binding site annotations
    """
    if len(sitesA) < 1 and len(sitesB) < 1:
        return 0.0
        # raise ValueError("cannot compute similarity of zero sites")
    sitesA = set(sitesA)
    sitesB = set(sitesB)
    intersection = sitesA.intersection(sitesB)
    size = len(sitesA) if len(sitesA) >= len(sitesB) else len(sitesB)
    return (len(intersection) / size)


def similarity(group, entries, grouping_keyword, limit_keyword, alignment_path, clustalw_command):
    used_ids = []

    if grouping_keyword == 'ec':
        if "-" in group or len(group.split(".")) != 4:
            return None
        sequences = [entry.sequence for entry in entries]
        alignment = multiple_alignment(sequences, alignment_path, 'ec_class_' + group, clustalw_command)
        alignment_dict = SeqIO.to_dict(alignment)

    if grouping_keyword == 'pfam':
        sequences = [entry.sequence for entry in entries]
        alignment = multiple_alignment(sequences, alignment_path, 'pfam_family_' + group, clustalw_command)
        alignment_dict = SeqIO.to_dict(alignment)

    if grouping_keyword == 'prosite':
        sequences = [entry.sequence for entry in entries]
        alignment = multiple_alignment(sequences, alignment_path, 'prosite_group_' + group, clustalw_command)
        alignment_dict = SeqIO.to_dict(alignment)

    binding_sites = []
    used_entries = []
    for i, entry in enumerate(entries):
        if entry.binding_site_id in used_ids:  # omit entries from already used uniprot ids
            continue
        used_ids.append(entry.binding_site_id)
        used_entries.append((entry.superfamily, entry.funfam, entry.id, entry.binding_site_id))
        if grouping_keyword == 'ec':
            entry.aligned_sequence_ec = alignment_dict.get(str(i))
        if grouping_keyword == 'pfam':
            entry.aligned_sequence_pfam = alignment_dict.get(str(i))
        if grouping_keyword == 'prosite':
            entry.aligned_sequence_prosite = alignment_dict.get(str(i))

        entry.map_binding_sites(grouping_keyword)

        if grouping_keyword == 'funfam' or grouping_keyword in ['funfam-on-ec-subset', 'funfam-on-pfam-subset', 'funfam-on-prosite-subset']:
            binding_sites.append(entry.mapped_sites_funfam)
        elif grouping_keyword == 'ec':
            binding_sites.append(entry.mapped_sites_ec)
        elif grouping_keyword == 'pfam':
            binding_sites.append(entry.mapped_sites_pfam)
        elif grouping_keyword == 'prosite':
            binding_sites.append(entry.mapped_sites_prosite)

    return (used_entries, multiple_similarity(binding_sites))

def read_used_entries(file):
    entries = set()
    with open(file, 'r') as f:
        for line in f:
            line = line[:-1] #remove line-break
            line_split = line.split(',')
            entries.add(tuple(line_split))
    print("number of entries read from file:",len(entries), "example entry:\n",list(entries)[0])
    return entries

def multiple_alignment(sequences, path, group_id, clustalw_command):
    records = [SeqRecord(sequences[x], id=str(x)) for x in range(0, len(sequences))]

    fasta_file = os.path.join(path, group_id + '.fasta')
    SeqIO.write(records, fasta_file, "fasta")

    aln_file = os.path.join(path, group_id + '.aln')
    clustalw_cline = ClustalwCommandline(clustalw_command,
                                         infile=fasta_file)

    clustalw_cline()
    try:
        align = AlignIO.read(aln_file, "clustal")
    except ValueError as e:
        print(e)
        print(aln_file)
    return (align)


def get_group_mapping(funfam_entries, groupby, limit, pfam_file, prosite_file, file_entries_to_use=None):
    mapping = defaultdict(list)
    used_uniprot_ids = defaultdict(list)
    used_limit_ids = defaultdict(list)

    if limit is not None:
        shuffle(funfam_entries)

    if groupby == 'funfam':
        for entry in funfam_entries:
            if limit is None:
                mapping[(entry.superfamily, entry.funfam)].append(entry)
            elif limit == 'ec':
                if entry.ec_ids is None or not entry.ec_ids:
                    continue
                # if entry.binding_site_id in used_uniprot_ids[(entry.superfamily, entry.funfam)]:
                #     continue
                filtered_ecs = [ec_id for ec_id in entry.ec_ids if not '-' in ec_id]
                if not filtered_ecs:
                    continue
                if [x for x in filtered_ecs if x in used_limit_ids[(entry.superfamily, entry.funfam)]]:
                    continue
                # used_uniprot_ids[(entry.superfamily, entry.funfam)].append(entry.binding_site_id)
                used_limit_ids[(entry.superfamily, entry.funfam)] += filtered_ecs
                mapping[(entry.superfamily, entry.funfam)].append(entry)

    elif groupby == 'ec':
        for entry in funfam_entries:
            if entry.ec_ids is None:
                continue
            for ec_id in entry.ec_ids:
                if "-" in ec_id or len(ec_id.split(".")) != 4:
                    continue
                if limit is None:
                    mapping[ec_id].append(entry)
                elif limit == 'funfam':
                    # if entry.binding_site_id in used_uniprot_ids[ec_id]:
                    #     continue
                    if (entry.superfamily, entry.funfam) in used_limit_ids[ec_id]:
                        continue
                    # used_uniprot_ids[ec_id].append(entry.funfam)
                    used_limit_ids[ec_id].append((entry.superfamily, entry.funfam))
                    mapping[ec_id].append(entry)

    elif groupby == 'pfam':
        mapping = get_pfam_mapping(pfam_file, funfam_entries, mapping)
        # uniprot_pfam_mapping = read_uniprot_pfam_mapping(pfam_file)
        # for entry in funfam_entries:
        #     pfam_ids_to_add = set()
        #     pfam_set = uniprot_pfam_mapping.get(entry.binding_site_id)
        #     if pfam_set is None:
        #         continue
        #     for start,end,pfam_id in pfam_set:
        #         start = int(start)
        #         end = int(end)
        #         if start >= entry.start and end <= entry.end:
        #             pfam_ids_to_add.add(pfam_id)
        #     for pfam_id in pfam_ids_to_add:
        #         mapping[pfam_id].append(entry)

    elif groupby == 'prosite':
        uniprot_prosite_mapping = read_uniprot_prosite_mapping(prosite_file)
        for entry in funfam_entries:
            prosite_ids_to_add = set()
            prosite_set = uniprot_prosite_mapping.get(entry.binding_site_id)
            if prosite_set is None:
                continue
            for prosite_id,start,end in prosite_set:
                start = int(start)
                end = int(end)
                if start >= entry.start and end <= entry.end:
                    prosite_ids_to_add.add(prosite_id)
            for prosite_id in prosite_ids_to_add:
                mapping[prosite_id].append(entry)

    if groupby == 'funfam-on-pfam-subset':
        used_for_pfam = read_used_entries(file_entries_to_use)
        for entry in funfam_entries:
            if (entry.superfamily, entry.funfam, entry.id) in used_for_pfam:
                mapping[(entry.superfamily, entry.funfam)].append(entry)

    if groupby == 'funfam-on-prosite-subset':
        used_for_prosite = read_used_entries(file_entries_to_use)
        for entry in funfam_entries:
            if (entry.superfamily, entry.funfam, entry.id) in used_for_prosite:
                mapping[(entry.superfamily, entry.funfam)].append(entry)

    if groupby == 'funfam-on-ec-subset':
        used_for_ec = read_used_entries(file_entries_to_use)
        for entry in funfam_entries:
            if (entry.superfamily, entry.funfam, entry.id) in used_for_ec:
                mapping[(entry.superfamily, entry.funfam)].append(entry)

    return mapping

def get_entries_to_remove(mapping):
    entries_to_drop = set()
    for key, entries in mapping.items():
        if len(entries) < 2:
            entries_to_drop.update(entries)
    return entries_to_drop

def consolidate_group_mappings(group_mapping1, group_mapping2):
    #mapping group_id -> list of entries
    #group_id = funfam id, ec number, pfam id, ...
    #1)get list of common entries
    #2)delete all that are not common from gm1 and gm2
    #3)check if there are groups with <2 members left in both gm1 and gm2
    #4)if yes, delete these members/groups from the other gm as well and goto 3)

    entries1 = set(chain(*group_mapping1.values()))
    entries2 = set(chain(*group_mapping2.values()))
    common_entries = entries1.union(entries2)

    #2)delete non commons from gm1
    to_delete1 = set()
    for key,entries in group_mapping1.items():
        for entry in entries:
            if entry not in common_entries:
                to_delete1.add((key,entry))
    for key, entry in to_delete1:
        group_mapping1[key].remove(entry)

    #2)delete non commons from gm2
    to_delete2 = set()
    for key,entries in group_mapping2.items():
        for entry in entries:
            if entry not in common_entries:
                to_delete2.add((key,entry))
    for key, entry in to_delete2:
        group_mapping2[key].remove(entry)

    for i in range(0,500):
        remove1 = get_entries_to_remove(group_mapping1)
        for entry_to_remove in remove1:
            delete_from_mapping(group_mapping1, entry_to_remove)
            delete_from_mapping(group_mapping2, entry_to_remove)
        remove2 = get_entries_to_remove(group_mapping2)
        for entry_to_remove in remove2:
            delete_from_mapping(group_mapping1, entry_to_remove)
            delete_from_mapping(group_mapping2, entry_to_remove)

    return (group_mapping1, group_mapping2)

def delete_from_mapping(mapping, entry):
    to_remove = []
    for group, entries in mapping.items():
        for entry2 in entries:
            if entry2 == entry:
                to_remove.append((group, entry2))
                break
    if to_remove:
        mapping[to_remove[0][0]].remove(to_remove[0][1])


def get_pfam_mapping(pfam_file, funfam_entries, mapping):
    uniprot_pfam_mapping = read_uniprot_pfam_mapping(pfam_file)
    for entry in funfam_entries:
        pfam_ids_to_add = set()
        pfam_set = uniprot_pfam_mapping.get(entry.binding_site_id)
        if pfam_set is None:
            continue
        for start, end, pfam_id in pfam_set:
            start = int(start)
            end = int(end)
            if start >= entry.start and end <= entry.end:
                pfam_ids_to_add.add(pfam_id)
        for pfam_id in pfam_ids_to_add:
            mapping[pfam_id].append(entry)

    return mapping


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

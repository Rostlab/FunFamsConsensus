"""this script provides utility functions for:
- processing FunFamEntries w.r.t. their binding site annotation
- computing similarity scores
- creating a sequence alignment

"""
import os
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline


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
        print("XXXXXXXX this should not happen XXXXXXX")
        return (0.0)
    return (sum(scores) / len(scores))


def pairwise_similarity(sitesA, sitesB):
    """
    returns a simple similarity score for two binding site annotations
    """
    if len(sitesA) < 1 and len(sitesB) < 1:
        print("XXXXXXXX this should not happen pairwise XXXXXXX")
        return (0.0)
    sitesA = set(sitesA)
    sitesB = set(sitesB)
    intersection = sitesA.intersection(sitesB)
    size = len(sitesA) if len(sitesA) >= len(sitesB) else len(sitesB)
    return (len(intersection) / size)


def similarity(group, entries, grouping_keyword, alignment_path, clustalw_command):
    relevant_entries = []
    used_ids = []

    if len(entries) < 2:
        return None

    if grouping_keyword == 'ec':
        if "-" in group or len(group.split(".")) != 4:
            return None
        sequences = [entry.sequence for entry in entries]
        alignment = multiple_alignment(sequences, alignment_path, 'ec_class_' + group, clustalw_command)
        alignment_dict = SeqIO.to_dict(alignment)

    for entry in entries:
        if entry.binding_site_id in used_ids:
            continue
        relevant_entries.append(entry)
        used_ids.append(entry.binding_site_id)

    if len(relevant_entries) < 2:
        return None


    binding_sites = []
    for i, entry in enumerate(entries):
        if grouping_keyword == 'ec':
            if entry not in relevant_entries:
                continue
            entry.aligned_sequence_ec = alignment_dict.get(str(i))

        entry.map_binding_sites(grouping_keyword)

        if grouping_keyword == 'funfam':
            binding_sites.append(entry.mapped_sites_funfam)
        elif grouping_keyword == 'ec':
            binding_sites.append(entry.mapped_sites_ec)

    return multiple_similarity(binding_sites)


def multiple_alignment(sequences, path, group_id, clustalw_command):
    # seqs = [Seq(x) for x in sequences if type(x) != Seq]
    records = [SeqRecord(sequences[x], id=str(x)) for x in range(0, len(sequences))]

    fasta_file = os.path.join(path, group_id + '.fasta')
    SeqIO.write(records, fasta_file, "fasta")

    #	clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
    # path + "\\" + group_id + ".fasta"
    # path + "\\" + group_id + ".aln"

    aln_file = os.path.join(path, group_id + '.aln')
    clustalw_cline = ClustalwCommandline(clustalw_command,
                                         infile=fasta_file)  # , outfile=path + "\\" + group_id + ".align")

    clustalw_cline()
    try:
        align = AlignIO.read(aln_file, "clustal")
    except ValueError as e:
        print(e)
        print(aln_file)
    return (align)


def get_group_mapping(funfam_entries, groupby, limit):
    mapping = defaultdict(list)
    used_uniprot_ids = defaultdict(list)
    used_limit_ids = defaultdict(list)

    if groupby == 'funfam':
        for entry in funfam_entries:
            if limit == 'none':
                mapping[(entry.superfamily, entry.funfam)].append(entry)
            elif limit == 'ec':
                if entry.ec_ids is None or not entry.ec_ids:
                    continue
                if entry.binding_site_id in used_uniprot_ids[(entry.superfamily, entry.funfam)]:
                    continue
                filtered_ecs = [ec_id for ec_id in entry.ec_ids if not '-' in ec_id]
                if not filtered_ecs:
                    continue
                if [x for x in filtered_ecs if x in used_limit_ids[(entry.superfamily, entry.funfam)]]:
                    continue
                used_uniprot_ids[(entry.superfamily, entry.funfam)].append(entry.binding_site_id)
                used_limit_ids[(entry.superfamily, entry.funfam)] += filtered_ecs
                mapping[(entry.superfamily, entry.funfam)].append(entry)

    elif groupby == 'ec':
        for entry in funfam_entries:
            if entry.ec_ids is None:
                continue
            for ec_id in entry.ec_ids:
                if "-" in ec_id or len(ec_id.split(".")) != 4:
                    continue
                if limit == 'none':
                    mapping[ec_id].append(entry)
                elif limit == 'funfam':
                    if entry.binding_site_id in used_uniprot_ids[ec_id]:
                        continue
                    if (entry.superfamily, entry.funfam) in used_limit_ids[ec_id]:
                        continue
                    used_uniprot_ids[ec_id].append(entry.funfam)
                    used_limit_ids[ec_id].append((entry.superfamily, entry.funfam))
                    mapping[ec_id].append(entry)

    return mapping

import os
import math
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict
from funfam_project.code_prediction.funfam import FunFam
from funfam_project.code_prediction.protein import Protein


def get_sequence_funfam_mapping(file):
    sequence_funfam_mapping = dict()
    funfam_superfamily_mapping = dict()

    for seq_record in SeqIO.parse(file, 'fasta'):
        header = seq_record.id.split('|')
        source = header[0]
        superfamily = header[1]
        funfam = header[2]
        uniprot_id = header[3]
        start = header[4].split('-')[0]
        end = header[4].split('-')[1]

        sites = header[8].split(':')[1].split(',')
        # make sure that only non-empty binding site annotations are mapped
        sites = list(map(int, sites)) if sites[0] else sites

        sequence_funfam_mapping[(uniprot_id, int(start), int(end))] = [superfamily, funfam, sites, seq_record.seq]
        funfam_superfamily_mapping[funfam] = superfamily

    return (sequence_funfam_mapping, funfam_superfamily_mapping)


def uniprot_funfam160_mapping(file):
    mapping = defaultdict(list)
    with open(file) as f:
        for line in f:
            line = line.strip()
            funfam = line.split("\t")[0]
            u_ids = line.split("\t")[1].split(",")
            for u_id in u_ids:
                mapping[u_id].append(funfam)

    return (mapping)


def get_valid_ids(file):
    with open(file) as f:
        content = f.readlines()
    valid_ids = [x.strip() for x in content]
    return (valid_ids)


def check_alignment_quality(alignment_statistics, seqlen_multiple, cov_cutoff):
    alignment_stats = pd.read_csv(alignment_statistics)
    if alignment_stats['num_seqs'][0] < seqlen_multiple * alignment_stats['seqlen'][0] or alignment_stats['perc_cov'][
        0] < cov_cutoff:
        return False
    else:
        return True


def get_focus_segment(final_outcfg, uniprot_id):
    out = []
    with open(final_outcfg) as f:
        segment = False
        i = 0
        for line in f:
            if line.startswith("  - " + uniprot_id):
                segment = True
                continue
            if segment:
                out.append(int(line.split(" ")[-1]))
                i = i + 1
            if i == 2:
                return (out)


def get_focus_region(final_outcfg):
    with open(final_outcfg) as f:
        for line in f:
            if line.startswith("focus_sequence:"):
                start, end = map(int, line.split("/")[1].split("-"))
                return ([start, end])


def get_evcouplings_positions(final_outcfg):
    positions = []
    with open(final_outcfg) as f:
        for line in f:
            if line.startswith("  - - ") or line.startswith("    - "):
                positions.append(int(line.split(" ")[-1]))
    return (positions)


def get_di_positions(di_file):
    if os.stat(di_file).st_size == 0:
        return (None)
    data = pd.read_csv(di_file, sep=' ', header=None)
    #	except pd.errors.EmptyDataError as e:
    data.columns = ['posA', 'aaA', 'posB', 'aaB', 'scoreAB', 'scoreBA']
    actual_start = data['posA'].iloc[0]
    actual_end = data['posB'].iloc[-1]
    temp_a = data['posA'].copy(deep=True).drop_duplicates()
    temp_a[len(temp_a) + 1] = actual_end
    di_positions = temp_a
    return (di_positions)


def check_region_sequence_length(start, end, aligned_sequence):
    len_region = end - start + 1
    len_seq = len([x for x in aligned_sequence if x != '-'])
    if len_region != len_seq:
        return False
    else:
        return True


def add_protein_to_funfam(basepath, funfams, consensus_cutoff, sites, superfamily, uniprot_id, aligned_sequence, funfam,
                          start_region, end_region, start_segment, end_segment, evc_positions, di_positions,
                          binding_annotation, start_seq, end_seq, cum_coup_cutoff, clust_coeff_cutoff):
    if funfam not in funfams.keys():
        ff = FunFam(funfam, consensus_cutoff)
        ff.superfamily = superfamily
        funfams[funfam] = ff
    else:
        ff = funfams.get(funfam)

    protein = Protein(uniprot_id, aligned_sequence, funfam, start_region, end_region, start_segment, end_segment,
                      evc_positions, di_positions, binding_annotation, start_seq, end_seq, cum_coup_cutoff,
                      clust_coeff_cutoff)
    if binding_annotation:
        protein.binding_sites(sites)
    protein.cum_scores(os.path.join(basepath, uniprot_id, uniprot_id + ".cum_scores"))
    protein.cluster_coeff(os.path.join(basepath, uniprot_id, uniprot_id + ".cluster_coeff"))
    if protein.no_predictions:
        return False
    protein.alignment_stats(os.path.join(basepath, uniprot_id, uniprot_id + "_alignment_statistics.csv"))

    ff.add_member(protein)
    return True


def get_funfam_sequence_info(file, uniprot_id):
    for seq_record in SeqIO.parse(file, "fasta"):
        if uniprot_id in seq_record.description:
            try:
                start, end = map(int, seq_record.id.split("/")[1].split("-"))
            except ValueError as e:
                print(e)
                print(uniprot_id, file)
            sequence = seq_record.seq
            return (start, end, sequence)

def standard_error(x):
    sd = np.std(x, ddof=1)
    se = sd / math.sqrt(len(x))
    return (se)

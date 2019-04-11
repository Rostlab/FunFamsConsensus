import os
import sys
import math
import pandas as pd
import numpy as np
from funfam_project.code_prediction.protein import Protein
from funfam_project.code_prediction.funfam import FunFam
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict


###################################
# (start|end)_range: the part of the uniprot sequence which is part of the FunFam and for which EVcouplings was run (FunFam sequence)
#
# (start|end)_segment: the part of the uniprot sequence with >70% column coverage. get_scores.py produces predictions only for these residues
#
# to easily map predictions onto the uniprot sequence, predictions of value 0.0 are added to the predictions such that every residue of the FunFam sequence has exactly one value
###################################

def standard_error(x):
    sd = np.std(x, ddof=1)
    se = sd / math.sqrt(len(x))
    return (se)


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
    protein.cum_scores(basepath + uniprot_id + "/" + uniprot_id + ".cum_scores")
    protein.cluster_coeff(basepath + uniprot_id + "/" + uniprot_id + ".cluster_coeff")
    if protein.no_predictions:
        return (False)
    protein.alignment_stats(basepath + uniprot_id + "/" + uniprot_id + "_alignment_statistics.csv")

    ff.add_member(protein)
    return (True)


def check_alignment_quality(uniprot_id, alignment_statistics, seqlen_multiple, cov_cutoff):
    alignment_stats = pd.read_csv(alignment_statistics)
    if alignment_stats['num_seqs'][0] < seqlen_multiple * alignment_stats['seqlen'][0] or alignment_stats['perc_cov'][
        0] < cov_cutoff:
        # print("\tnum_seqs:",alignment_stats['num_seqs'][0])
        # print("\tseqlen:",alignment_stats['seqlen'][0])
        # print("\tperc_cov:",alignment_stats['perc_cov'][0])
        return (False)
    else:
        return (True)


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


def get_focus_segment_alt(final_outcfg):
    with open(final_outcfg) as f:
        for line in f:
            if line.startswith("focus_sequence:"):
                start, end = list(map(int, line.split("/")[1].split("-")))
                break
    return ((start, end))


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
    #	out = []
    with open(final_outcfg) as f:
        # region = False
        # i = 0
        # for line in f:
        # 	if line.startswith("region:"):
        # 		region = True
        # 		continue
        # 	if region:
        # 		out.append(int(line.split(" ")[1]))
        # 		i = i + 1
        # 	if i == 2:
        # 		return(out)
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


def get_valid_ids(file):
    with open(file) as f:
        content = f.readlines()
    valid_ids = [x.strip() for x in content]
    return (valid_ids)


def check_region_sequence_length(start, end, aligned_sequence):
    len_region = end - start + 1
    len_seq = len([x for x in aligned_sequence if x != '-'])
    if len_region != len_seq:
        return (False)
    else:
        return (True)


def get_sequence_funfam_mapping(file):
    sequence_funfam_mapping = dict()
    funfam_superfamily_mapping = dict()
    # to_print = []
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
    # to_print.append([superfamily,funfam])

    # with open("/mnt/project/funfams/evcouplings/score_calculations/funfams/funfams.txt",'w') as f:
    # 	for entry in to_print:
    # 		f.write(entry[0]+'\t'+entry[1]+'\n')
    return ((sequence_funfam_mapping, funfam_superfamily_mapping))


def main():
    used_uniprot_ids = set()
    used_funfams = set()
    consensus_cutoff = float(sys.argv[1])
    cum_cutoff = float(sys.argv[2])
    clust_cutoff = float(sys.argv[3])

    # print("consensus_cutoff:",consensus_cutoff)
    # print("cum_cutoff:",cum_cutoff)
    # print("clust_cutoff:",clust_cutoff)

    basepath = "/mnt/project/funfams/evcouplings/score_calculations/"
    funfams = dict()

    # proteins for which evcouplings was run
    valid_ids = get_valid_ids(basepath + "evc_done.txt")
    #	print("valid ids:",len(valid_ids))

    # get funfam, binding sites, etc. information about proteins
    sequence_funfam_mapping, funfam_superfamily_mapping = get_sequence_funfam_mapping(
        basepath + 'trimmed_dataset_160_with_sites.csv')

    uniprot_id_funfam160_mapping = uniprot_funfam160_mapping(
        "/mnt/project/funfams/evcouplings/funfam160_uniprot_mapping.txt")
    uniprot_id_funfam160_mapping_whole = uniprot_funfam160_mapping(
        "/mnt/project/funfams/evcouplings/funfam_uniprot_mapping_160.csv")

    w_o_binding_site_annotation = 0

    # create FunFam and Protein objects for all data available
    for uniprot_id in valid_ids:
        uniprot_id_without_region = uniprot_id.split("_")[0]
        # print(uniprot_id)
        # print("\t",uniprot_id_without_region)

        # skip protein if alignment quality is insufficient
        if not check_alignment_quality(uniprot_id,
                                       "/mnt/project/funfams/evcouplings/output/" + uniprot_id + "/align/" + uniprot_id + "_alignment_statistics.csv",
                                       3, 0.7):
            # print("omitting",uniprot_id,"due to poor alignment quality")
            continue

        outcfg_file = "/mnt/project/funfams/evcouplings/output/" + uniprot_id + "_final.outcfg"

        # get start and end of protein range for which evcouplings was run
        start_region, end_region = get_focus_region(outcfg_file)

        # get start and end of protein segment for which evcouplings produced output
        start_segment, end_segment = get_focus_segment(outcfg_file, uniprot_id_without_region)

        # get all positions (columns) with sufficient coverage for evcouplings
        evc_positions = get_evcouplings_positions(outcfg_file)

        di_file = "/mnt/project/funfams/evcouplings/freecontact/" + uniprot_id + ".di"

        # get all positions with di scores
        di_positions = get_di_positions(di_file)
        if di_positions is None:
            # print("omitting",uniprot_id,"due to empty di scores file")
            # print("\t",di_file)
            continue

        funfams_list = uniprot_id_funfam160_mapping.get(uniprot_id_without_region)
        if funfams_list is None:
            funfams_list = uniprot_id_funfam160_mapping_whole.get(uniprot_id_without_region)
        # print("funfam from funfam160 mapping is None",uniprot_id,funfams_list)

        # get FunFam annotation for given uniprot_id + segment
        entry_info = sequence_funfam_mapping.get((uniprot_id_without_region, start_region, end_region))

        # protein is in trimmed_dataset_160_with_sites.csv and has binding site annotation
        if entry_info is not None:
            superfamily, funfam, sites, aligned_sequence = entry_info
            binding_annotation = True
            start_seq, end_seq = (start_region, end_region)
            if not check_region_sequence_length(start_region, end_region, aligned_sequence):
                # print("omitting",uniprot_id,"due to differing sequence and region length")
                # print("\tstart,end sequence:",start_seq,end_seq,end_seq-start_seq+1,"length sequence:",len([x for x in aligned_sequence if x != '-']))
                continue
            if funfam not in funfams_list:
                pass
            # print("error!",uniprot_id,"funfams from mapping != funfam from trimmed_dataset_160_with_sites",funfams_list,funfam)

            if not add_protein_to_funfam(basepath, funfams, consensus_cutoff, sites, superfamily, uniprot_id,
                                         aligned_sequence, funfam, start_region, end_region, start_segment, end_segment,
                                         evc_positions, di_positions, binding_annotation, start_seq, end_seq,
                                         cum_cutoff, clust_cutoff):
                # print("omitting",uniprot_id,"due to differing sequence and score length")
                continue
            used_uniprot_ids.add(uniprot_id_without_region)
            used_funfams.add(funfam)
        # protein is not part of the dataset with binding site annotations, use this protein only for consensus prediction
        else:
            # print("no binding site annotation for",uniprot_id)
            for funfam in funfams_list:
                superfamily = funfam_superfamily_mapping.get(funfam)
                #	print("\t"+uniprot_id,funfam,superfamily)
                start_seq, end_seq, aligned_sequence = get_funfam_sequence_info(
                    "/mnt/project/funfams/families/" + superfamily + "/" + funfam + ".core.nofrag.aln", uniprot_id)
                binding_annotation = False
                sites = None
                if not check_region_sequence_length(start_seq, end_seq, aligned_sequence):
                    # print("\tomitting",uniprot_id,"from",superfamily,funfam,"due to differing sequence and region length")
                    # print("\t",end_region,start_region,aligned_sequence)
                    continue

                if not start_region <= start_seq or not end_region >= end_seq:
                    # print("\tomitting, evcouplings segment is smaller than FunFam sequence region")
                    continue

                if not add_protein_to_funfam(basepath, funfams, consensus_cutoff, sites, superfamily, uniprot_id,
                                             aligned_sequence, funfam, start_region, end_region, start_segment,
                                             end_segment, evc_positions, di_positions, binding_annotation, start_seq,
                                             end_seq, cum_cutoff, clust_cutoff):
                    # print("omitting",uniprot_id,"due to differing sequence and score length")
                    continue
                used_uniprot_ids.add(uniprot_id_without_region)
                used_funfams.add(funfam)
            w_o_binding_site_annotation += 1

    i = 1
    j = 0
    eval_per_seq_entries = 0
    index = list(range(1, len(funfams.keys()) + 1))
    columns = ['prec_cum', 'cov_cum', 'F1_cum', 'prec_clust', 'cov_clust', 'F1_clust']
    evaluation_means = pd.DataFrame(index=index, columns=columns)
    evaluation_means_no_cons = pd.DataFrame(index=index, columns=columns)
    evaluation_means_consensus_annotation = pd.DataFrame(index=index, columns=["precision", "coverage", "F1"])
    num_predicted_binding_sites_cum = []
    num_predicted_binding_sites_clust = []
    fract_correct_cum = []
    fract_correct_clust = []
    no_correct_prediction_cum = []
    no_correct_prediction_clust = []
    len_alignment = []

    proteins = set()

    length = sum([len(funfam.members) for funfam in funfams.values()])
    # print("number of FunFam members:",length)

    evaluation_full = pd.DataFrame(index=index,
                                   columns=['FunFam', 'members', 'F1_cum', 'F1_cum_cons', 'F1_clust', 'F1_clust_cons',
                                            'Prec_cum', 'Prec_cum_const', 'Prec_clust', 'Prec_clust_cons', 'Cov_cum',
                                            'Cov_cum_cons', 'Cov_clust', 'Cov_clust_cons'])
    evaluation_per_seq = pd.DataFrame(index=list(range(1, length + 1)),
                                      columns=["FunFam", "protein", 'F1_cum_base', 'F1_clust_base', 'prec_cum',
                                               'cov_cum', 'F1_cum', 'prec_clust', 'cov_clust', 'F1_clust'])
    evaluation_new_consensus = pd.DataFrame(index=index, columns=['FunFam', 'members', 'prec_cum_base', 'cov_cum_base',
                                                                  'F1_cum_base', 'prec_cum_min', 'cov_cum_min',
                                                                  'F1_cum_min', 'prec_cum_mean', 'cov_cum_mean',
                                                                  'F1_cum_mean', 'prec_cum_max', 'cov_cum_max',
                                                                  'F1_cum_max', 'prec_clust_base', 'cov_clust_base',
                                                                  'F1_clust_base', 'prec_clust_min', 'cov_clust_min',
                                                                  'F1_clust_min', 'prec_clust_mean', 'cov_clust_mean',
                                                                  'F1_clust_mean', 'prec_clust_max', 'cov_clust_max',
                                                                  'F1_clust_max'])
    evaluation_consensus_annotation = pd.DataFrame(index=index,
                                                   columns=["FunFam", "members", "prec_cons_annot", "cov_cons_annot",
                                                            "F1_cons_annot"])
    # iterate through FunFam objects to compute evaluation metrics
    for ff_id, funfam in funfams.items():
        if len(funfam.members) < 1:
            continue
        funfam.binding_sites()
        funfam.predictions_cluster_coeff()
        funfam.predictions_cum_scores()
        funfam.evaluation()

        funfam.build_new_consensus()
        funfam.predictions_for_new_consensus(cum_cutoff, clust_cutoff)
        funfam.evaluate_new_consensus()

        #################################################################################
        funfam.build_annotation_consensus()
        funfam.evaluate_annotation_consensus()

        # if not np.array_equal(funfam.binding_sites["consensus"],funfam.consensus_annotation):
        # 	print(funfam.name)
        # 	print(list(funfam.binding_sites["consensus"]))
        # 	print()
        # 	print(funfam.consensus_annotation)
        # 	print()
        #################################################################################

        proteins.update([member.id for member in funfam.members])

        if not os.path.exists(basepath + "funfams/" + ff_id + "/"):
            os.makedirs(basepath + "funfams/" + ff_id + "/")

        funfam.binding_sites.to_csv(basepath + "funfams/" + ff_id + "/binding_sites.tsv", sep='\t')
        funfam.predictions_cluster_coeff.to_csv(basepath + "funfams/" + ff_id + "/predictions_cluster_coeff.tsv",
                                                sep='\t')
        funfam.predictions_cum_scores.to_csv(basepath + "funfams/" + ff_id + "/predictions_cum_scores.tsv", sep='\t')

        # write evaluation file for the FunFam
        try:
            os.remove(basepath + "funfams/" + ff_id + "/evaluation.tsv")
        except OSError:
            pass
        with open(basepath + "funfams/" + ff_id + "/evaluation.tsv", 'a') as f:
            funfam.evaluation.round(5).to_csv(f, sep='\t')
            f.write('\nconsensus:\n')
            funfam.evaluation_consensus.round(5).to_csv(f, sep='\t')

        # for consensus:
        if funfam.num_binding_members > 1:  # new, previously included all sequences ############################
            evaluation_means.loc[i] = funfam.evaluation_consensus.mean()
            evaluation_means_no_cons.loc[i] = funfam.evaluation.mean()
            evaluation_means_consensus_annotation.loc[i] = funfam.annotation_eval
            num_predicted_binding_sites_cum.append(sum(funfam.predictions_cum_scores['consensus']))
            num_predicted_binding_sites_clust.append(sum(funfam.predictions_cluster_coeff['consensus']))
            len_alignment.append(len(funfam.predictions_cluster_coeff['consensus']))
            if funfam.consensus_has_prediction_cum:
                fract_correct_cum.append(funfam.predictions_fract_correct_cum)
                no_correct_prediction_cum.append(sum(funfam.members_with_no_correct_prediction_cum) / len(
                    funfam.members_with_no_correct_prediction_cum))
            if funfam.consensus_has_prediction_clust:
                fract_correct_clust.append(funfam.predictions_fract_correct_clust)
                no_correct_prediction_clust.append(sum(funfam.members_with_no_correct_prediction_clust) / len(
                    funfam.members_with_no_correct_prediction_clust))
        # for baseline:
        # evaluation_means.loc[i] = funfam.evaluation.mean()

        # for full eval:
        #		print(funfam.name,len(funfam.members),len(funfam.evaluation['F1_cum']))
        for k, member in enumerate(funfam.members):
            try:
                #				print("\t",member.id,member.funfam)
                evaluation_per_seq.loc[eval_per_seq_entries] = [funfam.name, member.id, member.evaluation_values[2],
                                                                member.evaluation_values[5],
                                                                *member.evaluation_consensus]
                eval_per_seq_entries += 1
            except AttributeError as e:
                continue
        evaluation_full.loc[i] = [funfam.name, funfam.num_binding_members, funfam.evaluation['F1_cum'].mean(),
                                  funfam.evaluation_consensus['F1_cum'].mean(), funfam.evaluation['F1_clust'].mean(),
                                  funfam.evaluation_consensus['F1_clust'].mean(), funfam.evaluation['prec_cum'].mean(),
                                  funfam.evaluation_consensus['prec_cum'].mean(),
                                  funfam.evaluation['prec_clust'].mean(),
                                  funfam.evaluation_consensus['prec_clust'].mean(), funfam.evaluation['cov_cum'].mean(),
                                  funfam.evaluation_consensus['cov_cum'].mean(), funfam.evaluation['cov_clust'].mean(),
                                  funfam.evaluation_consensus['cov_clust'].mean()]
        evaluation_new_consensus.loc[i] = [funfam.name, funfam.num_binding_members,
                                           funfam.evaluation['prec_cum'].mean(), funfam.evaluation['cov_cum'].mean(),
                                           funfam.evaluation['F1_cum'].mean(), *funfam.cum_min_eval,
                                           *funfam.cum_mean_eval, *funfam.cum_max_eval,
                                           funfam.evaluation['prec_clust'].mean(),
                                           funfam.evaluation['cov_clust'].mean(), funfam.evaluation['F1_clust'].mean(),
                                           *funfam.clust_min_eval, *funfam.clust_mean_eval, *funfam.clust_max_eval]
        evaluation_consensus_annotation.loc[i] = [funfam.name, funfam.num_binding_members, *funfam.annotation_eval]
        j += len(funfam.evaluation['F1_cum'])
        i += 1

    # print("iterations:",i,"used:",j)
    # print('\t'.join(map(str,evaluation_means_no_cons.mean())))
    ###########################################################
    print('\t'.join(map(str, evaluation_means.mean())))
    print(standard_error(evaluation_means["prec_cum"]), standard_error(evaluation_means["cov_cum"]),
          standard_error(evaluation_means["F1_cum"]), standard_error(evaluation_means["prec_clust"]),
          standard_error(evaluation_means["cov_clust"]), standard_error(evaluation_means["F1_clust"]))
    print("cumulative couplings:\t", "F1 mean:", evaluation_means["F1_cum"].mean(), "F1 se:",
          standard_error(evaluation_means["F1_cum"]))
    print("\tbaseline:\t", "F1 mean:", evaluation_new_consensus["F1_cum_base"].mean(), "F1 se:",
          standard_error(evaluation_new_consensus["F1_cum_base"]))
    print("\tcoverage cum:", evaluation_means["cov_cum"].mean(), "se", standard_error(evaluation_means["cov_cum"]))
    print("clustering coefficents:\t", "F1 mean:", evaluation_means["F1_clust"].mean(), "F1 se:",
          standard_error(evaluation_means["F1_clust"]))
    print("\tbaseline:\t", "F1 mean:", evaluation_new_consensus["F1_clust_base"].mean(), "F1 se:",
          standard_error(evaluation_new_consensus["F1_clust_base"]))
    print("\tcoverage clust:", evaluation_means["cov_clust"].mean(), "se",
          standard_error(evaluation_means["cov_clust"]))
    print("number of predicted binding residues cum:", sum(num_predicted_binding_sites_cum), "clust:",
          sum(num_predicted_binding_sites_clust))
    print("averages per family cov:", sum(num_predicted_binding_sites_cum) / len(num_predicted_binding_sites_cum),
          "clust:", sum(num_predicted_binding_sites_clust) / len(num_predicted_binding_sites_clust))
    print("average alignment length:", sum(len_alignment) / len(len_alignment))
    print("num families used:", len(num_predicted_binding_sites_cum))
    print(sum(num_predicted_binding_sites_cum) / len(num_predicted_binding_sites_cum),
          sum(num_predicted_binding_sites_clust) / len(num_predicted_binding_sites_clust))
    num_predicted_binding_sites_cum = np.array(num_predicted_binding_sites_cum)
    num_predicted_binding_sites_clust = np.array(num_predicted_binding_sites_clust)
    print("num families with more than zero predicted binding sites; cum:",
          len(num_predicted_binding_sites_cum[num_predicted_binding_sites_cum > 0]), num_predicted_binding_sites_cum,
          "\nclust:", len(num_predicted_binding_sites_clust[num_predicted_binding_sites_clust > 0]),
          num_predicted_binding_sites_clust)
    print("fraction of funfams with zero predictions; cum:",
          len(num_predicted_binding_sites_cum[num_predicted_binding_sites_cum == 0]) / len(
              num_predicted_binding_sites_cum), "\nclust:",
          len(num_predicted_binding_sites_clust[num_predicted_binding_sites_clust == 0]) / len(
              num_predicted_binding_sites_clust))
    print("fraction of correct predictions; cum:", sum(fract_correct_cum) / len(fract_correct_cum), "clust:",
          sum(fract_correct_clust) / len(fract_correct_clust))
    print("fraction of sequences with no correct prediction; cum:",
          sum(no_correct_prediction_cum) / len(no_correct_prediction_cum), "clust:",
          sum(no_correct_prediction_clust) / len(no_correct_prediction_clust))

    print("number of used uniprot ids:", len(used_uniprot_ids))
    with open('/mnt/project/funfams/evcouplings/score_calculations/score_eval/used_uniprot_ids.txt', 'w') as f:
        for u_id in used_uniprot_ids:
            f.write(u_id+'\n')
        f.close()
    print("number of used funfams:", len(used_funfams))
    with open('/mnt/project/funfams/evcouplings/score_calculations/score_eval/used_funfams.txt', 'w') as f:
        for u_id in used_funfams:
            f.write(u_id+'\n')
        f.close()
    # print("fract_correct_cum:",fract_correct_cum)
    # print("fract_correct_clust:",fract_correct_clust)
    ###########################################################
    # print("number of used FunFams:",str(i))
    # print("number of Uniprot Ids:",str(len(proteins)),"(from total of",str(len(valid_ids))+",","without annotation:",str(w_o_binding_site_annotation)+")")
    evaluation_full.to_csv("/mnt/project/funfams/evcouplings/score_calculations/score_eval/evaluation_full_new.tsv",
                           sep="\t", index=False)
    evaluation_per_seq.to_csv("/mnt/project/funfams/evcouplings/score_calculations/score_eval/evaluation_per_seq.tsv",
                              sep="\t", index=False)
    evaluation_new_consensus.to_csv(
        "/mnt/project/funfams/evcouplings/score_calculations/score_eval/new_consensus_files/evaluation_new_consensus_" + str(
            cum_cutoff) + "_" + str(clust_cutoff) + ".tsv", sep="\t", index=False)
    evaluation_consensus_annotation.to_csv(
        "/mnt/project/funfams/evcouplings/score_calculations/score_eval/evaluation_consensus_annotation.tsv", sep="\t",
        index=False)


# print("\t".join(map(str,evaluation_means_consensus_annotation.mean())))


if __name__ == "__main__":
    main()

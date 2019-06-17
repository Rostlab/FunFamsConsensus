import os
import argparse
import pandas as pd
import numpy as np
from code_prediction.prediction_utils import get_sequence_funfam_mapping, uniprot_funfam160_mapping, \
    get_valid_ids, check_alignment_quality, get_focus_region, get_focus_segment, get_evcouplings_positions, \
    get_di_positions, check_region_sequence_length, add_protein_to_funfam, get_funfam_sequence_info, standard_error


def main():
    usage_string = 'python prediction.py'
    parser = argparse.ArgumentParser(description=__doc__, usage=usage_string)
    parser.add_argument("-consensus", dest="cons_cutoff", help="consensus cutoff", required=True)
    parser.add_argument("-cc", dest="clust_cutoff", help="clustering coefficient cutoff")
    parser.add_argument("-ccs", dest="cum_cutoff",
                        help="cumulative coupling score cutoff")
    parser.add_argument("-uniprot_ids", dest="uniprot_ids",
                        help="list of UNIPROT IDs to be processed", required=True)
    parser.add_argument("-mapping", dest="funfam_uniprot_mapping",
                        help="file with a mapping from FunFam to UNIPROT ID")
    parser.add_argument("-evc_info", dest="evc_info_dir",
                        help="directory with output files from EVcouplings, i.e. alignment statistics and outcfg files for the given UNIPROT IDs.")
    parser.add_argument("-funfam_data", dest="funfams_with_sites",
                        help="a file in FASTA FunFam format including mapped binding sites")
    parser.add_argument("-families", dest="funfams",
                        help="directory with FunFam data, FunFams sorted into subdirectories by their superfamilies")
    parser.add_argument("-out", dest="output_dir",
                        help="directory to which output files will be written")

    args = parser.parse_args()
    print("[ARGUMENTS]")
    print("\tconsensus cutoff:", args.cons_cutoff)
    print("\tclustering coefficient cutoff:", args.clust_cutoff)
    print("\tcumulative coupling scores cutoff:", args.cum_cutoff)
    print("\tUNIPROT ID file:", args.uniprot_ids)
    print("\tmapping file:", args.funfam_uniprot_mapping)
    print("\tevc output directory:", args.evc_info_dir)
    print("\tfunfam data file:", args.funfams_with_sites)
    print("[ARGUMENTS END]")

    funfams = dict()

    valid_ids = get_valid_ids(args.uniprot_ids)

    # get funfam, binding sites, etc. information about proteins
    sequence_funfam_mapping, funfam_superfamily_mapping = get_sequence_funfam_mapping(args.funfams_with_sites)

    # get funfam <-> uniprot id mapping
    # !check if smaller file has to be read too!
    uniprot_id_funfam160_mapping = uniprot_funfam160_mapping(args.funfam_uniprot_mapping)

    # create FunFam and Protein objects for all data available
    for uniprot_id in valid_ids:
        uniprot_id_without_region = uniprot_id.split("_")[0]
        #in general:
        # alignment_stats_file = os.path.join(args.evc_info_dir, uniprot_id, uniprot_id + '_alignment_statistics.csv')
        #di_file = os.path.join(args.evc_info_dir, uniprot_id,  uniprot_id + ".di")
        #outcfg_file = os.path.join(args.evc_info_dir, uniprot_id,  uniprot_id + "_final.outcfg")
        #our local structure:
        alignment_stats_file = os.path.join(args.evc_info_dir, 'output', uniprot_id, 'align', uniprot_id + '_alignment_statistics.csv')
        di_file = os.path.join(args.evc_info_dir, 'freecontact', uniprot_id + ".di")
        outcfg_file = os.path.join(args.evc_info_dir, 'output', uniprot_id + "_final.outcfg")

        # skip protein if alignment quality is insufficient
        COVERAGE_CUTOFF = 0.7
        SEQLEN_MULTIPLE = 3
        if not check_alignment_quality(alignment_stats_file, SEQLEN_MULTIPLE, COVERAGE_CUTOFF):
            continue

        # get start and end of protein range for which evcouplings was run
        start_region, end_region = get_focus_region(outcfg_file)

        # get start and end of protein segment for which evcouplings produced output
        start_segment, end_segment = get_focus_segment(outcfg_file, uniprot_id_without_region)

        # get all positions (columns) with sufficient coverage for evcouplings
        evc_positions = get_evcouplings_positions(outcfg_file)

        # get all positions with di scores
        di_positions = get_di_positions(di_file)
        if di_positions is None:
            #   di scores file empty
            continue

        funfams_list = uniprot_id_funfam160_mapping.get(uniprot_id_without_region)

        # get FunFam annotation for given uniprot_id + segment
        entry_info = sequence_funfam_mapping.get((uniprot_id_without_region, start_region, end_region))

        # protein is in trimmed_dataset_160_with_sites.csv and has binding site annotation
        if entry_info is not None:
            superfamily, funfam, sites, aligned_sequence = entry_info
            binding_annotation = True
            start_seq, end_seq = (start_region, end_region)
            if not check_region_sequence_length(start_region, end_region, aligned_sequence):
                # sequence and region length differ, omit sequence
                continue

            if not add_protein_to_funfam(args.evc_info_dir, funfams, args.cons_cutoff, sites, superfamily, uniprot_id,
                                         aligned_sequence, funfam, start_region, end_region, start_segment, end_segment,
                                         evc_positions, di_positions, binding_annotation, start_seq, end_seq,
                                         args.cum_cutoff, args.clust_cutoff):
                # sequence and score lengt differ, omit sequence
                continue

        # protein is not part of the dataset with binding site annotations, use this protein only for consensus prediction
        else:
            for funfam in funfams_list:
                superfamily = funfam_superfamily_mapping.get(funfam)
                start_seq, end_seq, aligned_sequence = get_funfam_sequence_info(
                    os.path.join(args.funfams, superfamily, funfam + ".core.nofrag.aln"), uniprot_id)
                binding_annotation = False
                sites = None
                if not check_region_sequence_length(start_seq, end_seq, aligned_sequence):
                    # sequence and region length differ, omit sequence
                    continue

                if not start_region <= start_seq or not end_region >= end_seq:
                    # evcouplings segment is smaller than FunFam sequence region, omit sequence
                    continue

                if not add_protein_to_funfam(args.evc_info_dir, funfams, args.cons_cutoff, sites, superfamily,
                                             uniprot_id,
                                             aligned_sequence, funfam, start_region, end_region, start_segment,
                                             end_segment, evc_positions, di_positions, binding_annotation, start_seq,
                                             end_seq, args.cum_cutoff, args.clust_cutoff):
                    # sequence and score lengt differ, omit sequence
                    continue

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
        funfam.predictions_for_new_consensus(args.cum_cutoff, args.clust_cutoff)
        funfam.evaluate_new_consensus()

        funfam.build_annotation_consensus()
        funfam.evaluate_annotation_consensus()

        proteins.update([member.id for member in funfam.members])

        if not os.path.exists(os.path.join(args.output_dir, "funfams", ff_id)):
            os.makedirs(os.path.join(args.output_dir, "funfams", ff_id))

        funfam.binding_sites.to_csv(os.path.join(args.output_dir, "funfam", ff_id, "binding_sites.tsv"), sep='\t')
        funfam.predictions_cluster_coeff.to_csv(os.path.join(args.output_dir, "funfams", ff_id, "predictions_cluster_coeff.tsv"),
                                                sep='\t')
        funfam.predictions_cum_scores.to_csv(os.path.join(args.output_dir, "funfams", ff_id, "predictions_cum_scores.tsv"), sep='\t')

        # write evaluation file for the FunFam
        try:
            os.remove(os.path.join(args.output_dir, "funfams", ff_id, "evaluation.tsv"))
        except OSError:
            pass
        with open(os.path.join(args.output_dir, "funfams", ff_id, "evaluation.tsv"), 'a') as f:
            funfam.evaluation.round(5).to_csv(f, sep='\t')
            f.write('\nconsensus:\n')
            funfam.evaluation_consensus.round(5).to_csv(f, sep='\t')

        # for consensus:
        if funfam.num_binding_members > 1:  # previously included all sequences
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

        for k, member in enumerate(funfam.members):
            try:
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

    evaluation_full.to_csv(os.path.join(args.output_dir, 'evaluation_full_new.tsv'),
                           sep="\t", index=False)
    evaluation_per_seq.to_csv(os.path.join(args.output_dir, 'evaluation_per_seq.tsv'),
                              sep="\t", index=False)
    evaluation_new_consensus.to_csv(
        os.path.join(args.output_dir, 'evaluation_new_consensus_' + str(
            args.cum_cutoff) + "_" + str(args.clust_cutoff) + ".tsv"), sep="\t", index=False)
    evaluation_consensus_annotation.to_csv(
        os.path.join(args.output_dir, 'evaluation_consensus_annotation.tsv'), sep="\t",
        index=False)


if __name__ == "__main__":
    main()

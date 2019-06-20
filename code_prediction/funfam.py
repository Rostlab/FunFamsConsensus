import time
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from sklearn.metrics import roc_auc_score


# clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"

class FunFam:
    '''object that represents a functional protein family
    a collection of Protein objects and methods to evaluate their collective binding site predictions'''

    def __init__(self, name, consensus_cutoff):
        self.name = name
        self.members = []
        self.consensus_cutoff = consensus_cutoff

    def superfamily(self, superfamily):
        self.superfamily = superfamily

    def add_member(self, protein):
        self.members.append(protein)

    def add_members(self, proteins):
        self.members.extend(proteins)

    def binding_sites(self):
        '''binding sites from trimmed_dataset_160_with_sites.csv are already mapped to alignment'''
        index = list(range(1, len(self.members[0].aligned_sequence) + 1))
        columns = [member.id for member in self.members if member.binding_annotation]
        data = [[True if x in member.binding_sites else False for x in index] for member in self.members if
                member.binding_annotation]
        self.binding_sites = pd.DataFrame(index=index, columns=columns, data=list(map(list, zip(*data))))
        self.binding_sites.insert(0, 'consensus', np.where(
            self.binding_sites.sum(axis=1) / len(self.binding_sites.columns) >= self.consensus_cutoff, True, False))

    def build_new_consensus(self):
        clust_scores = np.zeros((len(self.members), len(self.members[0].aligned_sequence)))
        cum_scores = np.zeros((len(self.members), len(self.members[0].aligned_sequence)))

        for i, member in enumerate(self.members):
            for j, position in enumerate(self.members[0].aligned_sequence):
                clust_scores[i, j] = member.cluster_coeffs[j]
                cum_scores[i, j] = member.cum_scores[j]

        self.consensus_clust_min = clust_scores.min(axis=0)
        self.consensus_clust_mean = clust_scores.mean(axis=0)
        self.consensus_clust_max = clust_scores.max(axis=0)

        self.consensus_cum_min = cum_scores.min(axis=0)
        self.consensus_cum_mean = cum_scores.mean(axis=0)
        self.consensus_cum_max = cum_scores.max(axis=0)

    def build_annotation_consensus(self):

        binding_per_position = np.zeros(len(self.members[0].aligned_sequence))

        for j, position in enumerate(self.members[0].aligned_sequence):
            num_binding = 0
            num_members = 0
            for i, member in enumerate(self.members):
                if not member.binding_annotation:
                    continue
                if j + 1 in member.binding_sites:  # binding sites start at 1, enumeration at 0 => 0 is not part of annotation
                    num_binding += 1
                num_members += 1
            binding_per_position[j] = num_binding

        self.num_binding_members = num_members
        self.consensus_annotation = (binding_per_position / num_members) > self.consensus_cutoff

    def evaluate_annotation_consensus(self):

        all_annotation_eval = []

        for member in self.members:
            if not member.binding_annotation:
                continue

            annotation_eval = self.compute_eval(
                self.map_from_alignment_to_sequence(member.aligned_sequence, self.binding_sites["consensus"]),
                self.map_from_alignment_to_sequence(member.aligned_sequence, self.binding_sites[member.id]))

            all_annotation_eval.append(annotation_eval)

        self.annotation_eval = np.mean(np.array(all_annotation_eval), axis=0)

    def predictions_for_new_consensus(self, cum_cutoff, clust_cutoff):
        self.predictions_clust_min = self.consensus_clust_min > clust_cutoff
        self.predictions_clust_mean = self.consensus_clust_mean > clust_cutoff
        self.predictions_clust_max = self.consensus_clust_max > clust_cutoff

        # if self.name == "295":
        # 	print("clust prediction min:",self.predictions_clust_min)
        # 	print("cut-off:",clust_cutoff)

        self.predictions_cum_min = self.consensus_cum_min > cum_cutoff
        self.predictions_cum_mean = self.consensus_cum_mean > cum_cutoff
        self.predictions_cum_max = self.consensus_cum_max > cum_cutoff

        if len(self.predictions_clust_min) != len(self.members[0].aligned_sequence):
            print("WTF XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

    def evaluate_new_consensus(self):

        all_clust_min_eval = []
        all_clust_mean_eval = []
        all_clust_max_eval = []

        all_cum_min_eval = []
        all_cum_mean_eval = []
        all_cum_max_eval = []

        for member in self.members:
            if not member.binding_annotation:
                continue

            #			print("XXXXXXXXXXXXXXXXXXXXX",len(member.aligned_sequence),len(self.predictions_clust_min))
            #			print(member.no_predictions)
            clust_min_eval = self.compute_eval(self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                   pd.Series(self.predictions_clust_min,
                                                                                             index=range(1, len(
                                                                                                 self.predictions_clust_min) + 1))),
                                               self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                   self.binding_sites[member.id]))
            clust_mean_eval = self.compute_eval(self.map_from_alignment_to_sequence(member.aligned_sequence, pd.Series(
                self.predictions_clust_mean, index=range(1, len(self.predictions_clust_mean) + 1))),
                                                self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                    self.binding_sites[member.id]))
            clust_max_eval = self.compute_eval(self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                   pd.Series(self.predictions_clust_max,
                                                                                             index=range(1, len(
                                                                                                 self.predictions_clust_max) + 1))),
                                               self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                   self.binding_sites[member.id]))

            all_clust_min_eval.append(clust_min_eval)
            all_clust_mean_eval.append(clust_mean_eval)
            all_clust_max_eval.append(clust_max_eval)

            cum_min_eval = self.compute_eval(self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                 pd.Series(self.predictions_cum_min,
                                                                                           index=range(1, len(
                                                                                               self.predictions_cum_min) + 1))),
                                             self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                 self.binding_sites[member.id]))
            cum_mean_eval = self.compute_eval(self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                  pd.Series(self.predictions_cum_mean,
                                                                                            index=range(1, len(
                                                                                                self.predictions_cum_mean) + 1))),
                                              self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                  self.binding_sites[member.id]))
            cum_max_eval = self.compute_eval(self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                 pd.Series(self.predictions_cum_max,
                                                                                           index=range(1, len(
                                                                                               self.predictions_cum_max) + 1))),
                                             self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                 self.binding_sites[member.id]))

            all_cum_min_eval.append(cum_min_eval)
            all_cum_mean_eval.append(cum_mean_eval)
            all_cum_max_eval.append(cum_max_eval)

        self.clust_min_eval = np.mean(np.array(all_clust_min_eval), axis=0)
        self.clust_mean_eval = np.mean(np.array(all_clust_mean_eval), axis=0)
        self.clust_max_eval = np.mean(np.array(all_clust_max_eval), axis=0)

        self.cum_min_eval = np.mean(np.array(all_cum_min_eval), axis=0)
        self.cum_mean_eval = np.mean(np.array(all_cum_mean_eval), axis=0)
        self.cum_max_eval = np.mean(np.array(all_cum_max_eval), axis=0)

    def map_to_seq2(self, aligned_sequence, predictions):
        if len(aligned_sequence) != len(predictions):
            print("ERROR sequence length != #predictions", self.name)
        if len(aligned_sequence) == 0 and len(predictions) == 0:
            print(self.name, "empty sequence & predictions")
            return (np.array([]))
        out = []
        for i, aa in enumerate(aligned_sequence):
            if aa != '-':
                out.append(predictions[i])
        return (np.array(out))

    def predictions_cluster_coeff(self):
        index = list(range(1, len(self.members[0].cluster_coeff_predictions) + 1))
        columns = [member.id for member in self.members]
        data = [member.cluster_coeff_predictions for member in self.members]
        self.predictions_cluster_coeff = pd.DataFrame(index=index, columns=columns, data=list(map(list, zip(*data))))
        self.predictions_cluster_coeff.insert(0, 'consensus', np.where(self.predictions_cluster_coeff.sum(axis=1) / len(
            self.predictions_cluster_coeff.columns) >= self.consensus_cutoff, True, False))
        self.predictions_fract_correct_clust = self.compute_fract_correct('clust')

    def predictions_cum_scores(self):
        index = list(range(1, len(self.members[0].cum_scores_predictions) + 1))
        columns = [member.id for member in self.members]
        data = [member.cum_scores_predictions for member in self.members]
        self.predictions_cum_scores = pd.DataFrame(index=index, columns=columns, data=list(map(list, zip(*data))))
        self.predictions_cum_scores.insert(0, 'consensus', np.where(
            self.predictions_cum_scores.sum(axis=1) / len(self.predictions_cum_scores.columns) >= self.consensus_cutoff,
            True, False))
        self.predictions_fract_correct_cum = self.compute_fract_correct('cum')

    def compute_fract_correct(self, method):
        fraction_correct_per_member = []
        if method == 'cum':
            consensus = self.predictions_cum_scores['consensus']
        elif method == 'clust':
            consensus = self.predictions_cluster_coeff['consensus']
        else:
            raise ValueError("method has to be cum or clust")

        if sum(consensus) == 0:
            if method == "cum":
                print(self.name, "no consensus prediction cum")
                self.consensus_has_prediction_cum = False
            elif method == "clust":
                print(self.name, "no consensus prediction cum")
                self.consensus_has_prediction_clust = False
            return ()
        else:
            if method == "cum":
                self.consensus_has_prediction_cum = True
                self.members_with_no_correct_prediction_cum = []
            elif method == "clust":
                self.consensus_has_prediction_clust = True
                self.members_with_no_correct_prediction_clust = []
        #		print(self.name)

        for member in self.members:
            if not member.binding_annotation:
                continue
            consensus_mapped = self.map_from_alignment_to_sequence(member.aligned_sequence, consensus)
            annotation = self.binding_sites[member.id]
            annotation_mapped = self.map_from_alignment_to_sequence(member.aligned_sequence, annotation)
            if len(annotation_mapped) != len(consensus_mapped):
                print("ERROR lenght mismatch annotation and consensus!")

            # print("\t",member.id)
            #			print("\tconsensus","\t",list(consensus))
            #			print("\tannotation","\t",list(annotation))
            #			print(list(consensus & annotation))
            # print(sum(consensus_mapped & annotation_mapped),sum(consensus_mapped))
            tps = sum(consensus_mapped & annotation_mapped)
            fraction_correct_per_member.append(tps / sum(consensus_mapped))
            if tps == 0:
                if method == "cum":
                    self.members_with_no_correct_prediction_cum.append(1)
                elif method == "clust":
                    self.members_with_no_correct_prediction_clust.append(1)
            else:
                if method == "cum":
                    self.members_with_no_correct_prediction_cum.append(0)
                elif method == "clust":
                    self.members_with_no_correct_prediction_clust.append(0)
        return (sum(fraction_correct_per_member) / len(fraction_correct_per_member))

    def map_score_to_sequence(self, aligned_sequence, scores):
        """returns the scores (float) mapped to the non-aligned sequence"""
        if len(aligned_sequence) != len(scores):
            raise IndexError("length of aligned sequence and scores differ!")
        else:
            length = len([x for x in aligned_sequence if x != '-'])
            output = np.zeros(length, dtype=float)
            j = 0
            for i, aa in enumerate(aligned_sequence):
                if aa != '-':
                    output[j] = scores[i + 1]
                    j += 1

            return (output)

    def compute_mean_auroc(self):
        values = []
        for i,member in enumerate(self.members):
            if not member.binding_annotation:
                continue
            annotation = self.map_from_alignment_to_sequence(member.aligned_sequence, self.binding_sites[member.id]).astype(int)
            cum_scores = self.map_score_to_sequence(member.aligned_sequence, self.predictions_cum_scores[member.id])
            clust_scores = self.map_score_to_sequence(member.aligned_sequence, self.predictions_cluster_coeff[member.id])
            #print(annotation,'\n')
            #print(cum_scores, '\n')
            #print(clust_scores, '\n')

            if sum(annotation) == 0:
                print("no binding sites left")
                continue
            cum_auroc = roc_auc_score(annotation, cum_scores)
            clust_auroc = roc_auc_score(annotation, clust_scores)
            #print('\n', cum_auroc, clust_auroc, '\n')
            values.append([cum_auroc, clust_auroc])
        out = pd.DataFrame(columns=['cum', 'clust'], index=range(0,len(values)), data=values)
        #print(out.head())
        #print(out.mean())
        return out.mean()

    def compute_eval(self, predictions, annotation, p=False, m=None, confusion_matrix=None, member_id=None):
        trues = sum(predictions)
        falses = sum((predictions == False))
        tp = sum(predictions & annotation)
        fp = trues - tp
        fn = sum((predictions == False) & annotation)
        tn = falses - fn

        if confusion_matrix is not None:
            print(member_id)
            if member_id is not None:
                confusion_matrix.append([self.name, member_id, fp, tp, fn, tn])
            else:
                confusion_matrix.append([fp, tp, fn, tn])

        prec = tp / trues if trues != 0 else 1
        cov = tp / (tp + fn) if (tp + fn) != 0 else 1
        F1 = 2 * (cov * prec) / (cov + prec) if (cov + prec) != 0 else 0
        acc = (tp + tn) / (tp + tn + fp + fn)

        prod = ((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
        if prod == 0:
            mcc = np.NaN    #no annotated and/or no predicted binding sites
        else:
            mcc = (tp * tn - fp * fn) / prod**(0.5) #if prod != 0 else 0
            #mcc = matthews_corrcoef(annotation, predictions)
        #auroc = roc_auc_score(annotation, predictions)

        if p:
            print(m, 'annotation:', len(annotation), annotation.astype(int))
            print(m, 'predictions:', len(predictions), predictions.astype(int))

        return ([prec, cov, F1, acc, mcc])

    def map_from_alignment_to_sequence(self, aligned_sequence, scores):
        # print(self.name)
        # print(aligned_sequence)
        # print(scores)
        # print()
        if len(aligned_sequence) != len(scores):
            raise IndexError("length of aligned sequence and scores differ!")
        else:
            # replace = {ord('-'):''}
            # length = len(aligned_sequence.translate(replace))
            length = len([x for x in aligned_sequence if x != '-'])
            output = np.zeros(length, dtype=bool)

            j = 0
            for i, aa in enumerate(aligned_sequence):
                if aa != '-':
                    # try:
                    output[j] = scores[i + 1]
                    j += 1
                # except KeyError as e:
                # 	print(e)
                # 	print(aligned_sequence)
                # 	print(scores)

            return (output)

    def evaluate_transferred_annotations(self):
        """transfer annotations within a FunFam to test how well the annotation of one sequence works for the others in that FunFam"""

        columns = ['prec', 'cov', 'F1', 'acc', 'mcc']
        data = []

        #take annotation of member1 and evaluate on annotation of member2
        for member1 in self.members:
            if not member1.binding_annotation:
                continue

            for member2 in self.members:
                if not member2.binding_annotation:
                    continue
                if member1.id == member2.id:
                    continue

                #take annotation of member1 mapped to sequence of member2 as "predictions" and annotation of member2 as actual annotation
                eval = self.compute_eval(
                    self.map_from_alignment_to_sequence(member2.aligned_sequence,
                                                        self.binding_sites[member1.id]),
                    self.map_from_alignment_to_sequence(member2.aligned_sequence, self.binding_sites[member2.id]))

                data.append(eval)

        df = pd.DataFrame(columns = columns, index = range(0,len(data)), data=data)
        mean_performances = df.mean()

        return mean_performances

    def evaluation(self):
        #########
        # evaluation on basis of aligned sequence, instead map to sequence without gaps/discard predictions for gaps
        # binding_sites are at alignment level!
        #########
        columns = ['prec_cum', 'cov_cum', 'F1_cum', 'acc_cum', 'mcc_cum', 'prec_clust', 'cov_clust', 'F1_clust', 'acc_clust', 'mcc_clust']
        index = self.binding_sites.columns[1:]

        confusion_matrix_cum = []
        confusion_matrix_clust = []
        confusion_matrix_cum_cons = []
        confusion_matrix_clust_cons = []

        for member in self.members:
            if not member.binding_annotation:
                continue

            # without mapping back to individual protein level:
            # eval_cum = self.compute_eval(self.predictions_cum_scores[member.id],self.binding_sites[member.id])
            # eval_clust = self.compute_eval(self.predictions_cluster_coeff[member.id],self.binding_sites[member.id])

            # eval_cum_cons = self.compute_eval(self.predictions_cum_scores['consensus'],self.binding_sites[member.id])
            # eval_clust_cons = self.compute_eval(self.predictions_cluster_coeff['consensus'],self.binding_sites[member.id])

            if self.name == '1999' and member.id == 'P9WHE9':
                p = True
                print(self.name, member.id)
                s = ''.join([x for x in member.aligned_sequence if x != '-'])
                print(len(s), s)
            else:
                p = False

            # mapping back to individual protein level:
            eval_cum = self.compute_eval(
                self.map_from_alignment_to_sequence(member.aligned_sequence, self.predictions_cum_scores[member.id]),
                self.map_from_alignment_to_sequence(member.aligned_sequence, self.binding_sites[member.id]), p, 'cum base', confusion_matrix=confusion_matrix_cum, member_id = member.id)
            eval_clust = self.compute_eval(
                self.map_from_alignment_to_sequence(member.aligned_sequence, self.predictions_cluster_coeff[member.id]),
                self.map_from_alignment_to_sequence(member.aligned_sequence, self.binding_sites[member.id]), p, 'clust base', confusion_matrix=confusion_matrix_clust)

            eval_cum_cons = self.compute_eval(
                self.map_from_alignment_to_sequence(member.aligned_sequence, self.predictions_cum_scores['consensus']),
                self.map_from_alignment_to_sequence(member.aligned_sequence, self.binding_sites[member.id]), p, 'cum consensus',  confusion_matrix=confusion_matrix_cum_cons)
            eval_clust_cons = self.compute_eval(self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                    self.predictions_cluster_coeff[
                                                                                        'consensus']),
                                                self.map_from_alignment_to_sequence(member.aligned_sequence,
                                                                                    self.binding_sites[member.id]), p, 'clust consensus',  confusion_matrix=confusion_matrix_clust_cons)

            # member.set_evaluation_values(eval_cum + eval_clust)
            # member.set_evaluation_consensus(eval_cum_cons + eval_clust_cons)
            member.evaluation_values = (eval_cum + eval_clust)
            # print(self.name,member.id,member.evaluation_values)
            member.evaluation_consensus = (eval_cum_cons + eval_clust_cons)

        data = [member.evaluation_values for member in self.members if member.binding_annotation]
        data_consensus = [member.evaluation_consensus for member in self.members if member.binding_annotation]

        self.evaluation = pd.DataFrame(index=index, columns=columns, data=data)
        self.evaluation_consensus = pd.DataFrame(index=index, columns=columns, data=data_consensus)

        return confusion_matrix_cum + confusion_matrix_clust + confusion_matrix_cum_cons + confusion_matrix_clust_cons

# def evaluation(self):
# 	columns = ['prec_cum','cov_cum','F1_cum','prec_clust','cov_clust','F1_clust']
# 	index = self.binding_sites.columns[1:]

# 	for member in self.members:
# 		trues_cum = sum(self.predictions_cum_scores[member.id])
# 		tp_cum = sum(self.predictions_cum_scores[member.id] & self.binding_sites[member.id])
# 		fp_cum = trues_cum - tp_cum
# 		fn_cum = sum((self.predictions_cum_scores[member.id]==False) & self.binding_sites[member.id])

# 		prec_cum = tp_cum / trues_cum if trues_cum != 0 else 0
# 		cov_cum = tp_cum / (tp_cum + fn_cum) if (tp_cum + fn_cum) != 0 else 0
# 		F1_cum = 2 * (cov_cum * prec_cum) / (cov_cum + prec_cum) if (cov_cum + prec_cum) != 0 else 0

# 		trues_clust = sum(self.predictions_cluster_coeff[member.id])
# 		tp_clust = sum(self.predictions_cluster_coeff[member.id] & self.binding_sites[member.id])
# 		fp_clust = trues_clust - tp_clust
# 		fn_clust = sum((self.predictions_cluster_coeff[member.id]==False) & self.binding_sites[member.id])

# 		prec_clust = tp_clust / trues_clust if trues_clust != 0 else 0
# 		cov_clust = tp_clust / (tp_clust + fn_clust) if (tp_clust + fn_clust) != 0 else 0
# 		F1_clust = 2 * (cov_clust * prec_clust) / (cov_clust + prec_clust) if (cov_clust + prec_clust) != 0 else 0

# 		trues_cum_cons = sum(self.predictions_cum_scores['consensus'])
# 		tp_cum_cons = sum(self.predictions_cum_scores['consensus'] & self.binding_sites[member.id])
# 		fp_cum_cons = trues_cum_cons - tp_cum_cons
# 		fn_cum_cons = sum((self.predictions_cum_scores['consensus']==False) & self.binding_sites[member.id])

# 		prec_cum_cons = tp_cum_cons / trues_cum_cons if trues_cum_cons != 0 else 0
# 		cov_cum_cons = tp_cum_cons / (tp_cum_cons + fn_cum_cons) if (tp_cum_cons + fn_cum_cons) != 0 else 0
# 		F1_cum_cons = 2 * (cov_cum_cons * prec_cum_cons) / (cov_cum_cons + prec_cum_cons) if (cov_cum_cons + prec_cum_cons) != 0 else 0

# 		trues_clust_cons = sum(self.predictions_cluster_coeff['consensus'])
# 		tp_clust_cons = sum(self.predictions_cluster_coeff['consensus'] & self.binding_sites[member.id])
# 		fp_clust_cons = trues_clust_cons - tp_clust_cons
# 		fn_clust_cons = sum((self.predictions_cluster_coeff['consensus']==False) & self.binding_sites[member.id])

# 		prec_clust_cons = tp_clust_cons / trues_clust_cons if trues_clust_cons != 0 else 0
# 		cov_clust_cons = tp_clust_cons / (tp_clust_cons + fn_clust_cons) if (tp_clust_cons + fn_clust_cons) != 0 else 0
# 		F1_clust_cons = 2 * (cov_clust_cons * prec_clust_cons) / (cov_clust_cons + prec_clust_cons) if (cov_clust_cons + prec_clust_cons) != 0 else 0

# 		member.evaluation_values([prec_cum,cov_cum,F1_cum,prec_clust,cov_clust,F1_clust])
# 		member.evaluation_consensus([prec_cum_cons,cov_cum_cons,F1_cum_cons,prec_clust_cons,cov_clust_cons,F1_clust_cons])

# 	data = [member.evaluation_values for member in self.members]
# 	data_consensus = [member.evaluation_consensus for member in self.members]

# 	self.evaluation = pd.DataFrame(index=index, columns=columns, data=data)
# 	self.evaluation_consensus = pd.DataFrame(index=index, columns=columns, data=data_consensus)

######data already aligned
# def align_member_sequences(self,clustalw_exe,temp_dir):
#    	timestamp = '_'.join(str(time.time()).split('.'))

#    	fasta_file = temp_dir+"/"+timestamp+".fasta"
#    	align_file = temp_dir+"/"+timestamp+".align"

#    	records = [SeqRecord(member.sequence,id=member.id) for member in self.members]
#    	SeqIO.write(records, fasta_file, "fasta")

#    	clustalw_cline = ClustalwCommandline(clustalw_exe, infile=fasta_file, outfile=align_file)
#    	std_out, std_err = clustalw_cline()

#    	align = AlignIO.read(align_file, "clustal")

#    	self.alignment = align
#    	index = [member.id for member in members]
#    	self.binding_sites = pd.DataFrame(index=index, columns=[1:len(align[0])])

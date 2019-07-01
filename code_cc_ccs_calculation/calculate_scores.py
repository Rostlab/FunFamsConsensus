import os
import sys
from code_cc_ccs_calculation.protein import Protein
from code_cc_ccs_calculation.file_manager import FileManager


class ScoreCalculator(object):
    def __init__(self, id_list, ec_param):
        self.query_proteins = dict()
        for name in id_list:
            self.query_proteins[name] = Protein(ec_param)
        self.fm = FileManager()

    def calc_scores(self, args):
        # 1) Read files
        print('1) Read files...')
        self.get_evc_files(args.evc_folder, args.ec)
        self.get_fasta_files(args.fasta_folder)
        # 2) Calculate scores
        print('2) Calculate scores for query proteins...')
        self.calculate_scores()
        # 3) Write results
        print('3) Write results...')
        self.write_scores_out(args.out_folder)

    def get_evc_files(self, folder, ec_param):
        """ For each id, determine the different evc files.
        :param folder: Folder where the folders with EVcouplings predictions are stored
        :param ec_param: Parameter specifying how EC scores were calculated (EVcouplings or Freecontact)
        :return:
        """

        for p in self.query_proteins:
            evc_folder = os.path.join(folder, p)
            if ec_param == 'evc':
                evc_file = os.path.join(evc_folder, (p + self.fm.ec_suffix))
            else:
                evc_file = os.path.join(evc_folder, (p + self.fm.di_suffix))
            align_file = os.path.join(evc_folder, (p + self.fm.align_suffix))

            self.query_proteins[p].evc_file = evc_file

            counter = 1
            with open(align_file) as af:
                for line in af:
                    if counter == 2:
                        splitted = line.split(",")
                        length = int(splitted[3])
                        length_cov = int(splitted[4])
                        self.query_proteins[p].length = length
                        self.query_proteins[p].length_cov = length_cov
                    counter += 1

    def get_fasta_files(self, folder):
        """ For each id, read fasta sequence to get length"""

        for p in self.query_proteins:
            fasta_file = os.path.join(folder, (p + '.fasta'))
            seq = '' 
            with open(fasta_file) as ff:
                for line in ff:
                    if '>' not in line:
                        seq += line.strip()

            seq_length = len(seq)
            # print('{p}\t{l}'.format(p=p, l=seq_length))
            self.query_proteins[p].length = seq_length

    def calculate_scores(self):
        """ Calculate cum scores, clustering coefficients, SNAP2/BLOSUM, EVmut/BLOSUM, Solv and Cons.
        :return
        """

        for p in self.query_proteins:
            self.query_proteins[p].calc_scores()

    def write_scores_out(self, folder):
        for p in self.query_proteins:
            ccs_file = os.path.join(folder, 'ccs', (p + self.fm.ccs_suffix))
            cs_file = os.path.join(folder, 'cc', (p + self.fm.cs_suffix))

            with open(ccs_file, 'w') as ccs_out:
                with open(cs_file, 'w') as cs_out:
                    scores = self.query_proteins[p].scores
                    for r in scores.keys():
                        cum_score = scores[r]['Cum']
                        cluster_score = scores[r]['Cluster']
                        ccs_out.write(str(r) + '\t' + str(cum_score) + '\n')
                        cs_out.write(str(r) + '\t' + str(cluster_score) + '\n')


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)

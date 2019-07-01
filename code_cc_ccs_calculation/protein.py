from collections import defaultdict
import os.path
import sys
import numpy


class Protein(object):
    """ A protein in the context of binding site prediction """

    def __init__(self, ec_param):
        self.evc_file = None
        self.length = 0
        self.length_cov = 0
        self.scores = defaultdict(dict)
        self.binding_res = []
        self.seq_dist = 5
        self.threshold = 0.0
        self.average = 0.0
        self.factor = 2

        self.ec_param = ec_param

    def calc_scores(self):
        """ Calculate cum scores, clustering coefficients, SNAP2/BLOSUM, EVmut/BLOSUM, Solv and Cons for a protein.
        If any of the needed files do not exist, the corresponding are not calculated for this protein
        :param:
        :return
        """

        if os.path.isfile(self.evc_file):
            # calculate cumulative coupling scores
            ec_scores = self.get_ec_scores()
            self.calc_thresh_avg(ec_scores)
            # filter by seq distance
            filtered_scores = self.filter_scores(ec_scores)
            cumulative_scores = self.calc_cum_scores(filtered_scores)

            # calculate clustering coefficients
            clustering_coefficients = self.calc_clustering_coefficients(filtered_scores)

            for r in range(1, self.length + 1):
                cum_score = cumulative_scores[r]
                cluster_score = clustering_coefficients[r]
                r = str(r)
                self.scores[r]['Cum'] = cum_score
                self.scores[r]['Cluster'] = cluster_score
        else:
            print(self.evc_file)

    def get_ec_scores(self):
        ec_dict = dict()
        with open(self.evc_file) as f:
            for line in f:
                splitted = line.split(" ")
                comb_pos = splitted[0].strip() + "/" + splitted[2].strip()
                if self.ec_param == 'evc':
                    score = float(splitted[5])
                else:
                    score = float(splitted[4])
                ec_dict[comb_pos] = score
        return ec_dict

    def calc_thresh_avg(self, ec_scores):
        num = self.length_cov * self.factor
        high_scores = list()
        counter = 0
        for key in ec_scores.keys():
            splitted = key.split("/")
            pos1 = int(splitted[0])
            pos2 = int(splitted[1])
            score = ec_scores[key]
            if abs(pos1 - pos2) > self.seq_dist:
                if len(high_scores) >= num:
                    el = high_scores[0]
                    if el < score:
                        del high_scores[0]
                        high_scores.append(score)
                else:
                    high_scores.append(score)
                high_scores.sort()
            else:
                counter += 1
        thresh = high_scores[0]
        avg = float(numpy.mean(high_scores))
        average = round(avg, 2)
        self.threshold = thresh
        self.average = average

    def filter_scores(self, ec_scores):
        filtered_scores = dict()

        for key in ec_scores.keys():
            score = ec_scores[key]
            key_parts = key.split("/")
            pos1 = int(key_parts[0])
            pos2 = int(key_parts[1])
            curr_seq_dist = abs(pos1 - pos2)
            if curr_seq_dist > self.seq_dist and score >= self.threshold:
                filtered_scores[key] = score
        return filtered_scores

    def calc_cum_scores(self, scores):
        """ Calc cumulative coupling scores from a given set of EC scores
        :param scores: list of ec scores
        :return cumulative coupling scores
        """
        cum_scores = dict()
        for i in range(1, self.length + 1):
            sum_ec = 0.0
            num = 0
            for j in range(1, self.length + 1):
                if i <= j:
                    key = str(i) + "/" + str(j)
                else:
                    key = str(j) + "/" + str(i)
                if key in scores.keys():
                    score = scores[key]
                    sum_ec += score
                    num += 1
            ec_strength = 0.0
            if num > 0:
                ec_strength = sum_ec / self.average
                ec_strength = round(ec_strength, 3)
            cum_scores[i] = ec_strength

        return cum_scores

    def calc_clustering_coefficients(self, scores):
        """ Calc clustering coefficients from a given set of EC scores
        :param scores: list of ec scores
        :return clustering coefficients
        """
        coefficients = dict()
        for i in range(1, self.length + 1):
            neighbourhood = set()
            for j in range(1, self.length + 1):
                if i <= j:
                    key = str(i) + "/" + str(j)
                else:
                    key = str(j) + "/" + str(i)
                # if key == '228/234':
                #    print(key + "\t" + str(key in scores.keys())) 
                if key in scores.keys():
                    neighbourhood.add(j)

            num_edges = 0
            for n1 in neighbourhood:
                for n2 in neighbourhood:
                    if n1 <= n2:
                        key = str(n1) + "/" + str(n2)
                    else:
                        key = str(n2) + "/" + str(n1)
                    if key in scores.keys():
                        num_edges += 1
            coeff = 0

            set_size = len(neighbourhood)
            if set_size > 1:
                coeff = num_edges / (set_size * (set_size - 1))
                coeff = round(coeff, 3)
            coefficients[i] = coeff

        return coefficients


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)

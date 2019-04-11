import numpy as np
from Bio.Seq import Seq


class Protein:
    '''object that represents a protein sequence and its annotated and predicted binding sites'''

    def __init__(self, id, aligned_sequence, funfam, start_region, end_region, start_segment, end_segment,
                 evc_positions, di_positions, binding_annotation, start_seq, end_seq, cum_coup_cutoff,
                 clust_coeff_cutoff):
        self.id = id
        self.aligned_sequence = aligned_sequence
        self.funfam = funfam
        self.start_region = start_region
        self.end_region = end_region
        self.start_segment = start_segment
        self.end_segment = end_segment
        self.evc_positions = evc_positions
        self.di_positions = di_positions
        self.cum_coup_cutoff = cum_coup_cutoff
        self.clust_coeff_cutoff = clust_coeff_cutoff
        self.no_predictions = False
        self.binding_annotation = binding_annotation
        if not binding_annotation:
            self.start_seq = start_seq
            self.end_seq = end_seq

    def binding_sites(self, sites):
        self.binding_sites = np.array(sites)

    def sequence(self):
        self.sequence = ''.join([str(x) for x in self.aligned_sequence if x != '-'])

    def di_scores(self, scores):
        '''to do if necessary'''
        pass

    def fill_missing_scores(self, scores):
        if self.start_region == self.start_segment and self.end_region == self.end_segment:
            return (scores)
        elif self.start_region == self.start_segment and self.end_region > self.end_segment:
            new_scores = np.zeros(self.end_region - self.end_segment)
            scores = np.append(scores, new_scores)
            return (scores)
        elif self.start_region < self.start_segment and self.end_region == self.end_segment:
            new_scores = np.zeros(self.start_segment - self.start_region)
            scores = np.append(new_scores, scores)
            return (scores)
        elif self.start_region < self.start_segment and self.end_region > self.end_segment:
            new_scores_start = np.zeros(self.start_segment - self.start_region)
            new_scores_end = np.zeros(self.end_region - self.end_segment)
            scores = np.append(new_scores_start, scores)
            scores = np.append(scores, new_scores_end)
            return (scores)
        elif self.start_region == self.start_segment and self.end_region < self.end_segment:
            # more scores than positions in sequence, but additional scores are at the end of the sequence
            scores = scores[:self.end_region - self.start_region + 1]
            # print(self.id,"cutting scores",len(scores))
            return (scores)
        else:
            print("scores not filled!", self.id)
            print("weird shift of start/end segments:")
            print("\tstart region:", self.start_region)
            print("\tend region:", self.end_region)
            print("\tstart segment:", self.start_segment)
            print("\tend segment:", self.end_segment)

            return (scores)

    def cum_scores(self, file):
        with open(file) as f:
            content = f.readlines()
        data = [x.strip() for x in content]
        cs = np.zeros(len(data))

        for elem in data:
            index, score = elem.split('\t')
            cs[int(index) - 1] = float(score)

        scores = self.process_scores(cs, "cumulative coupling")

        if len(scores) == 0:
            self.no_predictions = True
            return ()

        self.cum_scores = self.map_scores_to_alignment(scores, 'cum')
        self.cum_scores_predictions = self.cum_scores > self.cum_coup_cutoff  # cum_scores cutoff for binding residues as per Maria's thesis

    def process_scores(self, scores, kind):
        if self.binding_annotation:
            scores = self.fill_missing_scores(scores)
        else:
            scores = self.fill_missing_scores(scores)
            #			print("\t",self.id,"cutting",kind,"scores to conform to FunFam sequence segment")
            try:
                scores = scores[self.start_seq - 1:self.end_seq]
            except IndexError as e:
                print(e)
                print(len(scores), self.start_seq, self.end_seq)
        return (scores)

    def cluster_coeff(self, file):
        with open(file) as f:
            content = f.readlines()
        data = [x.strip() for x in content]
        cc = np.zeros(len(data))
        for elem in data:
            index, coeff = elem.split('\t')
            cc[int(index) - 1] = float(coeff)

        scores = self.process_scores(cc, "clustering coefficient")

        if len(scores) == 0:
            self.no_predictions = True
            return ()

        self.cluster_coeffs = self.map_scores_to_alignment(scores, 'clu')
        self.cluster_coeff_predictions = self.cluster_coeffs > self.clust_coeff_cutoff  # clustering_coeff cutoff for binding residues as per Maria's thesis

    def alignment_stats(self, file):
        with open(file) as f:
            content = f.readlines()
        data = [x.strip() for x in content]

        keys = data[0].split(',')
        values = data[1].split(',')

        self.alignment_stats = dict(zip(keys, values))

    def map_scores_to_alignment(self, scores, caller):
        # if caller == 'cum' and len(scores) != self.end_region - self.start_region + 1:
        # 	print("length of scores and length of sequence differ:",self.id,len(scores),str(self.end_region - self.start_region + 1))
        counter = 0
        mapped_scores = np.zeros(len(self.aligned_sequence))
        printed = False
        try:
            for i, aa in enumerate(self.aligned_sequence):
                if aa != '-':
                    mapped_scores[i] = scores[counter]
                    counter += 1
        except IndexError as e:
            if caller == 'cum':
                printed = True
                print(e)
                print("\t", self.id)
                print("\t", self.funfam)
                print("\t", self.start_region, self.end_region)
                print("\t", 'aligned sequence length:', len(self.aligned_sequence))
                print("\t", 'sequence length:', len([x for x in self.aligned_sequence if x != '-']))
                print("\t", 'scores length:', len(scores))
                print("\t", 'i:', i, '\tcounter:', counter)
                print("\t", "same number of positions:", len(self.di_positions) == len(self.evc_positions))
                print("\t", "di positions:", len(self.di_positions), list(self.di_positions), '\n')
                print("\t", "evc positions:", len(self.evc_positions), self.evc_positions, '\n')

        # if caller == 'cum' and not printed and len([x for x in self.aligned_sequence if x != '-']) != len(scores):
        # 	print("scores != sequence length")
        # 	print("\t",self.id)
        # 	print("\t",self.funfam)
        # 	print("\t",self.start_region,self.end_region)
        # 	print("\t",'aligned sequence length:',len(self.aligned_sequence))
        # 	print("\t",'sequence length:',len([x for x in self.aligned_sequence if x != '-']))
        # 	print("\t",'scores length:',len(scores))
        # 	print("\t",'i:',i,'\tcounter:',counter)
        # 	print("\t","same number of positions:",len(self.di_positions) == len(self.evc_positions))
        # 	print("\t","di positions:",len(self.di_positions),list(self.di_positions),'\n')
        # 	print("\t","evc positions:",len(self.evc_positions),self.evc_positions,'\n')
        return (mapped_scores)

    def set_evaluation_values(self, values):
        self.evaluation_values = values

    def set_evaluation_consensus(self, values):
        self.evaluation_consensus = values

# def map_sites_to_alignment(self,sites):
# 	out = np.zeros(len(self.aligned_sequence))
# 	for bs in sites:
# 		gaps = 0
# 		counter = 0
# 		for i,aa in enumerate(self.aligned_sequence):
# 			if aa == '-':
# 				gaps += 1
# 			else:
# 				counter += 1
# 			if counter == bs:
# 					break
# 			out[i] = bs + gaps
# 	return(out)

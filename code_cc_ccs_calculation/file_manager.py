""" Class specifying file endings """


class FileManager(object):

    @property
    def di_suffix(self):
        return ".di_scores"

    @property
    def ec_suffix(self):
        return '_ECs.txt'

    @property
    def align_suffix(self):
        return "_alignment_statistics.csv"

    @property
    def ccs_suffix(self):
        return ".cum_scores"

    @property
    def cs_suffix(self):
        return ".cluster_coeff"

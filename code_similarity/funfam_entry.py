"""
class that represents FunFam entries
"""


class FunFamEntry(object):

    def __init__(self, id, funfam, superfamily, start, end, uniprot_ids, ec_ids, sequence, sequence_aligned_by_funfam):
        self.id = id
        self.funfam = funfam
        self.superfamily = superfamily
        self.start = start
        self.end = end
        self.uniprot_ids = uniprot_ids
        self.ec_ids = ec_ids
        self.sequence = sequence
        self.aligned_sequence_funfam = sequence_aligned_by_funfam
        self.sites = []

    def __str__(self):
        return (self.id + ";\t" + str(self.uniprot_ids) + ";\t" + str(
            self.ec_ids) + ";\t" + self.funfam + ";\t" + self.superfamily + ";\t" + str(self.start) + ";\t" + str(
            self.end) + ";\n" + str(self.aligned_sequence_funfam) + ";\n" + str(self.sites) + ";\n" + str(
            self.mapped_sites_funfam))

    def __eq__(self, other):
        """
        checks if two funfam_entry objects represent the same domain sequence
        :param other: funfam_entry
        :return: boolean
        """
        return self.superfamily == other.superfamily and self.funfam == other.funfam and self.id == other.id

    def __hash__(self):
        return hash((self.id, self.superfamily, self.funfam, self.start, self.end))

    def add_binding_sites(self, uniprot_binding_site_mapping):
        for uniprot_id in self.uniprot_ids:
            self.sites.append(uniprot_binding_site_mapping.get(uniprot_id))

    def process_sites(self):
        i = 0
        index = -1
        relevant_sites = None
        for binding_sites in self.sites:
            index += 1
            if binding_sites is None:
                continue
            if not binding_sites:
                continue
            i += 1
            relevant_sites = binding_sites

        if i > 1:
            print("ambiguity processing sites:", self.id, self.funfam, self.superfamily)
            # print(self.sites)
            # raise ValueError("multiple binding site annotations for entry")
            return (False)
        if relevant_sites is None:
            self.sites = None
        else:
            self.sites = [site for site in relevant_sites if (site >= self.start and site <= self.end)]
            self.binding_site_id = self.uniprot_ids[index]

        return (True)

    def write_info(self, group):
        print(self.id, self.funfam, self.superfamily, self.start, self.end, self.binding_site_id)
        print("uniprot ids:", self.uniprot_ids)
        print("binding sites", self.sites)
        if group == "ec":
            print("mapped sites ec:", self.mapped_sites_ec)
        elif group == "funfam":
            print("mapped sites funfam:", self.mapped_sites_funfam)
        print("ec numbers:", self.ec_ids)

    def map_binding_sites(self, group):
        '''maps binding sites from the sequence to the ec alignment level'''
        if group == "funfam":
            sequence = self.aligned_sequence_funfam
        elif group == "ec":
            sequence = self.aligned_sequence_ec
        else:
            raise ValueError("group has to be funfam or ec")

        out = []
        for bs in self.sites:
            gaps = 0
            counter = self.start - 1
            for aa in sequence:
                if aa == '-':
                    gaps = gaps + 1
                else:
                    counter = counter + 1
                if counter == bs:
                    break
            out.append(bs + gaps - self.start + 1)

        if group == "funfam":
            self.mapped_sites_funfam = out
        elif group == "ec":
            self.mapped_sites_ec = out

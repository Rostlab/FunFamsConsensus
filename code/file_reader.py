"""this script provides file reading for multiple input formats
"""

from .funfam_entry import FunFamEntry
from Bio import SeqIO

GAP_CHARACTER = '-'

class FunFamReader(object):

    """reads a FunFam file in .aln format and creates a list of FunFamEntry objects from it
    note: individual entries are omitted if they have no UNIPROT ID in the header
    """

    def __init__(self, filepath, superfamily, funfam):
        self.path = filepath
        self.superfamily = superfamily
        self.funfam = funfam
        self.entries = list()
        self.num_rejected_entries = 0

    def read(self):
        for seq_record in SeqIO.parse(self.path, 'fasta'):
            try:
                uniprot_ids = self.get_uniprot_ids(seq_record.description)
                ec_ids = self.get_ec_ids(seq_record.description)
                unique_id = seq_record.id.split('|')[2].split('/')[0]
                start,end = self.get_sequence_range(seq_record.id)
                sequence_aligned_by_funfam = seq_record.seq
                sequence = seq_record.seq.ungap(GAP_CHARACTER)
                if len(sequence) != end - start + 1:
                    # omit records where sequence length does not match annotated length
                    raise ValueError("annotated sequence length does not match real sequence length")

            except ValueError as e:
                self.num_rejected_entries += 1
                #print('[FunFamReader]:\t rejected an entry when reading!', e)
                continue

            #print('[FunFamReader]:\t','accepted an entry!')
            entry = FunFamEntry(unique_id, self.funfam, self.superfamily, start, end, uniprot_ids, ec_ids, sequence, sequence_aligned_by_funfam)
            self.entries.append(entry)

    def get_uniprot_ids(self, seq_record_description):
        """
        parses list of UNIPROT IDs from SeqRecord.description
        """
        descriptors = seq_record_description.split('UNIPROT=')
        if len(descriptors) == 1:
            raise ValueError("not UNIPROT IDs annotated")
        uniprot_ids = descriptors[1].split(';')[0].split(',')
        return uniprot_ids

    def get_ec_ids(self, seq_record_description):
        """
        parses the EC numbers from SeqRecord.description
        None if no EC numbers are annotated
        """
        descriptors_id_list = seq_record_description.split('EC=')
        if len(descriptors_id_list) > 1:
            ec_ids_all = descriptors_id_list[1].split(';')[0].split(',')
            ec_ids = [ec_id for ec_id in ec_ids_all if not '-' in ec_id]
            if not ec_ids: ec_ids = None
        else:
            ec_ids = None
        return ec_ids

    def get_sequence_range(self, seq_record_id):
        """
        parses start and end of sequence from SeqRecord.id
        """
        seq_range = seq_record_id.split('/')[1]
        if '_' in seq_range or seq_range.startswith('-'):
            # raise error for unexpected characters
            raise ValueError("failed to parse sequence range")
        start, end = map(int, seq_range.split('-'))
        return start,end

def read_uniprot_binding_site_mapping(file):
    uniprot_binding_site_mapping = dict()
    with open(file, "r") as f:
        for line in f:
            line = line.strip()
            line_split = line.split("\t")
            uniprot_binding_site_mapping[line_split[0]] = list(map(int, line_split[1:]))

    return uniprot_binding_site_mapping
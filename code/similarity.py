"""
launch the similarity computation of binding residues
"""

import os
import argparse
import numpy as np
import pickle
from itertools import chain
from collections import defaultdict
from funfam_project.code.utils import process_funfam_entries, similarity, get_group_mapping
from funfam_project.code.file_reader import FunFamReader, read_uniprot_binding_site_mapping


def main():
    groups = ['funfam', 'ec']
    usage_string = 'python similarity.py'
    parser = argparse.ArgumentParser(description=__doc__, usage=usage_string)
    parser.add_argument("-families", dest="family_dir", help="directory with FunFam families", required=True)
    parser.add_argument("-sites", dest="binding_site_file", help="file with UNIPROT ID to binding site mapping",
                        required=True)
    parser.add_argument("-groupby", dest="grouping_keyword",
                        help="keyword by which sequences are grouped for similarity calculation", choices=groups,
                        required=True)
    parser.add_argument("-limit", dest="limit_keyword",
                        help="keyword by which sequences are separated for similarity calculation", choices=groups)
    parser.add_argument("-align", dest="alignment_path",
                        help="directory to which sequence alignments will be written")
    parser.add_argument("-clustalw", dest="clustalw_command",
                        help="command to call clustalw on this system")

    args = parser.parse_args()
    print("[ARGUMENTS]")
    print("\tdirectory with FunFam data:", args.family_dir)
    print("\tfile with binding site mapping:", args.binding_site_file)
    print("\tgroupby:", args.grouping_keyword)
    print("\tlimitby:", args.limit_keyword)
    print("\talign:", args.alignment_path)
    print("\tclustalw:", args.clustalw_command)
    print("[ARGUMENTS END]")

    uniprot_binding_site_mapping = read_uniprot_binding_site_mapping(args.binding_site_file)
    funfam_entries = list()
    rejected_entries = 0

    print("reading FunFam data...")

    # i = 0
    # j = 0
    # for superfamily in os.listdir(args.family_dir):
    #     print('\tsuperfamily:', i, superfamily)
    #     if ".DS_Store" in superfamily:
    #         continue
    #     for funfam_file in os.listdir(os.path.join(args.family_dir, superfamily)):
    #         # print('\tj', j, funfam_file)
    #         if ".DS_Store" in funfam_file:
    #             continue
    #         funfam = funfam_file.split(".")[0]
    #         reader = FunFamReader(os.path.join(args.family_dir, superfamily, funfam_file), superfamily, funfam)
    #         reader.read()
    #         rejected_entries += reader.num_rejected_entries
    #         # print('\tnumber of read entries:', len(reader.entries))
    #         new_entries = process_funfam_entries(reader.entries, uniprot_binding_site_mapping)
    #         # print('\tnumber of accepted entries:', len(new_entries))
    #         funfam_entries += new_entries
    #         j += 1
    #         # if len(funfam_entries) >= 100:
    #         #     break
    #     i += 1
    #     # if len(funfam_entries) >= 100:
    #     #     break
    #
    # print('superfamilies:', i, 'funfams:', j, 'found entries:', len(funfam_entries), 'unique:',
    #       len(set(funfam_entries)))
    # print('rejected entries:', rejected_entries)

    pickle_file = 'C:\\Users\\Linus\\LRZ Sync+Share\\UniversitätMünchen\\Bioinformatik\\6. Semester\\Bachelorarbeit\\funfam_project\\data\\ff2.p'
    # print('start serializing funfams')
    # pickle.dump(funfam_entries, open(pickle_file, 'wb'))
    # print('done serializing funfams')

    print("reading pickle...")
    funfam_entries = pickle.load(open(pickle_file, 'rb'))
    print("done reading.")

    group_mapping = get_group_mapping(funfam_entries, args.grouping_keyword, args.limit_keyword)

    similarities = []
    num_used_entries = 0

    print('number of groups:', len(group_mapping), '\nnumber of entries:', len(set(chain(*group_mapping.values()))))

    for group, entries in group_mapping.items():
        if len(entries) < 2:
            continue

        sim = similarity(group, entries, args.grouping_keyword, args.limit_keyword, args.alignment_path,
                         args.clustalw_command)
        similarities.append((group, sim))

    scores = np.array([x[1][1] for x in similarities if x is not None])
    num_members = np.array([len(x[1][0]) for x in similarities if x is not None])
    num_used_entries = sum((len(x[1][0]) for x in similarities))

    print('similarity:', scores[num_members != 1].mean(), '\nnumber of groups:', len(scores[num_members != 1]))
    print('used entries:', num_used_entries)


if __name__ == '__main__': main()

"""
launch the similarity computation of binding residues
"""

#to use precalculated representation of the funfam data:
#-pickle
#"C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\funfam_project\data\ff2.p"
#-limit
#"funfam"
#-pickle
#"C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\FunFams\data\ff3.p"

import os
import argparse
import numpy as np
import pickle
from itertools import chain
from collections import defaultdict
from code_similarity.utils import process_funfam_entries, similarity, get_group_mapping, standard_error
from code_similarity.file_reader import FunFamReader, read_uniprot_binding_site_mapping


def main():
    groups = ['funfam', 'ec', 'pfam', 'prosite', 'funfam-on-ec-subset', 'funfam-on-pfam-subset', 'funfam-on-prosite-subset']
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
    parser.add_argument("-pickle", dest="pickle_file", help="for debugging: read FunFam data from pickle file")
    parser.add_argument("-pfam_map", dest="uniprot_pfam_file", help="if groupby == pfam")
    parser.add_argument("-prosite_map", dest="uniprot_prosite_file", help="if groupby == prosite")
    parser.add_argument("-entries", dest="file_entries_to_use", help="if only some entries should be used")
    parser.add_argument("-write_used_entries", dest="write_used_entries",help="print the used entries")

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
    uniprot_ids = set()

    if args.pickle_file is None:
        print("reading FunFam data...")
        i = 0
        j = 0
        for superfamily in os.listdir(args.family_dir):
            print('\tsuperfamily:', i, superfamily)
            if ".DS_Store" in superfamily:
                continue
            for funfam_file in os.listdir(os.path.join(args.family_dir, superfamily)):
                # print('\tj', j, funfam_file)
                if ".DS_Store" in funfam_file:
                    continue
                funfam = funfam_file.split(".")[0]
                reader = FunFamReader(os.path.join(args.family_dir, superfamily, funfam_file), superfamily, funfam)
                reader.read()
                rejected_entries += reader.num_rejected_entries
                # print('\tnumber of read entries:', len(reader.entries))
                new_entries = process_funfam_entries(reader.entries, uniprot_binding_site_mapping)
                # print('\tnumber of accepted entries:', len(new_entries))
                funfam_entries += new_entries
                for entry in new_entries:
                    for u_id_entry in entry.uniprot_ids:
                        uniprot_ids.add(u_id_entry)
                j += 1
                # if len(funfam_entries) >= 100:
                #     break
            i += 1
            # if len(funfam_entries) >= 100:
            #     break

        print('superfamilies:', i, 'funfams:', j, 'found entries:', len(funfam_entries), 'unique:',
              len(set(funfam_entries)))
        print('rejected entries:', rejected_entries)

    else:
        print("reading pickle...")
        funfam_entries = pickle.load(open(args.pickle_file, 'rb'))
        print("done reading.")

    # print('start serializing funfams')
    # pickle.dump(funfam_entries, open('data\\ff4.p', 'wb'))
    # print('done serializing funfams')

    # with open('used_uniprot_ids_similarity', 'a') as the_file:
    #     for u_id in uniprot_ids:
    #         the_file.write(u_id+'\n')


    group_mapping = get_group_mapping(funfam_entries, args.grouping_keyword, args.limit_keyword, args.uniprot_pfam_file, args.uniprot_prosite_file, file_entries_to_use=args.file_entries_to_use)

    similarities = []
    num_used_entries = 0

    print('number of groups:', len(group_mapping), '\nnumber of entries:', len(set(chain(*group_mapping.values()))))

    for group, entries in group_mapping.items():
        if len(entries) < 2:
            continue

        sim = similarity(group, entries, args.grouping_keyword, args.limit_keyword, args.alignment_path,
                         args.clustalw_command)
        similarities.append((group, sim))

    print(similarities[0])

    scores = np.array([x[1][1] for x in similarities if x is not None])
    num_members = np.array([len(x[1][0]) for x in similarities if x is not None])
    num_used_entries = sum((len(x[1][0]) for x in similarities))
    num_used_entries_new = sum((len(x[1][0]) for x in similarities if len(x[1][0]) is not 1))

    print('similarity:', scores[num_members != 1].mean(),'+-', standard_error(scores[num_members != 1]), '\nnumber of groups:', len(scores[num_members != 1]))
    print('used entries:', num_used_entries)
    print('used entries new:', num_used_entries_new)

    if args.write_used_entries is not None and args.write_used_entries == 'True':
        with open(os.path.join(args.alignment_path, 'used_entries_'+args.grouping_keyword+'.txt'), 'w') as f:
            for data in [x[1] for x in similarities]:
                if len(data[0]) < 2:
                    #print('small funfam')
                    continue
                for superfamily, funfam, e_id, uniprot_id in data[0]:
                    f.write(superfamily+','+funfam+','+e_id+','+uniprot_id+'\n')

if __name__ == '__main__': main()

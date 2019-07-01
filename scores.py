"""
DESCRIPTION:

Development and training of bindPredict to predict binding site residues
"""

import argparse
import sys

from code_cc_ccs_calculation.calculate_scores import ScoreCalculator


def main():
    usage_string = 'python get-scores.py evcouplings_res/ id_list'

    # parse command line options
    parser = argparse.ArgumentParser(description=__doc__, usage=usage_string)
    parser.add_argument('evc_folder', help='Folder with protein folders containing EVcouplings results')
    parser.add_argument('fasta_folder', help='Folder containing fasta sequences')
    parser.add_argument('id_file', help='File containing a list of ids to calculate scores for')
    parser.add_argument('ec', help='Method used to calculate ECs [evc|freecontact]')
    parser.add_argument('out_folder', help='Folder to write output to')
    args = parser.parse_args()

    id_list = []
    with open(args.id_file) as id_in:
        for line in id_in:
            id_list.append(line.strip())

    sc = ScoreCalculator(id_list, args.ec)
    sc.calc_scores(args)


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)


if __name__ == "__main__":
    main()

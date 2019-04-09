This repository provides functionality to compute the binding residue similarity of sequences in the FunFam dataset.

# Usage
## Computing binding residue similarity
The python program `similarity.py` has the following command line parameters:
* `-families`
   the path to a folder containing the FunFam dataset (with one sub-directory per superfamily)
* `-sites`
   the path to a file with a mapping of UNIPROT IDs to binding site annotation (see /data)
* `groupby` [funfam|ec]
   the way in which the sequences should be grouped for similarity computations.
* `limit` [funfam|ec]   optional
   the groups which should not occurr multiple times within the group specified by `groupby`
* `align`               optional
   the path to a directory in which data for the generation of multiple sequence alignments can be stored
* `clustalw`            optional
   the command to call the clustalw MSE tool, necessary only if `groupby` == ec
   
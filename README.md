This repository provides functionality to compute the binding residue similarity of sequences in the FunFam dataset.

# Usage
## Compute binding residue similarity
`similarity.py` has the following command line parameters:
* `-families`  
    the path to a directory containing the FunFam dataset (with one sub-directory per superfamily).
* `-sites`  
    the path to a file with a mapping of UNIPROT IDs to binding site annotation (see /data).
* `-groupby` [funfam|ec]  
    the way in which the sequences should be grouped for similarity computations.
* `-limit` [funfam|ec]   optional  
    the groups which should not occur multiple times within the group specified by `groupby`.
* `-align`               optional  
    the path to a directory in which data for the generation of multiple sequence alignments can be stored.
* `-clustalw`            optional  
    the command to call the external clustalw MSE tool, necessary only if `groupby` == ec.

## Build and evaluate consensus prediction
`prediction.py` has the following command line parameters:
* `-consensus`  
    the consensus cut-off at positions are classified as binding.
* `-cc`  
    the cut-off above which a position is classified as binding by its clustering coefficient.
* `-ccs`  
    the cut-off above which a position is classified as binding by its cumulative coupling score.
* `-uniprot_ids`  
    path to a file containing all UNIPROT IDs for which data is available, one id per line.
* `-mapping`  
    path to a file with a mapping of FunFams to UNIPROT IDs.
* `-evc_info`  
    path to a directory with output files from [EVcouplings](http://evfold.org/evfold-web/newmarkec.do) (\_final.outcfg, .alignment_statistics.csv), [FreeContact](https://rostlab.org/owiki/index.php/FreeContact) (.di) as well as [bindPredict](https://github.com/Rostlab/bindPredict) (.cum_scores, .cluster_coeff). Data for each UNIPROT ID ought to be in a seperate subdirectory.
* `-funfam_data`  
    path to a file in FASTA FunFam format including mapped binding sites for each entry.
* `-families`  
    the path to a directory containing the FunFam dataset (with one sub-directory per superfamily).
* `-out`  
    the path to a directory to which output files will be written.

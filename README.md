# Multi-Omics NETwork inference (MONET)
A Bayesian nonparametric integration approach for functional regulatory network inference in haematopoietic development from multi-omics data

By: Laia Carolina Alcalà

Supervised by: Anas Rana and Jean-Baptiste Cazier

ESCI-UPF and University of Birmingham (Centre for Computational Biology)


# Repository content:
## Readme:
Information about the repository and folders content.

## Report
----

## Data:
It contains all the preprocessed data in which the model will be applied
- ATAC-seq data:
  - ATAC_Sp1hyp_binary_table.txt: Binarized measure of chromatin state over the different stages of hematopoietic development in mice. A 0 value represents a closed status and a 1 value an open status.
  - ATAC_Sp1hyp_time.txt
  - ATAC_WT-time.txt
  - ATAC_WT_binary_table.txt
- ChIP-seq data:
  - ES_ChIP_binary_table.txt
  - ES_chip_binding_sp1_sp3_binary_matrix.txt
  - ES_fpkm_WT_Sp1hyp.txt
- RNA-seq data:
  - ge_comb_Sp1WT_Sp1hyp_fpkm.txt
  - ge_full_Sp1WT.txt
  - ge_full_Sp1hyp.txt

- binary_tf.txt

## Scripts
- Model.stan : Stan implementation of the Bayesian non parametric model

## Figures
All the figures created by me used in the report
- Experiment_scheme.png : A visual representation of the experiment on which the model will be applied
- Model_scheme.png : A visual representation of the model structure

## Results:
---


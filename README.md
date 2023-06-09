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
  - ATAC_Sp1hyp_binary_table.txt: Binary table containing the chromatin state over the 5 different time-points with hypomorph Sp1 TF. Rows are genes and columns are the time-points. A 0 value represents a closed status and a 1 value an open status.
  - ATAC_Sp1hyp_time.txt: Binarized measure of the chromatin state over the 4 time-points (it does not contain t=0) with hypomorph Sp1 TF. Rows are genes and columns are the time-points. A 0 value represents a closed status and a 1 value an open status.
  - ATAC_WT-time.txt: Binarized measure of the chromatin state over the 4 time-points (it does not contain t=0) with Wild Type Sp1 TF. Rows are genes and columns are the time-points. A 0 value represents a closed status and a 1 value an open status.
  - ATAC_WT_binary_table.txt: Binary table containing the chromatin state over the 5 different time-points with  Wild Type Sp1 TF. Rows are genes and columns are the time-points. A 0 value represents a closed status and a 1 value an open status.

- ChIP-seq data:
  - ES_ChIP_binary_table.txt: Binary table containing the interaction of the different TFs at ES time-point (t=0) with the genes. A 0 value means non-interaction and a 1 value means interaction.
  - ES_chip_binding_sp1_sp3_binary_matrix.txt: Binarized measure of the binding of Sp1 and Sp3 TFs for the different genes at ES time-point (t=0). A 0 value means non-binding and a 1 value means binding.
  - binary_tf.txt: Binary table containing the interaction of the different TFs at the 5 time-points with genes. A 0 value means non-interaction and a 1 value means interaction.
 
- RNA-seq data:
  - ge_full_Sp1WT.txt: Gene expression of genes over the 5 different time-points with Wild Type Sp1 TF (measured as FPKM).
  - ge_full_Sp1hyp.txt: Gene expression of genes over the 5 different time-points with hypomorph Sp1 TF (measured as FPKM).
  - ge_comb_Sp1WT_Sp1hyp_fpkm.txt: Table containing the gene expression of genes over the 5 different time-points with both Wild Type and hypomorph Sp1 TF (measured as FPKM).
  - ES_fpkm_WT_Sp1hyp.txt: Table containing the gene expression of genes at ES time-point (t=0) for both Wild Type and hypomorph Sp1 TF (measured as FPKM).


## Scripts
- Model.stan : Stan implementation of the Bayesian non parametric model
- Monet-sim: Monet data simulation function
- monet-sim-run.R : Data simulation with 2000 genes and 30 TF using the Monet package and model fitting to this data.
- comparison.R : Statictical analysis and comparison of the processed and simulated dataset
- monet_sim_yg : Exctraction of the estimated parameters of the simulated data, computation of the estimated gene expression (yg) of the simulated data and data visualization. y0 and ag are taken from the simulated data and bgk, xk, wk and ek from the point estimates of the monet fitted model.
- monet-obs_data_run.R: Model fitting to the preprocessed (observed data), exctraction of the estimated parameters data visualization.

## Figures
It contains all the figures created by me used in the report and plots used in the scientific report.
-

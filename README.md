The repository contains scripts showing the pipeline to analyze metabolomics measurements (from a Nightingale NMR panel) of plasma samples collected during a meal challenge test. The data is assumed to
be arranged as a third-order tensor with modes: subjects, time and metabolites, and analyzed using a CANDECOMP/PARAFAC (CP) tensor factorization model. 
- script_CP_NMR is the main function. It shows the preprocessing steps and fits the model using cp_wopt. It also calls the CP_replicability function to assess the replicability and 
shows how to relate the subject factor matrix to various meta variables of interest such as body mass index (BMI).
- functions/CP_replicability checks the replicability of the factors extracted by an R-component CP model.

The scripts are based on the assumption that the data is stored as a dataset object (https://eigenvector.com/software/dataset-object/). They also make use of functions from the Tensor Toolbox 
(https://www.tensortoolbox.org/) and Poblano toolbox (https://github.com/sandialabs/poblano_toolbox).

While the data we have used is not available publicly, the same pipeline can be used to analyze similar time-resolved data sets. For more details about the data and the analysis, please see the full paper:
-  S. Yan, L. Li, D. Horner, P. Ebrahimi, B. Chawes, L. O. Dragsted, M. A. Rasmussen, A. K. Smilde,  E. Acar. Characterizing human postprandial metabolic response using multiway data analysis, 2023
https://doi.org/10.1101/2023.08.31.555521

  

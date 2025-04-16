# Multivariate Analyses

This directory contains data, code, and supporting information for various multivariate analyses appearing in the manuscript **"Continent-wide Differentiation of Fitness Traits and Patterns of Climate Adaptation among European Populations of _Drosophila melanogaster_."**.

## Data and Code

Below is an overview of what each of the different data files (found within /MultivariateAnalyses/Data) and scripts (found within /MultivariateAnalyses/Code) do:

- The data file *all_models_compound_coefs_forMultivariateAnalysis.csv* contains the trait values estimated from linear models. The file is used as in input for the phenotypic PCA which can be reproduced with the R script *mv_01_Phenotype_PCA.R*.
    - This script generates the output files *F9_drosEU.Rdata*, *FmaxP_drosEU.Rdata*, and *M9_drosEU.Rdata*, as well as the tables *F9_drosEU_PCcoords.csv*, *FmaxP_drosEU_PCcoords.csv*, and *M9_drosEU_PCcoords.csv*. These files are also located in /MultivariateAnalyses/Data, as they are used as input for subsequent scripts.

- The file *all_models_coefsDietControl_forMultivariateAnalysis.csv* contains trait values estimated using a subset of the data (from labs where the protein to carbohydrate ratio was similar). This data file is the input for script *mv_02_Phenotype_PCA_DietControl.R*.

- The three \**_PCcoords.csv files* are used as inputs for scripts *mv_03A_Climate_PCA_30days_ClimatePhenotypeAssociation.R* and *mv_03B_Climate_PCA_30years_ClimatePhenotypeAssociation.R*, which contain the code necessary to carry out climate PCAs, as well as analyses correlating phenotypic principal components and climate principal components. 

- Scripts mv_03A and mv_03B generate the data files *all_traits_30d_WS.Rdata* and *all_traits_30y_WS.Rdata* which are used as inputs for script *mv_03C_ClimatePhenotype_Permutation.R*. This script tests whether correlations between phenotypic PCs and climatye PCs were greater than expected by chance by using a permutation-based approach.

- The three \**_drosEU.Rdata files* and two *all_traits_30\*_WS.Rdata* files are used as input for script *mv_03D_ClimatePhenotypeAssociation_figures.R* which generates figures showing how phenotypic PCs correlate with climate PCs.

## Supporting Information

- Further files relating to the phenotypic PCAs (with all data and diet controlled data), the Discriminant Function Analysis (DFA), the climate PCA, and analyses correlating phenotypic PCs with climate PCs can be found in /MultivariateAnalyses/SI. Many of these files are R data objects, figures, or tables that appear in the main DrosEU_PhenotypingWG.html.

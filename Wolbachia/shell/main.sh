# Main shell script for project DrosEU_Wolbachia

cd /media/inter/mkapun/projects/DrosEU_PhenotypingWG

## copy Phenotype data
cp PhenotypingWG/InfoTables/Dros*.csv \
  Wolbachia/data

## Isolate the phenotyping data from the BIG R object
Rscript Wolbachia/scripts/GetData.R

## Perform Wolbachia-specific analysis for all phenotypes
Rscript Wolbachia/scripts/WolbachiaAnalysis.R

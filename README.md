### Integration and analysis of Lung macrophage cells from patients with COVID-19, Idiopathic pulmonary fibrosis,  COPD. 
The following notebooks described the data integration, data exploration and *proximity analysis* steps of Methods as implemented in

**Wendisch, Dietrich, Mari *et al.* (2021). (*In press*).**

Scripts in `notebooks` indicate the main routines including additional data from Budinger *et al.*. Main processing steps and analyses are divided into three main steps:
  1. Data Preparation of Adams/Morse/Bharat and our work. H5AD are generated for each dataset using the public raw data, and then deployed in `data`. Pre-processed H5AD files can be downloaded from **this directory**.
  2. Data Integration of all datasets using (scVI), and patient/sample IDs as batches.
  3. Downstream analyses with the integrated object analyses. Mainly:
  3.1. Gene expression quantification post-processing.
  3.2. Quantification of associations between macrophages from COVID-19 patients and IPF, using proximities in the kNN graph from the integrated cells as a reference for mapping conditions.

##### Dependencies
- scVI
- SCANPY
- matplotlib, numpy, seaborn.
*environment and installation files coming soon*.

##### Questions and Troubleshooting:
Please open an [issue](https://github.com/theislab/covid_macrophages_integration/issues).

**License**: MIT.
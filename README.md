### Integration and analysis of Lung macrophage cells from patients with COVID-19, Idiopathic pulmonary fibrosis,  COPD. 
The following notebooks described the data integration, data exploration and *proximity analysis* steps of Methods as implemented in

[**Wendisch, Dietrich, Mari, Stillfried *et al.* (2021). (in press)**](https://www.sciencedirect.com/science/article/pii/S0092867421013830)

Jupyter notebooks located in the directory `notebooks` describe steps for data of external samples and BAL (this study).

### Processed Data
Pre-processed H5AD and the integrated object are zipped and available in [**Dropbox (~15 GB)**](https://www.dropbox.com/s/v2fdj5u4svxmj4j/data.zip?dl=0).
Please update the `data` directory with those.

### Pre-processing and integration
Main processing steps are:
  1. Data Preparation of Adams/Morse/Bharat and our work (filtering, normalization). H5AD files are generated for each dataset using the public raw data, and then deployed in `data`.
  2. Data Integration of all datasets using (scVI), and patient/sample IDs as batches.

### Post integration analyses
  3.1. UMAP visualization and gene expression analyses using macrophage gene groups post-processing.
  3.2. Counting and quantification of link between macrophages from COVID-19 patients and IPF, using distances in the kNN graph from the integrated cells as a reference for mapping conditions (i.e. proximity analysis).

#### Installation
To install environment and relevant dependencies, please clone and execute the following command.
```
conda env create -f environment.yml
```
**Option 1)** Install `jupyter` in the same environment.
```
conda activate covid_macrophages_integration
conda install -c conda-forge jupyterlab
```
**Option 2)** Add kernel to jupyter running in another environment.
```
conda activate covid_macrophages_integration
conda install -c anaconda ipykernel
python -m ipykernel install --name covid_macrophages_integration
```


Troubleshooting:
Please open an [issue](https://github.com/theislab/covid_macrophages_integration/issues).

**License**: MIT.

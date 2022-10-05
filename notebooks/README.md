# Analysis Pipelines

This README describes the analysis pipelines currently implemented by the notebooks in this repo.

## General naming conventions

Notebooks for running analyses are prefixed with a number based on their function, as described below.

`1_` – notebooks for downloading files
`2_` – notebooks for performing embedding of datasets
`3_` – notebooks for analyzing cells of a single species and comparing to the same cells in a new embedding
`4_` – notebooks for analyzing cells of multiple species in the same joint embedding

## Vertebrate adult brain analysis (Drer, Mmus, Xlae) using OrthoFinder
Analysis of cells from three adult brain scRNA-Seq datasets from zebrafish, mouse, and frog (_Xenopus laevis_) using an OrthoFinder embedding.

0. Make sure there's an empty `output/` folder in the home `glial-origins/` directory.  
> This is needed to store downloaded and modified output files.  
> Need to patch this out so it automatically makes it for you.  

1. Run the downloading scripts at:
- `Drer_adultbrain/1_Drer_adultbrain_downloading.ipynb`  
- `Mmus_adultbrain/1_Mmus_adultbrain_downloading.ipynb`  
- `Xlae_adultbrain/1_Xlae_adultbrain_downloading.ipynb` 
> This downloads the three datasets and generates necessary files for the next step in the analysis.  
> Note: You'll need to manually copy the `transdecoder.pep` file into the `output/` folder for the specific species at the end of this analysis.  
> Working on patching this still.

2. Run the OrthoFinder analysis at:  
- `DrerMmusXlae_adultbrain_OrthoFinder/2_DrerMmusXlae_adultbrain_runOrthoFinder.ipynb`  
> This runs OrthoFinder on peptides from the three datasets.

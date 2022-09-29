# Analysis Pipelines

This README describes the analysis pipelines currently implemented by the notebooks in this repo.

## Vertebrate adult brain analysis (Drer, Mmus, Xlae) using OrthoFinder
Analysis of cells from three adult brain scRNA-Seq datasets from zebrafish, mouse, and frog (_Xenopus laevis_) using an OrthoFinder embedding.

1. Run the downloading scripts at:
- `Drer_adultbrain/1_Drer_adultbrain_downloading.ipynb`  
- `Mmus_adultbrain/1_Mmus_adultbrain_downloading.ipynb`  
- `Xlae_adultbrain/1_Xlae_adultbrain_downloading.ipynb` 
> This downloads the three datasets and generates necessary files for the next step in the analysis. 

2. Run the OrthoFinder analysis at:  
- `DrerMmusXlae_adultbrain_OrthoFinder/1_DrerMmusXlae_adultbrain_runOrthoFinder.ipynb`
> This runs OrthoFinder on peptides from the three datasets.

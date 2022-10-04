# Analysis Pipelines

This README describes the analysis pipelines currently implemented by the notebooks in this repo.

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
- `DrerMmusXlae_adultbrain_OrthoFinder/1_DrerMmusXlae_adultbrain_runOrthoFinder.ipynb`  
> This runs OrthoFinder on peptides from the three datasets.

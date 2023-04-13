## Pipeline organization

### General naming conventions
Notebooks for running analyses are prefixed with a number based on their function, as described below.

`1_` – notebooks for downloading genome sequences, annotations, and genes x cells matrices 
    `a_` – notebooks for downloading cell type annotation data
`2_` – notebooks for performing embedding of datasets  
`3_` – notebooks for analyzing cells of a single species and comparing to the same cells in a new embedding  
    `c_` – notebooks for generating plots of gene expression vs. feature abundance
`4_` – notebooks for analyzing cells of multiple species in the same joint embedding  
    `b_` – notebooks for generating interactive plots of joint embedding spaces using Plotly

### Vertebrate adult brain analysis (Drer, Mmus, Xlae) using **OrthoFinder** and **AlphaFold → FoldSeek**
Analysis of cells from three adult brain scRNA-Seq datasets from zebrafish, mouse, and frog (_Xenopus laevis_) using an Orthogroup or Structural scluster embedding.  
Our analysis can be reproduced by running the provided notebooks in the following order.  

0. Make sure there's an empty `output/` folder in the home `glial-origins/` directory.  
> This should automatically be created when you first run a script.

1. Download necessary genome FASTA, genome annotation GFF, and genes x cells matrix TXT files:
- [1_downloading/1_Drer_adultbrain_downloading.ipynb](1_downloading/1_Drer_adultbrain_downloading.ipynb)
- [1_downloading/1_Mmus_adultbrain_downloading.ipynb](1_downloading/1_Mmus_adultbrain_downloading.ipynb)
- [1_downloading/1_Xlae_adultbrain_downloading.ipynb](1_downloading/1_Xlae_adultbrain_downloading.ipynb)
> This downloads the three datasets and generates necessary files for the next steps in the analysis.

2. Download and format cell type annotation files:
- [1_downloading/a_Drer_adultbrain_cellannot_downloading.ipynb](1_downloading/a_Drer_adultbrain_cellannot_downloading.ipynb)
- [1_downloading/a_Mmus_adultbrain_cellannot_downloading.ipynb](1_downloading/a_Mmus_adultbrain_cellannot_downloading.ipynb)
- [1_downloading/a_Xlae_adultbrain_cellannot_downloading.ipynb](1_downloading/a_Xlae_adultbrain_cellannot_downloading.ipynb)
> This downloads cell type annotations for three datasets and generates necessary files for the next steps in the analysis.

3. Run the OrthoFinder and FoldSeek analyses:
- [2_feature-embedding/2_DrerMmusXlae_adultbrain_runOrthoFinder.ipynb](2_feature-embedding/2_DrerMmusXlae_adultbrain_runOrthoFinder.ipynb)
> This runs OrthoFinder on peptides from the three datasets.  

- [2_feature-embedding/2_DrerMmusXlae_adultbrain_runFoldSeek.ipynb](2_feature-embedding/2_DrerMmusXlae_adultbrain_runFoldSeek.ipynb)
> This downloads AlphaFold structures for all available proteins for each species, then clusters structures using FoldSeek.

4. Visualize cells from each species and collect differentially expressed Orthogroups/ Structural clusters:
- [3_single-species-exploration/3_Drer_adultbrain_exploration-FoldSeek.ipynb](3_single-species-exploration/3_Drer_adultbrain_exploration-FoldSeek.ipynb)
- [3_single-species-exploration/3_Drer_adultbrain_exploration-OrthoFinder.ipynb](3_single-species-exploration/3_Drer_adultbrain_exploration-OrthoFinder.ipynb)
- [3_single-species-exploration/3_Mmus_adultbrain_exploration-FoldSeek.ipynb](3_single-species-exploration/3_Mmus_adultbrain_exploration-FoldSeek.ipynb)
- [3_single-species-exploration/3_Mmus_adultbrain_exploration-OrthoFinder.ipynb](3_single-species-exploration/3_Mmus_adultbrain_exploration-OrthoFinder.ipynb)
- [3_single-species-exploration/3_Xlae_adultbrain_exploration-FoldSeek.ipynb](3_single-species-exploration/3_Xlae_adultbrain_exploration-FoldSeek.ipynb)
- [3_single-species-exploration/3_Xlae_adultbrain_exploration-OrthoFinder.ipynb](3_single-species-exploration/3_Xlae_adultbrain_exploration-OrthoFinder.ipynb)
> These notebooks visualize and analyze gene expression/ feature abundance for each specices for each embedding space.  
> The notebooks generate plots for visualization purposes and also generate lists of differentially expressed Orthogroups or Structural clusters for the next step.  

- [3_single-species-exploration/c_Drer_adultbrain_supplementary_plots.ipynb](3_single-species-exploration/c_Drer_adultbrain_supplementary_plots.ipynb)
> This notebook generates plots of gene expression vs. feature abundance for select Zebrafish genes.  
> The notebook should be easily modifiable to show similar plots for the other species in the analysis.  

5. Visualize cells from all three species in joint embedding spaces:
- [4_multi-species-exploration/4_DrerMmusXlae_adultbrain_exploration-Foldseek.ipynb](4_multi-species-exploration/4_DrerMmusXlae_adultbrain_exploration-Foldseek.ipynb)
- [4_multi-species-exploration/4_DrerMmusXlae_adultbrain_exploration-OrthoFinder.ipynb](4_multi-species-exploration/4_DrerMmusXlae_adultbrain_exploration-OrthoFinder.ipynb)
> These notebooks visualize and analyze feature abundance for all three species in either Orthogroup or Structural cluster feature space.  

- [4_multi-species-exploration/b_DrerMmusXlae_FoldSeek_plotly-testing.ipynb](4_multi-species-exploration/b_DrerMmusXlae_FoldSeek_plotly-testing.ipynb)
- [4_multi-species-exploration/b_DrerMmusXlae_FoldSeek_plotly-testing.ipynb](4_multi-species-exploration/b_DrerMmusXlae_FoldSeek_plotly-testing.ipynb)
> These notebooks build Plotly interactive HTML plots for each of the two embedding spaces.  

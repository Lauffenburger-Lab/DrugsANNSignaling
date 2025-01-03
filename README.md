[![DOI](https://zenodo.org/badge/678252886.svg)](https://doi.org/10.5281/zenodo.14057134)
# Inferring Off-Target effects of drugs on cellular signaling using Interactome-Based deep learning
Github repository of the study:
> Inferring Off-Target effects of drugs on cellular signaling using Interactome-Based deep learning <br>
> Nikolaos Meimetis<sup>1</sup>, Douglas A. Lauffenburger<sup>1</sup>, Avlant Nilsson<sup>1,2,3*</sup>
> 1) Department of Biological Engineering, Massachusetts Institute of Technology, Cambridge, MA 02139, USA
> 2) Department of Cell and Molecular Biology, SciLifeLab, Karolinska Institutet, Sweden
> 3) Department of Biology and Biological Engineering, Chalmers University of Technology, Gothenburg, SE 41296, Sweden
> * Corresponding author, avlant.nilsson@ki.se

doi: https://doi.org/10.1016/j.isci.2024.109509

This repository is administered by @NickMeim. For questions contact meimetis@mit.edu

**Trained models of this study are too big to be uploaded here and are available upon reasonable request. Supplementary Data File 1.xlsx and Supplementary Data File 2.xlsx is in the results folder**

**Ensembles of 50 models are trained for 33 cell lines in the L1000 dataset and are available here:**
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14057298.svg)](https://doi.org/10.5281/zenodo.14057298)

Many diseases emerge from dysregulated cellular signaling, and drugs are often designed to target specific nodes in cellular networks e.g. signaling proteins, or transcription factors. However, off-target effects are common and may ultimately result in failed clinical trials. Computational modeling of the cell’s transcriptional response to drugs could improve our understanding of their mechanisms of action. Here we develop such an approach based on ensembles of artificial neural networks, that simultaneously infer drug-target interactions and their downstream effects on intracellular signaling. Applied to gene expression data from different cell lines, it outperforms basic machine learning approaches in predicting transcription factors’ activity, while recovering most known drug-target interactions and inferring many new, which we validate in an independent dataset. As a case study, we explore the inferred interactions of the drug Lestaurtinib and its effects on downstream signaling. Beyond its intended target (FLT3) the model predicts an inhibition of CDK2 that enhances downregulation of the cell cycle-critical transcription factor FOXM1, corroborating literature findings. Our approach can therefore enhance our understanding of drug signaling for therapeutic design.

The current repository contains code for:
1. Initial evaluation of the quality and preprocessing of the data.
2. Training and fitting of ANN and other models.
3. Evaluation of the predictions of various models.
4. Network construction of the MoA of off-target effects of drugs.
5. Drug-target interaction inference.
6. Code to re-create the results of the research article.

## User case studies
To run your own case study follow the instructions in each folder (there are user friendly scripts explained in the README files of each folder) :
1. First visit the preprocessing folder.
2. Then visit the learning folder.
3. Then visit the postprocessing folder.
4. Finally visit the MoA folder.

## Data
The transcriptomic signatures (level 3 and level 5 profiles) of the L1000 CMap resource[^1] are used for this study, together with data from the Bioconductor resource[^2].

The transcriptomic profiles were generated by measuring 978 important (landmark) genes in cancer with a Luminex bead-based assay and computationally inferring the rest[^1]. 

**Details on how to access these data can be found in the data folder**, but generally the main resources can be accessed in GEO: [GSE92742](https://www-ncbi-nlm-nih-gov.libproxy.mit.edu/geo/query/acc.cgi?acc=GSE92742)

## Folder structure
1. **article_supplementary_info** : Folder containing code to re-create the supplementary figures and tables of the article
2. **data** : Folder that should contain the retrieved raw data of the study.
3. **figures** : Folder containing the scripts to produce the figures of the study.
5. **learning** : Folder containing deep learning and machine learning algorithms and models.
6. **preprocessing** : Folder containing scripts to pre-process the raw data and evaluate their quality.
	* **preprocessed_data** : Here the pre-processed data to be used in the subsequent analysis are stored.
7. **results** : Here the results of a subsequent analysis should be stored. Here you can also find all the inferred interactions in Supplementary Data File 1.xlsx 
8. **postprocessing** : Folder containing scripts to evaluate models' results and predictions.
9. **MoA** : Folder containing code and data to construct the MoA of off-target effects.

## Installation
The study utilizes multiple resources from the Python and R programming languages.

**Important Note:**
* **This installation has been validated to work in Unix-based, macOS, and WINDOWS operating systems.** 
* **For a Linux installation there might be needed some manual installation of external dependencies (especially) for tidyverse. Please check libraries' documentation online**
* **Please note that macOS are not compatible with the GPU components of this installation guide (which are not necessary though!).**

**Python installation**
```bash
# After installing anaconda create a conda environment:
conda create -n DTLembas
conda activate DTLembas
conda install -c conda-forge rdkit
conda install -c conda-forge scikit-learn 
pip install networkx
# For general (CPU) pytorch version run the following
# Otherwise for GPU installation run for your own cuda version this command: conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
conda install pytorch torchvision torchaudio -c pytorch
conda install captum -c pytorch
```

**R installation**
Install R studio, open it, and run:
```bash
> if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
> BiocManager::install(c("cmapR","rhdf5","dorothea","org.Hs.eg.db","hgu133a.db"))
> if (!require("tidyverse", quietly = TRUE))
install.packages("tidyverse")
> if (!require("ggplot2", quietly = TRUE))
install.packages("ggplot2")
> install.packages("ggrepel")
> install.packages("ggpubr")
> install.packages("doRNG")
> install.packages("doFuture")
```
Alternatively, use conda and always use R from the terminal:
```bash
conda create -n DTLembas_r_env
conda activate DTLembas_r_env
conda install -c r r-essentials
conda install r-BiocManager
conda install conda-forge::r-ggrepel
conda install r-ggpubr
conda install r-doRNG
conda install r-doFuture
R()
BiocManager::install(c("cmapR","rhdf5","dorothea","org.Hs.eg.db","hgu133a.db"))
```

**R dependencies**: 
You can check the list below and manually install your preferences.

In a quick overview, the following R libraries and versions (**although any version of the following libraries is appropriate**) were/are used to produce the figures and results of the study:
1. [R](https://cran.r-project.org/bin/windows/base/) version 4.1.2
2. [tidyverse](https://www.tidyverse.org/packages/) 1.3.1
3. [BiocManager](https://www.bioconductor.org/install/) 1.30.16
4. [cmapR](https://bioconductor.org/packages/release/bioc/html/cmapR.html) 1.4.0
5. [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) 3.13.0
6. [rhdf5](https://bioconductor.org/packages/release/bioc/html/rhdf5.html) 2.36.0
7. [doFuture](https://cran.r-project.org/web/packages/doFuture/index.html) 0.12.0
8. [doRNG](https://cran.r-project.org/web/packages/doRNG/index.html) 1.8.2
9. [ggplot2](https://ggplot2.tidyverse.org/) 3.3.5
10. [ggpubr](https://www.rdocumentation.org/packages/ggpubr/versions/0.4.0) 0.4.0
11. [GeneExpressionSignature](https://www.bioconductor.org/packages/release/bioc/html/GeneExpressionSignature.html) 1.38.0
12. [caret](https://cran.r-project.org/web/packages/caret/index.html) 6.0-94
13. [ggpubr](https://rpkgs.datanovia.com/ggpubr/) 0.6.0
14. [ggpattern](https://coolbutuseless.github.io/package/ggpattern/) 1.1.0
15. [ggridges](https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html) 0.5.4
16. [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html) 0.9.3
17. [rstatix](https://cran.r-project.org/web/packages/rstatix/index.html) 0.7.2
18. [patchwork](https://patchwork.data-imaginist.com/) 1.1.2.9000
19. [dorothea](https://saezlab.github.io/dorothea/) 1.4.2
20. [AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html) 1.54.1
21. [PharmacoGx](https://bioconductor.org/packages/release/bioc/html/PharmacoGx.html) 2.4.0
22. [GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html) 2.60.0
23. [hgu133a.db](https://bioconductor.org/packages/release/data/annotation/html/hgu133a.db.html) 3.13.0
24. [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) 3.48.3
25. [affy](https://www.bioconductor.org/packages/release/bioc/html/affy.html) 1.70.0
26. [dbparser](https://cran.r-project.org/web/packages/dbparser/vignettes/dbparser.html) 2.0.1

**Python dependencies**: 
First, install conda (anaconda) environment on your computer, and then you can use the commands **in a bash-terminal** after the list of libraries.

In a quick overview, the following Python libraries and versions (**although different versions are POSSIBLY also appropriate**) were/are used:
1. [python](https://www.python.org/downloads/) 3.8.8
2. [seaborn](https://seaborn.pydata.org/installing.html) 0.11.2 (version does not matter for this library)
3. [numpy](https://numpy.org/install/) 1.20.3 (version does not matter for this library)
4. [pandas](https://pandas.pydata.org/docs/getting_started/install.html) 1.3.5 (version does not matter for this library)
5. [matplotlib](https://anaconda.org/conda-forge/matplotlib) 3.5.1 (version does not matter for this library)
6. [scipy](https://anaconda.org/anaconda/scipy) 1.7.3
7. [scikit-learn](https://scikit-learn.org/stable/install.html) 1.0.2
8. [networkx](https://networkx.org/documentation/stable/install.html) 2.6.3
9. [rdkit](https://www.rdkit.org/docs/index.html) 2021.03.5
10. [captum](https://captum.ai/docs/getting_started) 0.5.0
11. [pytorch](https://pytorch.org/get-started/locally/) 1.12.0


## References
[^1]: Subramanian, Aravind, et al. "A next generation connectivity map: L1000 platform and the first 1,000,000 profiles." Cell 171.6 (2017): 1437-1452.
[^2]: Gentleman, Robert C., et al. "Bioconductor: open software development for computational biology and bioinformatics." Genome biology 5.10 (2004): 1-16.

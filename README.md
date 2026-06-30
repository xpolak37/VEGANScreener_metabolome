# VEGANScreener: metabolomics evaluation

<div style="font-size: larger;">
Petra Polakovicova,  Anna Ouradova, Helena Pelantova,  ...., Selma Kronsteiner-Gicevic, Jan Gojda, Marek Kuzma, Monika Cahova
</div>

## General information

This repository provides a comprehensive report of the study **The serum metabolome reflects VEGANScreener-assessed diet quality across five European vegan populations**

All reported results can be reproduced using the code in this repository. Feel free to contact Petra Polakovicova by [petra.polakovicova@ikem.cz](petra.polakovicova@ikem.cz) if you have any questions about the computational part of the study.

📚 **Citation** 

If you find this code and report helpful, cite the original publication:

(TO BE ADDED)

💾 **Data Availability**

Metabolomic data for this study have been uploaded to MetaboLights database under accession number (TO BE ADDED).


## Report info

This project analyses serum metabolomic data from five European cohorts (Belgium, Switzerland, Czech Republic, Germany, Spain) as part of an European project [VEGANScreener](https://www.veganscreener.eu/). 357 subjects are included who completed the VEGANScreener, the general questionnaire, and the 4-day prospective diet record. 

The analysis is divided into 5 main questions:

- **Q0 – General**: What are the main demographic, anthropometric, and dietary determinants of variation in the serum metabolome across the cohort?
- **Q1 – Country**: To what extent does country of residence shape the serum metabolome, and can metabolomic profiles discriminate participants by country?
- **Q2 – VEGANScreener**: Does the VEGANScreener score capture a biologically meaningful, metabolomically detectable signal of plant-based diet quality, independent of other covariates?
- **Q3 – PDI**: Is the overall plant-based diet index (PDI) associated with the serum metabolome, and does this association hold after adjusting for covariates?
- **Q4 – UPF**: Is ultra-processed food intake (%TE) associated with distinct serum metabolomic signatures, and how does its predictive power compare to PDI and Veganscreener score?
- **Q5** – Which individual PDI food groups (e.g., whole grains, fruits, vegetables ..) are associated with specific serum metabolites, and do these associations reveal which dietary components drive the overall PDI–metabolome relationship?

**Project Structure**

Below is an overview of the folder structure:

- **analysis/**  
  - `scripts/` - scripts for data analysis
	- `preprocessing/` – data cleaning and integration of LC-MS and NMR datasets
	- `custom_functions_vegans.R` – shared helper functions used across analysis scripts
	- `VSmetabo_Figures.qmd` – generates main and supplementary figures
	- `VSmetabo_overfitting_check.qmd` – label-shuffling check for model overfitting
	- `VSmetabo_Q0.qmd` – metabolome vs. general covariates
	- `VSmetabo_Q1.qmd` – metabolome vs. country
	- `VSmetabo_Q1h.qmd` – metabolome vs. country, refined sub-analysis
	- `VSmetabo_Q2.qmd` – metabolome vs. VEGANScreener score
	- `VSmetabo_Q2minus.qmd` – metabolome vs. VGscore negative items
	- `VSmetabo_Q2plus.qmd` – metabolome vs. VGscore positive items
	- `VSmetabo_Q3_PDI.qmd` – metabolome vs. overall PDI
	- `VSmetabo_Q3_hPDI.qmd` – metabolome vs. healthful PDI
	- `VSmetabo_Q3_uDI.qmd` – metabolome vs. unhealthful PDI
	- `VSmetabo_Q4_UPF.qmd` – metabolome vs. ultra-processed food intake
	- `VSmetabo_Q5_PDIgroups.qmd` – metabolites vs. individual PDI food groups

  - `results` – results generated directly via provided scripts

## Methodology

For detailed methodology, see the original publication. 

**Data collection** 

Dietary quality was assessed using the VEGANScreener, a 29-item questionnaire scoring key plant-based food groups and supplements on a traffic-light system (total range 0–66), alongside plant-based diet indices (overall PDI, healthful hPDI, unhealthful uPDI) and ultra-processed food intake (%TE, NOVA classification), all derived from up to four non-consecutive diet records collected via the Nutrixo platform and energy-adjusted using the residual method. Serum metabolomics was performed using two complementary platforms: untargeted LC-MS (HPLC separation with amino acid column, high-resolution MS in positive/negative mode, metabolite identification against HMDB, MoNA, and MassBank) and untargeted NMR (CPMG experiments on a 600 MHz spectrometer, identification via Chenomx and HMDB). The two metabolomic datasets were preprocessed separately (filtering, imputation, inverse normal transformation, scaling), then integrated either through sparse PCA (sPCA, first 20 components) or by simple concatenation.

**Statistical analysis**

Statistical analyses were conducted in R, combining exploratory correlation analysis (Spearman, PCs vs. covariates), multivariate testing (PERMANOVA on Euclidean distances, adjusting for sex, BMI, age, dietary scores, and country), and univariate linear models per metabolite (FDR-corrected), visualized through bipartite networks and forest plots. A metabolite-based index was built via cross-validated linear regression using metabolites significant in univariate analysis, correlated against each dietary score. Random forest models (ranger package, bootstrap resampling with OOB validation, ROC/AUC via pROC) assessed the predictive power of metabolomic principal components for dietary scores and country of residence, with permutation testing confirming models were not overfit (shuffled-label AUC ≤0.506).

## Results 

The code with reported results can be found:

**Exploratory analysis**:
- [Q0_analysis](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q0.html): metabolome vs. general covariates

**Country-related associations**:
- [Q1_analysis](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q1.html): metabolome vs. country of residence

**Veganscreener evaluation**:
- [Q2_analysis VEGANSCREENER](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q2.html): metabolome vs. overall VEGANScreener score
- [Q2_analysis VEGANSCREENER (+)](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q2plus.html): metabolome vs. VGscore positive items
- [Q2_analysis VEGANSCREENER (-)](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q2minus.html): metabolome vs. VGscore negative items

**PDI evaluation**:
- [Q3_analysis PDI](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q3_PDI.html): metabolome vs. overall PDI
- [Q3_analysis hPDI](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q3_hPDI.html): metabolome vs. healthful PDI
- [Q3_analysis uPDI](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q3_uPDI.html): metabolome vs. unhealthful PDI

**UPF evaluation**:
- [Q4_analysis](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q4_UPF.html): metabolome vs. ultra-processed food intake

**PDI groups**:
- [Q5_analysis](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q5_PDIgroups.html): metabolites vs. individual PDI food groups

**ML overfitting check**:
- [ML_overfitting_check](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_overfitting_check.html): label-shuffling validation of RF model performance

**Figures included in the original publication**:
- [Figures](https://github.com/xpolak37/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Figures.html): code generating the manuscript's main and supplementary figures


---------------------------------------------------------------------------------------------------

## Acknowledgment

(TO BE ADDED)



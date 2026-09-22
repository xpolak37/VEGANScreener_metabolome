# VEGANScreener: metabolomics evaluation

<div style="font-size: larger;">
Petra Polakovicova<sup>1,2,3†</sup>, Anna Ouradova<sup>4†</sup>, Doris Janousova<sup>3,5</sup>, Dominika Sindelarova<sup>3,6</sup>,
Helena Pelantova<sup>3</sup>, Selma Kronsteiner-Gicevic<sup>7</sup>, Jan Gojda<sup>4†</sup>, Marek Kuzma<sup>3†</sup>, Markus Keller<sup>8†</sup>,
Leonie H. Bog<sup>l9†</sup>, Isabelle Herter-Aeberli<sup>10†</sup>, Stefaan De Henauw<sup>11†</sup>, Maira Bes-Rastrollo<sup>12,13,14,15†</sup>, Eva S. Schernhammer<sup>7,16†</sup>, Monika Cahova<sup>1†</sup>

</div>

<br>

<sup>†</sup> These authors have contributed equally to this work and share first/last authorship 

<sup>1</sup> Institute for Clinical and Experimental Medicine, Department of Hepatogastroenterology, Prague, CR, Czech Republic  
<sup>2</sup> Faculty of Science, Charles University, Prague, Czech Republic  
<sup>3</sup> Institute of Microbiology, The Czech Academy of Sciences, Prague, Czech Republic  
<sup>4</sup> Department of Internal Medicine, University Hospital Kralovske Vinohrady, Third Faculty of Medicine, Charles University, Prague, Czech Republic  
<sup>5</sup> Faculty of Food and Biochemical Technology, University of Chemistry and Technology, Prague, Czech Republic  
<sup>6</sup> Faculty of Chemical Technology, University of Chemistry and Technology, Prague, Czech Republic  
<sup>7</sup> Department of Epidemiology, Center for Public Health, Medical University of Vienna, Vienna, Austria  
<sup>8</sup> Research Institute for Plant-Based Nutrition, Giessen, Germany  
<sup>9</sup> Department of Health Professions, Division of Nutrition and Dietetics, Bern University of Applied Sciences, Bern, Switzerland  
<sup>10</sup> Laboratory of Nutrition and Metabolic Epigenetics, Institute of Food, Nutrition and Health, ETH Zürich, Zürich, Switzerland  
<sup>11</sup> Department of Public Health and Primary Care, Faculty of Medicine and Health Sciences, Ghent University, Ghent, Belgium  
<sup>12</sup> Department of Preventive Medicine and Public Health, University of Navarra, Pamplona, Spain  	
<sup>13</sup> Institute for Nutrition and Health (INS), University of Navarra, Pamplona, Spain  
<sup>14</sup> Navarra Health Research Institute (IdiSNA), Pamplona, Spain  
<sup>15</sup> CIBERobn, Instituto de Salud Carlos III, Madrid, Spain  
<sup>16</sup> Department of Epidemiology, Harvard T.H. Chan School of Public Health, Boston, MA, USA  

---------------------------------------------------------------------------------------------------

## General information

This repository provides a comprehensive report of the study **The serum metabolome reflects VEGANScreener-assessed diet quality in vegans from five European countries**

All reported results can be reproduced using the code in this repository. Feel free to contact Petra Polakovicova by [petra.polakovicova@ikem.cz](petra.polakovicova@ikem.cz) if you have any questions about the computational part of the study.

📚 **Citation** 

If you find this code and report helpful, cite the original publication:

(TO BE ADDED)

💾 **Data Availability**

Metabolomic data for this study data have been deposited in the ASEP repository with the dataset
identifier xxx (TO BE ADDED). 

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

## Methodology

For detailed methodology, see the original publication. 

**Data collection** 

Dietary quality was assessed using the VEGANScreener, a 29-item questionnaire scoring key plant-based food groups and supplements on a traffic-light system (total range 0–66), alongside plant-based diet indices (overall PDI, healthful hPDI, unhealthful uPDI) and ultra-processed food intake (%TE, NOVA classification), all derived from up to four non-consecutive diet records collected via the Nutrixo platform and energy-adjusted using the residual method. Serum metabolomics was performed using two complementary platforms: untargeted LC-MS (HPLC separation with amino acid column, high-resolution MS in positive/negative mode, metabolite identification against HMDB, MoNA, and MassBank) and untargeted NMR (CPMG experiments on a 600 MHz spectrometer, identification via Chenomx and HMDB). The two metabolomic datasets were preprocessed separately (filtering, imputation, inverse normal transformation, scaling), then integrated either through sparse PCA (sPCA, first 20 components) or by simple concatenation.

**Statistical analysis**

Statistical analyses were conducted in R, combining exploratory correlation analysis (Spearman, PCs vs. covariates), multivariate testing (PERMANOVA on Euclidean distances, adjusting for sex, BMI, age, dietary scores, and country), and univariate linear models per metabolite (FDR-corrected), visualized through bipartite networks and forest plots. A metabolite-based index was built via cross-validated linear regression using metabolites significant in univariate analysis, correlated against each dietary score. Random forest models (ranger package, bootstrap resampling with OOB validation, ROC/AUC via pROC) assessed the predictive power of metabolomic principal components for dietary scores and country of residence, with permutation testing confirming models were not overfit (shuffled-label AUC ≤0.506).

## Results 

The code with reported results can be found:

**Exploratory analysis**:
- [Q0_analysis](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q0.html): metabolome vs. general covariates

**Country-related associations**:
- [Q1_analysis](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q1.html): metabolome vs. country of residence

**Veganscreener evaluation**:
- [Q2_analysis VEGANSCREENER](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q2.html): metabolome vs. overall VEGANScreener score
- [Q2_analysis VEGANSCREENER (+)](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q2plus.html): metabolome vs. VGscore positive items
- [Q2_analysis VEGANSCREENER (-)](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q2minus.html): metabolome vs. VGscore negative items

**PDI evaluation**:
- [Q3_analysis PDI](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q3_PDI.html): metabolome vs. overall PDI
- [Q3_analysis hPDI](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q3_hPDI.html): metabolome vs. healthful PDI
- [Q3_analysis uPDI](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q3_uPDI.html): metabolome vs. unhealthful PDI

**UPF evaluation**:
- [Q4_analysis](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q4_UPF.html): metabolome vs. ultra-processed food intake

**PDI groups**:
- [Q5_analysis](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Q5_PDIgroups.html): metabolites vs. individual PDI food groups

**ML overfitting check**:
- [ML_overfitting_check](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_overfitting_check.html): label-shuffling validation of RF model performance

**Figures included in the original publication**:
- [Figures](https://xpolak37.github.io/VEGANScreener_metabolome/analysis/scripts/VSmetabo_Figures.html): code generating the manuscript's main and supplementary figures

---------------------------------------------------------------------------------------------------

## Funding

HDHL-INTIMIC: Standardized measurement, monitoring and/or biomarkers to study food intake, physical activity and health (STAMIFY 2021) in 2022–25. The VEGANScreener study is supported by ERA-Net HDHL-INTIMIC, the European Union’s Horizon 2020 Research and Innovation Programme under grant agreement No. 727565. The following institutes provide support: the Austrian Research Promotion Agency FFG (Austria; project no. FO999890542) and the Austrian Federal Ministry of Education, Science and Research (BMBWF); the Research Foundation Flanders FWO (Belgium, project no. G0G5121N); the Czech Ministry of Education, Youth, and Sports (Czech Republic, project no. 8F22003, Numbers of contracts: MSMT-88/2021-29/2 and MSMT-88/2021-29/3); MH CZ - DRO („Institute for Clinical and Experimental Medicine – IKEM, IN 00023001“); the Ministry of Health of the Czech Republic in cooperation with the Czech Health Research Council under project No. NW26A-CARDIA; Federal Ministry of Education and Research (BMBF) (Germany, project no. 01EA2202), and the Spanish Ministry of Science, Innovation and Universities within the framework of Next Generation EU (Spain, project no. AC21\_2/00015). Occident Foundation additionally supported the Spanish team (Research Award 2023).  The Swiss arm of the VEGANScreener study is supported by the Federal Food Safety and Veterinary Office (FSVO), Vontobel Foundation and Foundation for the encouragement of Nutrition Research in Switzerland (SFEFS).



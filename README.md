## HPV and p16INK4a in Penile Squamous Cell Carcinomas

### Description
This dataset contains information on 59 patients with penile lesions, including penile intraepitheial neoplasia (PeIN) and squamous cell carcinoma (SCC). Formalin-fixed and paraffin-embedded tissue samples were used to build 2 tissue microarrays (TMA). Clinicopathologic and outcome information was obtained from medical records. This repository contains the following files, pages, and folders:

* [Results](https://github.com/alcideschaux/PenisHPVp16/blob/master/hpvp16penis.md) This page contains the results of the analysis.
* [Codebook](https://github.com/alcideschaux/PenisHPVp16/blob/master/CODEBOOK.md) This page contains the description of the variables with labels and levels.
* [Dataset](https://github.com/alcideschaux/PenisHPVp16/blob/master/PenisDataset.csv) This file corresponds to the dataset that it is used for analysis.
* [Figures](https://github.com/alcideschaux/PenisHPVp16/tree/master/figure) This folder contains the figures (survival curves) in PNG format.

### Aim
The aim of this project is to evaluate human papillomavirus (HPV) status and p16INK4a immunohistochemical expression, accessing the relationship between HPV and p16INK4a and clinicopathologic and outcome features, and building prognostic models of penile cancer recurrence, progression, and outcome.

### Material & Method
#### Patient cohort
Fifty-nine in situ and invasive penile squamous cell carcinomas from fifty-three different patients diagnosed between 1985 and 2013 were retrieved from our surgical pathology archives. Histologic subtyping was carried out in whole tissue sections using the morphologic criteria presented in the Atlas of Tumor Pathology of the Armed Forces Institute of Pathology. Other pathologic features of the tumor were also re-assessed to include tumor grade, tumor thickness, presence of angiolymphatic and perineural invasion, extent of invasion, urethral involvement, presence of metastasis, pathologic and clinical stage.

Representative formalin-fixed, paraffin-embedded archival blocks were used for the construction of two high-density tissue microarrays (TMAs). Tumors and paired nonneoplastic tissue were spotted 3 to 5 times each using 1.6-mm cores.

Electronic medical records were reviewed for pertinent clinical information including age, gender, anatomical site of involvement, date of surgery, type of treatment, recurrence, progression, death status, cause of death and last follow up date. Recurrence was defined as development of a new tumor in the same location as the prior tumor (local recurrence). Progression was defined as the presence of local or distant metastasis after the primary treatment.

#### Immunohistochemistry for p16INK4a expression
Immunohistochemistry for p16INK4a was performed using the automated Ventana XT system (Ventana Medical Systems, Inc., Tucson, AZ). The reaction was developed using streptavidin-HRP detection I-View kit (Ventana Medical Systems). All sections were then counterstained with hematoxylin, dehydrated, and cover-slipped. Expression for p16INK4a was categorized in each TMA spot as: I) negative p16INK4a expression (complete absence of p16INK4a staining or only focal positivity); II) positive p16INK4a expression (strong and diffuse nuclear and cytoplasmic p16INK4a positivity in majority of tumor cells). A case was classified as “p16INK4a overexpressed” if it showed at least 1 TMA spot with positive p16INK4a expression.

#### High-risk HPV detection by in situ hybridization
HPV in situ hybridization was performed using the automated benchmark XT system (Ventana Medical Systems, Inc.). INFORM HPV III family 16 probe (B) cocktail, with affinity for high-risk (HR) genotypes (16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, and 66), was applied and the reaction was developed using hybrid ready detection kit (reference number 780-4295 and reference number 254-2210-2116, respectively, Ventana Medical Systems, Inc.). A case was considered positive when unequivocal punctuated blue nuclear staining was observed in tumor cells.

#### Statistical analysis
Associations between variables were evaluated using the Mann-Whitney U test and the Fisher's exact test. Predictor variables were grouped in clinical features (age and race), pathologic features (histologic subtype, anatomical site, anatomical level, histologic grade, tumor thickness, tumor invasion of penile urethra, lymphovascular invasion, perineural invasion, pT stage, pN stage, clinical stage, and inguinal lymph node metastasis), HPV status, and p16INK4a overexpression. For time-to-event analysis endpoints included tumor recurrence, tumor progression, overall mortality, and cancer-related mortality.

Survival curves were built using the Kaplan-Meier method and compared using the log-rank test. Odds ratios and hazard ratios were estimated using unconditional logistic regression and Cox's proportional hazard regression, respectively, for each of the aforementioned endpoints. A 2-tailed P < 0.05 was required for statistical significance.

Data were analyzed using R Version 3.1.1 “Sock it to Me” (R Foundation for Statistical Computing, Vienna, Austria).

### Researchers
* Stephania M. Bezerra<sup>1</sup>
* Alcides Chaux<sup>1,2</sup>
* Mark W. Ball<sup>3</sup>
* Sheila F. Faraj<sup>1</sup>
* Enrico Munari<sup>1</sup>
* Nilda Gonzalez-Roibon<sup>1</sup>
* Rajni Sharma<sup>1</sup>
* Trinity J. Bivalacqua<sup>3</sup>
* Arthur L. Burnett<sup>3</sup>
* George J. Netto<sup>1,2,3</sup>

<sup>1</sup>Department of Pathology, The Johns Hopkins Medical Institutions, Baltimore, MD 21231  
<sup>2</sup>Norte University, Office of Scientific Research, Asunción, Paraguay  
<sup>3</sup>Department of Urology, The Johns Hopkins Medical Institutions, Baltimore, MD 21231  
<sup>4</sup>Department of Oncology, The Johns Hopkins Medical Institutions, Baltimore, MD 21231
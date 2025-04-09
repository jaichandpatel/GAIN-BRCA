# GAIN-BRCA: A Graph-Based Explainable AI Framework for Breast Cancer Subtype Classification Based on Multi-Omics

## Overview

•	GAIN-BRCA leverages multi-omics datasets to classify TCGA breast cancer subtypes based on PAM50.
•	The framework integrates gene expression, DNA methylation, and microRNA data by utilizing their biological interactions. 
•	Specifically, it considers the regulatory relationships where gene expression is influenced by methylation and miRNA interactions. 
•	This graph-based explainable AI framework combines these interactions for improved multi-omics integration and prediction.

## GAIN-BRCA Workflow

![image](https://github.com/user-attachments/assets/f7a5bedb-94e3-47f9-bc6e-584084d4704a)

## File Structure
GAIN-BRCA/ ├── GAIN_BRCA.py
└── dataset/ 
└── Input/ 
├── mRNA_NormCount.csv 
├── miRNA_NormCount.csv 
├── methyl_NormBeta.csv 
├── miRNA_mRNA_interaction.csv 
├── CpG_mRNA_interaction.csv 
└── dependent_variables.csv
├── delta_integration.py 
├── GAIN_BRCA_ANN.py 

Requirements
•	Tensorflow (>= 2.9.1)
•	Keras (>= 2.9.0)
•	Python (>= 3.8.10)
•	Scikit-learn (>=1.2.1)
•	Basic Python packages: numpy, pandas, math
Ensure these libraries are installed.
Data sources: 
1.	Clone or Download the Repository:
Clone this repository or download the source code files.
2.	Prepare the Input Data:
Place the breast cancer multi-omics datasets in the dataset/Input/ folder.
The required files are:
o	mRNA_NormCount.csv
o	miRNA_NormCount.csv
o	methyl_NormBeta.csv
o	miRNA_mRNA_interaction.csv
o	CpG_mRNA_interaction.csv
o	dependent_variables.csv
The datasets can be found here “https://zenodo.org/records/15175435”
3.	Edit File Paths (if necessary):
Ensure that the directory path in GAIN_BRCA.py matches the location of your input files.
4.	Run the Main Script:
Execute the following command from the root directory:
python GAIN_BRCA.py
The script will:
•	Import and integrate the multi-omics datasets.
•	Compute integrated expression values using the provided interaction data.
•	Train an ANN model using stratified k-fold cross-validation.
•	Output the final prediction accuracy and save prediction probabilities to a CSV file.
Code Documentation
•	GAIN_BRCA.py:
Serves as the main execution script, orchestrating data import, integration, and model training.
•	delta_integration.py:
Contains functions to integrate genomic data:
•	int_Delta1(miRNA, mRNA, Int1): Integrates miRNA and mRNA expression data.
•	int_Delta2(methyl, mRNA, Int2): Integrates methylation data and mRNA expression data.
•	GAIN_BRCA_ANN.py:
Implements the ANN model using TensorFlow/Keras.
This module handles data scaling, cross-validation, model training, evaluation, and saving prediction outputs.

Contact
Please contact japatel@unmc.edu if you have any inquiries or issues.


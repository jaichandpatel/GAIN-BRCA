<a name="table-of-contents"></a>
## Table of Contents

- [Overview](#overview)
- [Workflow](#workflow)
- [File Structure](#file-structure)
- [Requirements](#requirements)
- [Data Sources](#data-sources)
  - [Clone or Download](#clone-or-download)
  - [Prepare Input Data](#prepare-input-data)
  - [Configure Paths](#configure-paths)
  - [Run the Pipeline](#run-the-pipeline)
- [Code Documentation](#code-documentation)
- [Contact](#contact)

---

<a name="overview"></a>
## Overview

- **GAIN-BRCA** leverages multi‑omics datasets to classify TCGA breast cancer subtypes based on PAM50.
- Integrates gene expression, DNA methylation, and microRNA data via biological interactions.
- Models regulatory relationships where gene expression is influenced by methylation and miRNA.
- Employs a graph‑based explainable AI framework for enhanced integration and prediction.

<a name="workflow"></a>
## Workflow

![image](https://github.com/user-attachments/assets/f7a5bedb-94e3-47f9-bc6e-584084d4704a)
*Figure: GAIN-BRCA analysis pipeline.*

<a name="file-structure"></a>
## File Structure

```bash
GAIN-BRCA/
├── GAIN_BRCA.py
├── delta_integration.py
├── GAIN_BRCA_ANN.py
└── dataset/
    └── Input/
        ├── mRNA_NormCount.csv
        ├── miRNA_NormCount.csv
        ├── methyl_NormBeta.csv
        ├── miRNA_mRNA_interaction.csv
        ├── CpG_mRNA_interaction.csv
        └── dependent_variables.csv
```

<a name="requirements"></a>
## Requirements

- Python ≥ 3.8.10  
- TensorFlow ≥ 2.9.1  
- Keras ≥ 2.9.0  
- scikit-learn ≥ 1.2.1  
- numpy  
- pandas  
- math

<a name="data-sources"></a>
## Data Sources

<a name="clone-or-download"></a>
### Clone or Download
```bash
git clone https://github.com/GudaLab/GAIN-BRCA
```

<a name="prepare-input-data"></a>
### Prepare Input Data

Place the following CSV files into `dataset/Input/`:

- `mRNA_NormCount.csv`
- `miRNA_NormCount.csv`
- `methyl_NormBeta.csv`
- `miRNA_mRNA_interaction.csv`
- `CpG_mRNA_interaction.csv`
- `dependent_variables.csv`

Datasets available at: https://zenodo.org/records/15175435

<a name="configure-paths"></a>
### Configure Paths

Edit file paths in `GAIN_BRCA.py` to match your local directory structure.

<a name="run-the-pipeline"></a>
### Run the Pipeline
```bash
python GAIN_BRCA.py
```
This will:

- Import and integrate multi‑omics datasets.  
- Compute integrated expression values using interaction data.  
- Train an ANN model using stratified k‑fold cross-validation.  
- Output prediction accuracy and save probabilities to a CSV file.

<a name="code-documentation"></a>
## Code Documentation

- **GAIN_BRCA.py**  
  Main script orchestrating data import, integration, and model training.
- **delta_integration.py**  
  - `int_Delta1(miRNA, mRNA, Int1)`: Integrates miRNA and mRNA expression.  
  - `int_Delta2(methyl, mRNA, Int2)`: Integrates methylation data with mRNA expression.
- **GAIN_BRCA_ANN.py**  
  Implements the ANN model using TensorFlow/Keras: data scaling, cross-validation, training, evaluation, and saving outputs.

<a name="contact"></a>
## Contact

For questions or issues, please contact:  
Jai Patel — japatel@unmc.edu

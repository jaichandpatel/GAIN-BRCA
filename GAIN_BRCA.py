"""
# ================================
Module: GAIN_BRCA
# Authors: Jai Chand Patel, Sushil Kumar Shakyawar, Sahil Sethi, Chittibabu Guda
Description:
    This script performs data integration and prediction for breast cancer (BRCA) subtype analysis.
    It imports mRNA, miRNA, methylation, and interaction data from the specified directory, 
    integrates the data using functions from the delta_integration module, and then uses an 
    ANN model (from the GAIN_BRCA_ANN module) to predict the outcome. The final prediction 
    accuracy is printed to the console.
    # =================================
"""

# --- Import Required Libraries ---
import pandas as pd
import math
import numpy as np
import os

# --- Import Custom Modules ---
from GAIN_BRCA_modules import delta_integration         # Handles interaction-based integration
from GAIN_BRCA_modules import GAIN_BRCA_ANN				# Contains ANN training and evaluation pipeline

# Set the number of loops for cross-validation (recommended value is between 50 and 100)
n_loops = 50  # Recommended value is 50-100

# --- Set Working Directory ---
# Ensure this path points to the folder containing the input CSV files
os.chdir("dataset/Input/")

# --- Load Input Data ---
# Gene expression profiles and interaction networks
# Rounding numerical values to 3 decimal places for uniformity
mRNA = round(pd.read_csv(r"mRNA_NormCount.csv", header=0, index_col=0), 3)
miRNA = round(pd.read_csv(r"miRNA_NormCount.csv", header=0, index_col=0), 3)
methyl = round(pd.read_csv(r"methyl_NormBeta.csv", header=0, index_col=0), 3)

# Interaction networks between omics layers
interaction1 = pd.read_csv(r"miRNA_mRNA_interaction.csv", header=0, index_col=0)
interaction2 = pd.read_csv(r"CpG_mRNA_interaction.csv", header=0, index_col=0)
print("data_upload")

# --- Integration Step ---
# Generate Delta1 matrix: Interaction-aware integration of miRNA and mRNA
delta1 = delta_integration.int_Delta1(miRNA, mRNA, interaction1)
print("delta1_done")

# Generate Delta2 matrix: Interaction-aware integration of methylation and mRNA
delta2 = delta_integration.int_Delta2(methyl, mRNA, interaction2)
print("delta2_done")

# Final integrated feature matrix X by averaging delta1 and delta2
X = (delta1 + delta2) / 2

# Class labels for each sample (e.g., cancer subtype)
y = pd.read_csv(r'dependent_variables.csv', sep=",", header=0, index_col=0)

# --- Model Training and Evaluation ---
# Train ANN multiple times and compute average prediction accuracy
prediction_accuracy = GAIN_BRCA_ANN.ANN(X, y, n_loops)

# Output performance
print(prediction_accuracy)

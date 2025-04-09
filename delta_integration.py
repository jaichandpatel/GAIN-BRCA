"""
# ================================================================
# GAIN-BRCA Delta Integration Module
Description:
    This module provides functions to integrate genomic data by combining gene expression 
    with regulatory interactions. Specifically, it computes delta values by integrating 
    miRNA expression with mRNA expression (int_Delta1) and methylation data with mRNA 
    expression (int_Delta2) using provided interaction data.

Functions:
    int_Delta1(miRNA, mRNA, Int1)
        Integrates miRNA and mRNA expression data based on miRNA-mRNA interaction information.
    
    int_Delta2(methyl, mRNA, Int2)
        Integrates methylation data and mRNA expression data based on CpG-mRNA interaction information.
"""
# ================================================================

# --- Import Required Libraries ---
import pandas as pd
import math
import numpy as np

# ---------------------------------------------------------------
def int_Delta1(miRNA, mRNA, Int1):
    """
    Integrates miRNA expression and mRNA expression data based on miRNA-mRNA interactions.
    
    The function iterates through each gene in the mRNA dataset. For each gene, if there is an 
    associated miRNA interaction (as given by Int1), the corresponding miRNA expression is 
    transformed and combined with the mRNA expression. If no interaction is found, the mRNA 
    expression is used directly.
    
    Parameters:
        miRNA (pd.DataFrame): DataFrame containing miRNA expression data with gene names as index.
        mRNA (pd.DataFrame): DataFrame containing mRNA expression data with gene names as index.
        Int1 (pd.DataFrame): DataFrame containing miRNA-mRNA interaction data with mRNA genes as index.
    
    Returns:
        delta1 (pd.DataFrame): DataFrame containing the integrated expression values.
    """
    # Extract list of mRNA gene names from the index
    mrna_list = mRNA.index.tolist()

    # Create empty DataFrames to hold intermediate results and the final integrated delta1 values
    df = pd.DataFrame()
    delta1 = pd.DataFrame()

    ##########################################################################################
    # Iterate over each gene in the mRNA dataset
    for i in mrna_list:
        # Extract gene expression values for the current mRNA gene
        mrna_temp = mRNA.loc[i]
    
        ######################################################################################
        # Check if the current gene has a corresponding miRNA interaction
        if i in Int1.index.tolist():
            print("Yes miRNA FOR ", i)
            ##################################################################################
            # Determine if the interaction data is a single string or a list and create a miRNA list
            if isinstance(Int1.loc[i]["Mirna"], str) == True:
                mir_list = []
                mir_list.append(Int1.loc[i]["Mirna"])
            else:
                mir_list = Int1.loc[i]["Mirna"].tolist()
    
            # For each miRNA in the interaction list, check if it exists in the miRNA dataset and integrate
            for j in mir_list:
                if j in miRNA.index.tolist():
                    mirna = miRNA.loc[j]
                    # Apply a transformation function to miRNA expression and multiply elementwise with mRNA expression
                    mirna_exp = pd.Series(mirna.apply(lambda x: 1/math.e ** (1.5*(math.log(x+1,2))))).mul(mrna_temp)
                    # Concatenate the transformed data and compute the mean across the miRNA interactions
                    df = pd.concat([df, mirna_exp], axis=1).mean(axis=1)
                else:
                    None
        else:
            print("NO miRNA FOR ", i)
            # If no interaction exists, use the original mRNA expression data
            df = mrna_temp
        # Append the result for the current gene to the delta1 DataFrame
        delta1 = pd.concat([delta1, df], axis=1)
        # (Optional debug: print(delta1.shape))
    
    # Set column names for delta1 based on the list of mRNA genes
    delta1.columns = [mrna_list]
    # Optionally, the integrated data can be saved to a file (currently commented out)
    # delta1.to_csv(r"transformed_mirna.tsv", sep="\t")
    return delta1

# ---------------------------------------------------------------
def int_Delta2(methyl, mRNA, Int2):
    """
    Integrates methylation and mRNA expression data based on CpG-mRNA interactions.
    
    The function processes each gene in the mRNA dataset. For each gene, if a CpG interaction 
    exists (from Int2), the corresponding methylation data is transformed and integrated with 
    the mRNA expression values. If no interaction is found, the mRNA expression is used directly.
    
    Parameters:
        methyl (pd.DataFrame): DataFrame containing methylation (beta values) with CpG sites as index.
        mRNA (pd.DataFrame): DataFrame containing mRNA expression data with gene names as index.
        Int2 (pd.DataFrame): DataFrame containing CpG-mRNA interaction data with mRNA genes as index.
    
    Returns:
        delta2 (pd.DataFrame): DataFrame containing the integrated methylation and expression values.
    """
    # Extract list of mRNA gene names from the index
    mrna_list = mRNA.index.tolist()

    # Create empty DataFrames to hold intermediate results and the final integrated delta2 values
    df = pd.DataFrame()
    delta2 = pd.DataFrame()

    ##########################################################################################
    # Iterate over each gene in the mRNA dataset
    for i in mrna_list:
        # Extract gene expression values for the current gene
        mrna_temp = mRNA.loc[i]
    
        ######################################################################################
        # Check if the current gene has an associated CpG interaction
        if i in Int2.index.tolist():
            print("Yes CPG FOR ", i)
            ##################################################################################
            # Determine if the interaction data is a single string or a list and create a CpG list
            if isinstance(Int2.loc[i]["CpG"], str) == True:
                CpG_list = []
                CpG_list.append(Int2.loc[i]["CpG"])
            else:
                CpG_list = Int2.loc[i]["CpG"].tolist()
    
            # For each CpG site in the list, check if it exists in the methylation dataset and integrate
            for j in CpG_list:
                if j in methyl.index.tolist():
                    meth = methyl.loc[j]
                    # Apply a transformation to methylation data and multiply elementwise with mRNA expression
                    meth_exp = pd.Series(meth.apply(lambda x: 1/math.e ** (1.5*(math.log(x+1,2))))).mul(mrna_temp)
                    # Concatenate the transformed data and compute the mean across the CpG interactions
                    df = pd.concat([df, meth_exp], axis=1).mean(axis=1)
                else:
                    None
        else:
            print("NO CPG FOR", i)
            # If no CpG interaction exists, use the original mRNA expression data
            df = mrna_temp
        # Append the result for the current gene to the delta2 DataFrame
        delta2 = pd.concat([delta2, df], axis=1)
        # (Optional debug: print(i), print(delta2.shape))
    
    # Set column names for delta2 based on the list of mRNA genes
    delta2.columns = [mrna_list]
    return delta2

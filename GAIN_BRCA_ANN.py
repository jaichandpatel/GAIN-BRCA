"""
# ================================================================
Module: GAIN_BRCA_ANN
GAIN-BRCA Artificial Neural Network Classifier
Description:
    This module implements an Artificial Neural Network (ANN) for predicting outcomes using 
    integrated genomic data. The ANN function preprocesses data, trains a Keras model using 
    stratified k-fold cross-validation, evaluates performance, and saves prediction probabilities.
    
Functions:
    ANN(X, y, n_loops)
        Trains an ANN model with cross-validation and returns the maximum achieved accuracy.
# ================================================================
"""

import tensorflow as tf
from tensorflow import keras
import pandas as pd
from sklearn.model_selection import StratifiedKFold
import numpy as np
import statistics
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, classification_report
from sklearn.metrics import cohen_kappa_score, roc_auc_score, confusion_matrix

def ANN(X, y, n_loops):
    """
    Trains an Artificial Neural Network (ANN) using stratified k-fold cross-validation.
    
    This function scales the feature data using a MinMax scaler, then repeatedly splits 
    the data into training and testing sets (using 10-fold stratified splitting) and trains a 
    sequential neural network model. It tracks and saves the predictions from the fold that 
    achieves the highest accuracy. Finally, it prints out a classification report along with 
    precision, recall, and F1 scores, and saves the predicted probabilities to a CSV file.
    
    Parameters:
        X (pd.DataFrame): Feature matrix with integrated genomic data.
        y (pd.DataFrame or pd.Series): Target labels.
        n_loops (int): Number of iterations for performing cross-validation.
    
    Returns:
        final_acc (float): The highest accuracy achieved across all folds and iterations.
    """
    # Determine the number of features (input dimension)
    feature_dim = len(X.axes[1])

    # Data scaling: normalize features using MinMaxScaler
    row = X.index.tolist()
    header = X.columns.values.tolist()
    X = pd.DataFrame(MinMaxScaler().fit_transform(X)).set_axis(row, axis=0).set_axis(header, axis=1).to_numpy()

    # Convert target labels to a 1D numpy array
    y = y.to_numpy().ravel()
    # Determine the number of classes from unique labels
    n_classes = len(np.unique(y))

    # Lists to track accuracy across folds and iterations
    acc_list_master = []
    avg_acc_list = []
    a = []

    # Variables to track the best-performing fold
    max_accuracy = 0
    saved_true_labels = None
    saved_predicted_probs = None
    saved_y_pred = None

    # Repeat the cross-validation process for the specified number of loops
    for i in range(n_loops):
        kfold = StratifiedKFold(n_splits=10, shuffle=True)
        cvscores = []
        test_acc_list = []

        # Perform 10-fold cross-validation
        for train, test in kfold.split(X, y):
            # Define a sequential model with three hidden layers and one output layer
            model = keras.Sequential([
                keras.layers.Dense(300, input_shape=(feature_dim,), activation="relu"),
                keras.layers.Dense(100, activation="relu"),
                keras.layers.Dense(25, activation="relu"),
                keras.layers.Dense(n_classes, activation="softmax")
            ])

            # Compile the model using Adam optimizer and sparse categorical crossentropy loss
            model.compile(optimizer="adam", loss="sparse_categorical_crossentropy", metrics=["accuracy"])
            # Train the model on the training set
            model.fit(X[train], y[train], epochs=20, batch_size=10, verbose=2)

            # Evaluate the model on the testing set
            scores = model.evaluate(X[test], y[test], verbose=0)
            print("%s: %.2f%%" % (model.metrics_names[1], scores[1] * 100))
            cvscores.append(scores[1] * 100)
            print("%.2f%% (+/- %.2f%%)" % (np.mean(cvscores), np.std(cvscores)))

            # Obtain prediction probabilities and predicted class labels
            yhat_probs = pd.DataFrame(model.predict(X[test], verbose=2))
            y_pred = yhat_probs.idxmax(axis=1)

            # Compute accuracy for the current fold
            Accuracy = accuracy_score(y[test], y_pred)
            test_acc_list.append(Accuracy)
            acc_list_master.append(Accuracy)

            # Save the fold with the maximum accuracy for further evaluation
            if Accuracy > max_accuracy:
                max_accuracy = Accuracy
                saved_true_labels = y[test]
                saved_predicted_probs = yhat_probs.values
                saved_y_pred = y_pred

        # Compute and store the average accuracy for the current loop
        avg_test_acc = statistics.mean(test_acc_list)
        test_acc_list.append(avg_test_acc)
        a.append(test_acc_list)

    # Output a detailed classification report for the best performing fold
    print("\n==== Classification Report ====")
    print(classification_report(saved_true_labels, saved_y_pred, digits=4))

    # Calculate additional performance metrics (weighted)
    precision = precision_score(saved_true_labels, saved_y_pred, average='weighted')
    recall = recall_score(saved_true_labels, saved_y_pred, average='weighted')
    f1 = f1_score(saved_true_labels, saved_y_pred, average='weighted')

    print(f"Precision (weighted): {precision:.4f}")
    print(f"Recall (weighted):    {recall:.4f}")
    print(f"F1 Score (weighted):  {f1:.4f}")

    # Save true labels and prediction probabilities to a CSV file for further analysis
    saved_true_labels = np.array(saved_true_labels).reshape(-1, 1)
    saved_predicted_probs = np.array(saved_predicted_probs)

    roc_df = pd.DataFrame(np.hstack((saved_true_labels, saved_predicted_probs)))
    col_names = ['TrueLabel'] + [f'Prob_Class_{i}' for i in range(saved_predicted_probs.shape[1])]
    roc_df.columns = col_names

    roc_df.to_csv("predicted_probs_selected_fold.csv", index=False)

    # Return the highest accuracy observed
    final_acc = max(acc_list_master)
    return final_acc

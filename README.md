This repository contains code and data for predicting coupling constants using machine learning models. The project explores different regression and classification techniques to evaluate their performance in capturing the relationship between molecular geometry and coupling constants.

Models Implemented

Logistic Regression
Random Forest Classifier
XGBoost Classifier
Ridge regression including Kernel Ridge Regression (KRR)
Among the models tested, XGBoost and KRR achieved the highest accuracy in predicting coupling constants, demonstrating its effectiveness in handling complex, non-linear relationships in the dataset.

Repository Structure

clean.py – Preprocessing and cleaning of the dataset
generate_batch_final.py – Script for generating training batches
train_coupling_ml.py – Training script for different ML models
coup-distance-plot.py – Visualization of coupling constants vs. distances
*.csv – Processed datasets used for training and testing
geometryfile/ and Z-0-geometry/ – Input molecular geometry data
Results

XGBoost outperformed logistic regression, lasso regression, and random forest in terms of predictive accuracy.
Visualization plots (distance_vs_coupling.png, pred_vs_actual.png) show strong correlation between predicted and actual coupling constants for XGBoost.
Model performance summary is provided in model_performance_summary.csv.

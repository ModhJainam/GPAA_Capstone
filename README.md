# GPAA Capstone Project
## Designing Explainable Machine Learning Methods for RNA-Sequencing Analysis of Atherosclerosis
Jainam Modh (jhm3ab)  
Advisors: Andrew Warren, Chunghong Mao  
UVA Biocomplexity Institute
### Summary
Atherosclerosis and coronary artery disease (CAD) account for at least 30% of all deaths globally. As such, there is a tremendous need to develop personalized medicine to prevent and reverse the development of atherosclerosis in young adults. The primary goal of this project is to identify and understand the relationship between different feature selection and machine learning methods which differentiate between early and late stages of atherosclerosis and the biological features for which they select. The selected features will expand upon the list of known genes and proteins that represent the full spectrum of molecular events underlying atherosclerosis pathogenesis
### Methods
- Data Collection and Tissue-Level RNA-Sequencing
	- Subjects: Coroner's autopsies of young adults (n=128)
	- Tissue samples: abdominal aorta (n=243) and left anterior descending artery (n=89)
	- Pathology scoring: normal, fatty streak, fibrous plaque, complex fibrous plaque
- Data Preprocessing
	- Assign binary labels: normal (early) and diseased (late)
	- Compute and standardize FPKM values and gene counts
- Feature Selection: 5-fold cross validation using 80/20 stratified train/test splits.
	- Regularized Linear Regression (LASSO)
	- Random Forest
	- Recursive Feature Elimination (RFE)
	- Differential Expression Analysis
- Training Shallow learners and XGBoost model: Use same dataset splits as the feature selection methods
	- Shallow Learners: Logistic Regression, Random Forest, SVM, Decision Tree, Naive Bayes
- Enrichment Analysis
	- Gene Sets: GO Biological Processes
- UMAP Dimensionality Reduction and Supervised Clustering
	- Correlate model results with pathology scores


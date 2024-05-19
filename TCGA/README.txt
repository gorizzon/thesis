The Cancer Genome Atlas (TCGA) is a cancer genomics program that characterized over 20,000 primary cancer samples and matched normal samples. The data used in this analysis is from a single cohort. There are 347 tumor samples and RNA gene expression data for 645 genes, DNA methylation data for 574 probes, and reverse phase protein array data for 171 proteins. The cancer samples comprise four cancer subtypes. The analysis objective is to identify these type classes from the three omics data sets.

The raw omics data was obtained from another study: https://github.com/ttriche/bayesCC/tree/master
The raw class data was obtained from the original study: https://www.nature.com/articles/nature11412.
Details on how to obtain the files are in "R scripts/raw_data/README.txt".

This study received ethical approval from the Ethics Review Board of the Faculty of Social & Behavioural Sciences at Utrecht University. FETC approval 24-0804.

Instructions:
- Run the script "data_tocsv.R" in the R scripts folder to pre-process the raw omics data
- Run the script "sample_class.R" in the R scripts folder to create a csv with the tumor class data and define the test set
- Run the script "oversampling.R" in the R scripts folder to balance the classes in the training data using SMOTE, and then write the csv files used for further analysis
- Run the script "ran_par_search.ipynb" (Jupyter Notebook) in the MoGCN scripts folder to perform random search for MoGCN hyperparameter optimization
- Run the script "run MoGCN TCGA.ipynb" (Jupyter Notebook) in the MoGCN scripts folder to run MoGCN on the data
- This will produce the "results" folder in the MoGCN scripts folder
- Run the script "glmnet.R" in the R scripts folder to perform Ridge Regression and calculate its accuracy
- Copy the file "GCN_predicted_data.csv" from the "MoGCN scripts/results" folder to the "R script/output" folder
- Run the script "gcn_accu.R" in the R scripts folder to calculate MoGCN accuracy
- Run the script "EDA.R" in the R scripts folder to obtain the plots in the folder "images"

In case of any issues please reach out to g.orizzonte@uu.nl.

Ownership of the data is with The Cancer Genome Atlas. Ownership of the MoGCN code is with its developing team (https://github.com/Lifoof/MoGCN/tree/master). The R scripts and the randomized hyperparameter search code can be reproduced according to the respository License.
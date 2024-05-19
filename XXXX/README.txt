This study, referred to as "XXXX," is currently not authorized for public dissemination. In order to uphold the confidentiality of sensitive information, all identifying details pertaining to this study have been removed. This includes any references to specific organizations, research topics, or data contents, which have been anonymized. This means any files containing potentially sensitive information are not part of this repository.

In order to access the data it must be necessary to obtain approval from UMC Utrecht and the organization carrying out the study. The data can only be shared through the secure server of UMC Utrecht or the organization, and must be analyzed from an approved laptop or server. In principle, only researchers contributing to the study can be approved for access. For further information about the acess permission please reach out to Dr. Said el Bouhaddani, S.elBouhaddani-2@umcutrecht.nl.

The omics data comprises three platforms with data from 78 patients. Platform 1 has 368 variables, Platform 2 has 57, and Platform 3 has 70. The outcome data contains two classes of patients. There are exactly 39 patients in each class.

This study received ethical approval from the Ethics Review Board of the Faculty of Social & Behavioural Sciences at Utrecht University. FETC approval 24-0804.

Instructions:
- Run the script "data_tocsv.R" in the R scripts folder to pre-process the raw omics data, create the training data set, and define the test set
- Run the script "ran_par_search.ipynb" (Jupyter Notebook) in the MoGCN scripts folder to perform random search for MoGCN hyperparameter optimization
- Run the script "run MoGCN TCGA.ipynb" (Jupyter Notebook) in the MoGCN scripts folder to run MoGCN on the data
- This will produce the "results" folder in the MoGCN scripts folder
- Run the script "glmnet.R" in the R scripts folder to perform Ridge Regression
- This will produce the file "result_ridge.csv" in the "R scripts/results" folder
- Copy the file "GCN_predicted_data.csv" from the "MoGCN scripts/results" folder to the "R script/results" folder
- Run the script "results.R" in the R scripts folder to calculate accuracy and ROC AUC of both models

In case of any issues please reach out to g.orizzonte@uu.nl.

Ownership of the MoGCN code is with its developing team (https://github.com/Lifoof/MoGCN/tree/master). The R scripts and the randomized hyperparameter search code can be reproduced according to the respository License.
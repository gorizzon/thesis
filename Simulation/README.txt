This simulation study compares the performance of MoGCN and Ridge Regression under 15 different conditions. Each condition represents once combination of number of observations, N in {100, 250, 500, 1000, 5000}, and number of features, F1 = F2 = F3 in {100, 500, 1000}. The excel table "Simulation scenarios" in the folder R scripts is a reference guide to the characteristics of each condition.

I created three sets of omics-like data using PO2PLS, a factor model designed specifically for multi-omics integration. The structural equations produce three separate omics data sets using one underlying factor structure. Additionally, they are used to generate two patient classes. 

Aside from the two parameters above which define the different conditions, more parameters were held constant to generate the data. The noise of the generated data was set to 0.7 and number of total components used to generate the data was set to F/10. Half of the total number of components was set to be shared components and the remaining half was set to be orthogonal components. In PO2PLS the shared components capture shared variance across omics platforms, while the orthogonal components capture shared variance within a single platform. Additionally, Gaussian noise was added to the generated outcome variable, which was subsequently dichotomized to replicate the categorical outcomes typically found in clinical research.

The simulated data and results folders are too large to be pushed to GitHub. The zip files can be downloaded from https://filesender.surf.nl/?s=download&token=e25a7452-1100-45cb-85d3-f78c71b3be86 . The link will remain active until 02/06/2024. To access the files after this date, please email g.orizzonte@uu.nl.

The data folder has the following structure:

data/
├── data_A/
│   ├── oma_perm01.csv
│   ├── oma_perm02.csv
│   ├── ...
│   ├── oma_perm15.csv
│   ├── omb_perm01.csv
│   ├── omb_perm02.csv
│   ├── ...
│   ├── omb_perm15.csv
│   ├── omc_perm01.csv
│   ├── omc_perm02.csv
│   ├── ...
│   ├── omc_perm15.csv
│   ├── test_perm01.csv
│   ├── test_perm02.csv
│   ├── ...
│   ├── test_perm15.csv
│   ├── class_perm01.csv
│   ├── class_perm02.csv
│   ├── ...
│   ├── class_perm15.csv
├── data_B/
│   ├── oma_perm01.csv
│   ├── oma_perm02.csv
│   ├── ...
│   ├── oma_perm15.csv
│   ├── omb_perm01.csv
│   ├── omb_perm02.csv
│   ├── ...
│   ├── omb_perm15.csv
│   ├── omc_perm01.csv
│   ├── omc_perm02.csv
│   ├── ...
│   ├── omc_perm15.csv
│   ├── test_perm01.csv
│   ├── test_perm02.csv
│   ├── ...
│   ├── test_perm15.csv
│   ├── class_perm01.csv
│   ├── class_perm02.csv
│   ├── ...
│   ├── class_perm15.csv
├── ...
├── data_T/
│   ├── oma_perm01.csv
│   ├── oma_perm02.csv
│   ├── ...
│   ├── oma_perm15.csv
│   ├── omb_perm01.csv
│   ├── omb_perm02.csv
│   ├── ...
│   ├── omb_perm15.csv
│   ├── omc_perm01.csv
│   ├── omc_perm02.csv
│   ├── ...
│   ├── omc_perm15.csv
│   ├── test_perm01.csv
│   ├── test_perm02.csv
│   ├── ...
│   ├── test_perm15.csv
│   ├── class_perm01.csv
│   ├── class_perm02.csv
│   ├── ...
│   ├── class_perm15.csv


The results folder has the following structure (file names formatted for readability):

results/
├── data_A/
│   ├── result_01/
│   │   ├── Autoencoder training loss.png
│   │   ├── Autoencoder latent data.csv
│   │   ├── top omics 1.csv
│   │   ├── top omics 2.csv
│   │   ├── top omics 3.csv
│   │   ├── SNF fused matrix.csv
│   │   ├── SNF cluster map.png
│   │   ├── GCN predicted data.csv
│   ├── result_02/
│   │   ├── Autoencoder training loss.png
│   │   ├── Autoencoder latent data.csv
│   │   ├── top omics 1.csv
│   │   ├── top omics 2.csv
│   │   ├── top omics 3.csv
│   │   ├── SNF fused matrix.csv
│   │   ├── SNF cluster map.png
│   │   ├── GCN predicted data.csv
│   ├── ...
│   ├── result_15/
│   │   ├── Autoencoder training loss.png
│   │   ├── Autoencoder latent data.csv
│   │   ├── top omics 1.csv
│   │   ├── top omics 2.csv
│   │   ├── top omics 3.csv
│   │   ├── SNF fused matrix.csv
│   │   ├── SNF cluster map.png
│   │   ├── GCN predicted data.csv
├── data_B/
│   ├── result_01/
│   │   ├── Autoencoder training loss.png
│   │   ├── Autoencoder latent data.csv
│   │   ├── top omics 1.csv
│   │   ├── top omics 2.csv
│   │   ├── top omics 3.csv
│   │   ├── SNF fused matrix.csv
│   │   ├── SNF cluster map.png
│   │   ├── GCN predicted data.csv
│   ├── result_02/
│   │   ├── Autoencoder training loss.png
│   │   ├── Autoencoder latent data.csv
│   │   ├── top omics 1.csv
│   │   ├── top omics 2.csv
│   │   ├── top omics 3.csv
│   │   ├── SNF fused matrix.csv
│   │   ├── SNF cluster map.png
│   │   ├── GCN predicted data.csv
│   ├── ...
│   ├── result_15/
│   │   ├── Autoencoder training loss.png
│   │   ├── Autoencoder latent data.csv
│   │   ├── top omics 1.csv
│   │   ├── top omics 2.csv
│   │   ├── top omics 3.csv
│   │   ├── SNF fused matrix.csv
│   │   ├── SNF cluster map.png
│   │   ├── GCN predicted data.csv
├── ...
├── data_T/
│   ├── result_01/
│   │   ├── Autoencoder training loss.png
│   │   ├── Autoencoder latent data.csv
│   │   ├── top omics 1.csv
│   │   ├── top omics 2.csv
│   │   ├── top omics 3.csv
│   │   ├── SNF fused matrix.csv
│   │   ├── SNF cluster map.png
│   │   ├── GCN predicted data.csv
│   ├── result_02/
│   │   ├── Autoencoder training loss.png
│   │   ├── Autoencoder latent data.csv
│   │   ├── top omics 1.csv
│   │   ├── top omics 2.csv
│   │   ├── top omics 3.csv
│   │   ├── SNF fused matrix.csv
│   │   ├── SNF cluster map.png
│   │   ├── GCN predicted data.csv
│   ├── ...
│   ├── result_15/
│   │   ├── Autoencoder training loss.png
│   │   ├── Autoencoder latent data.csv
│   │   ├── top omics 1.csv
│   │   ├── top omics 2.csv
│   │   ├── top omics 3.csv
│   │   ├── SNF fused matrix.csv
│   │   ├── SNF cluster map.png
│   │   ├── GCN predicted data.csv


This study received ethical approval from the Ethics Review Board of the Faculty of Social & Behavioural Sciences at Utrecht University. FETC approval 24-0803.

Instructions:
- Run the script "data_gen_code.R" in the R scripts folder to generate the synthetic data
- This will produce the "data" folder available at the link above
- Run the script "run mogcn.ipynb" (Jupyter Notebook) in the MoGCN scripts folder to run MoGCN on one iteration of the simulation (with 15 conditions)
- The script needs to run 20 times for each of the 20 iterations
- This will produce the "results" folder available at the link above
- Run the script "glmnet.R" in the R scripts folder to perform Ridge Regression
- The MoGCN results folder needs to be in the R scripts folder, like "R script/results/...", OR the path file needs to be specified in the code
- Run the script "AUC.R" in the R scripts folder to calculate accuracy and ROC AUC of the 20 * 15 MoGCN and Ridge predictions
- Run the script "results_writing.R" in the R scripts folder to summarize the 20 * 15 accuracy and AUC scores into long-format tables and calculate their averages per condition across iterations
- The script above also produces the plot in the images folder

In case of any issues please reach out to g.orizzonte@uu.nl.

Ownership of the MoGCN code is with its developing team (https://github.com/Lifoof/MoGCN/tree/master). The synthetic data and R scripts can be reproduced according to the respository License.

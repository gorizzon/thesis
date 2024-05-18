These are the processed data files used for analysis.

- All features were rescaled to 0-1 by column.
- The files ending in "train_bal"contain only the training observations after balancing.
- The files ending only in "_bal" contain all observations after balancing.
- The files starting with "Meth", "RPPA", "mRNA" contain the features of the three omics platforms, respectively: DNA Methylation, reverse phase protein array, and RNA gene expression.
- The files starting with "outcome" contain the classes of the tumor samples.
- The file "test_sample.csv" contains the IDs of the observations in the test set.
# Equivalent Imputation Models for Handling Missing Data in Compositional Geochemical Databases of Geothermal Fluids

The data and code corresponds to the paper with the same name as this repository, were we are analyzing compositional geothermal fluid data from different geothermal boreholes/wells. The purpose is to impute missing data in the whole dataset so it can be used for multiple purposes such as predicting the bottom-hole temperature of a geothermal resource.

The content of the repository is mainly:

    1. Data (recollected and used for experimentation)
    2. Notebooks (developed in Python and R)

The content in the Data folder is the following: 

    - WCGDb.csv (A comma-separated values file containing the Working Compositional Geothermal Database)
    - simplifiedWCGDb.csv (a simplified version of the WCGDb)

The content of the notebooks is the following:

    - EquivImpGFD_2021.ipynb (the main notebook)
    - ImputatorTester_alg_541a.ipynb (imputations using single imputation algorithms WITHOUT parameters tunning)
    - ImputatorTester_alg_541b.ipynb (imputations using single imputation algorithms WITH parameters tunning)
    - ImputatorTester_alg_521.ipynb (imputations using the MICE)
    - ImputatorTester_alg_522.ipynb (imputations using the multiple imputation algorithms WITH parameters tunning)
    - ImputationQualityAssessment.ipynb (imputations evaluation)
    - ImputationQualityAssessment1.md (fancier notebook developed in R for showing imputations evaluation)



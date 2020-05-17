# LURK-Vecchia

Simultaneous variable seletion and estimation of LUR models with spatiotemporally correlated errors that is scalable for big data

# Citation
coming soon ...

# Installation 
Based on the package GPvecchia (https://cran.r-project.org/web/packages/GPvecchia/index.html) and
ncvreg (https://cran.r-project.org/web/packages/ncvreg/index.html). 
Install the following R packages: GPvecchia,ncvreg, fields, Matrix

Download LURK_Functions.R to your local machine and use R to source:
```
source("/Location/of/your/functions/LURK_Functions.R") 
``` 

# Data

## NO<sub>2</sub>
The NO<sub>2</sub> analysis data are stored in a single csv file in the [data](https://github.com/NIEHS/LURK-Vecchia/tree/master/data) subfolder. 

```
US_ST_NO2_Data.csv
```

Columns 1 - 6 are latitude, longitude, date, NO2 concentration, and time in fractional years, respectively. Columns 7 - 145 make up the design matrix of covariates. The covariate names are followed by the distance hyperparameter (i.e. buffer, decay range).
To see all of the field names, call: 
```
colnames(US_ST_NO2_Data)
```

## Simulations

The simulation analysis data are stored in Excel files in the [data](https://github.com/NIEHS/LURK-Vecchia/tree/master/data) subfolder. 
There is an Excel file for each spatiotemporal parameter that is varied: nugget-to-sill ratio, total variance (sill), spatial range, and temporal range.
For each spatiotemporal parameter there is a training set and test set, stored in separate files. The training and test sets are based on the same
simulation.

```
Vecchia_ST_Simulation_NUG2SILL5_20191107.xlsx
Vecchia_ST_Simulation_NUG2SILL5_test_20191107.xlsx
Vecchia_ST_Simulation_VARIANCE5_20191107.xlsx
Vecchia_ST_Simulation_VARIANCE5_test_20191107.xlsx
Vecchia_ST_Simulation_SPATIAL_RANGES5_20191107.xlsx
Vecchia_ST_Simulation_SPATIAL_RANGES5_test_20191107.xlsx
Vecchia_ST_Simulation_TEMPORAL_RANGES5_20191107.xlsx
Vecchia_ST_Simulation_TEMPORAL_RANGES5_test_20191107.xlsx
```
Each Excel file contains a sheet for each simulation scenario - that is there are 20 sheets for 20 variations of the spatiotemporal parameter. 
Each sheet is simply labeled "1", "2",..., etc. for the simulation scenario. 
The first sheet also contains data for the true observations without measurement error and the covariates. 
Hence, in each simulation you see differences for the first sheet. 

As an example, we read in the first sheet for the nugget-to-sill simulations: 
```
data1 <- as.data.frame(read_excel("Vecchia_ST_Simulation_NUG2SILL5_20191107.xlsx", 
                                  sheet = "1", col_names = FALSE))
```

# Running the Simulation Analysis 

# Running the NO<sub>2</sub> Analysis 

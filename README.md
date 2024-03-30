**Time Series Analysis Code**
This repository contains R code for conducting time series analysis on chronic lifestyle diseases data. The analysis includes various tests, such as MAN-KENDALL and COX-STUART tests, calculation of autocorrelation and partial autocorrelation, periodogram analysis, and seasonal decomposition using LOESS and TBATS methods.

**Instructions:**
Installation of Required Packages:

**The following R packages need to be installed before running the code:**
trend
TSA
forecast
**Loading Required Libraries:**

**The code loads necessary libraries including trend, TSA, and forecast.**
Importing Dataset:

**The code prompts to select the dataset file (CSV format) which contains the time series data.**
Handling Missing Values:

**Missing values in the dataset are handled using the na.omit() function.**
Conducting Tests:

**MAN-KENDALL Test**:
Tests for trend detection in the data.

**COX-STUART Test**:
Tests for the presence of cyclical patterns in the data.

**Plotting Autocorrelation and Partial Autocorrelation:**
Autocorrelation and Partial Autocorrelation plots are generated for each column in the dataset.

**Calculating Periodicity:**
Periodogram analysis is performed to identify the dominant frequencies and corresponding periods in the data.

**Seasonal Decomposition:**
Seasonal decomposition is conducted using both LOESS and TBATS methods to identify seasonal patterns in the data.

**Usage:**
Clone or download the repository.
Ensure R and RStudio are installed on your system.
Install required packages using install.packages().
Run the provided R script in RStudio.
Follow the instructions in the R script to select the dataset file.
Analyze the results generated by the code.

**Note:**
Make sure the dataset is in CSV format and follows the required structure.
Ensure the dataset does not contain any missing or erroneous values for accurate analysis.

**Citation:**
If you use this code or the analysis results in your research, please cite the original paper:

**Patel JC, Khurana P, Sharma YK, Kumar B, Ragumani S (2018) Chronic lifestyle diseases display seasonal sensitive comorbid trend in human population evidence from Google Trends. PLOS ONE 13(12): e0207359. https://doi.org/10.1371/journal.pone.0207359**

**Patel JC, Khurana P, Sharma YK, Kumar B, Sugadev R. Google trend analysis of climatic zone based Indian severe seasonal sensitive population. BMC Public Health. 2020 Mar 12;20(1):306. doi: 10.1186/s12889-020-8399-0. PMID: 32164654; PMCID: PMC7069044.**

**Patel JC, Gupta A, Kumar P, Waidha KM, Deep A, Kumar A, Katare DP, Sharma AK. Cardiovascular diseases display etiological and seasonal trend in human population: Evidence from seasonal cardiovascular comorbid diseases (SCCD) index. Am J Hum Biol. 2023 Jun;35(6):e23867. doi: 10.1002/ajhb.23867. Epub 2023 Jan 18. PMID: 36651684.**

**License:**
This code is provided under the MIT License. Feel free to use, modify, or distribute it according to your needs.

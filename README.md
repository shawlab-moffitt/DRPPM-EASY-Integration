# DRPPM - EASY - Integration

# Introduction



# Installation

* Download ZIP file from https://github.com/shawlab-moffitt/DRPPM-EASY-Integration
* Unzip and load into directory as a project in R Studio
* Open the ‘App.R’ script and write in user input files and options as directed at the top of the script
  * ‘App.R’ script begins with example files available on the front page and within the ExampleData Folder
* Press ‘Run App’ button in R Studio to run in application or browser window and enjoy!
  * The app script will install any missing packages that the user may not have locally


# Requirments

* `R` - https://cran.r-project.org/src/base/R-4/
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

# R Dependencies

|  |  |  |  |  |
| --- | --- | --- | --- | --- |
| shiny_1.6.0 | shinythemes_1.2.0 | shinyjqui_0.4.0 | shinycssloaders_1.0.0 | stringi_1.7.6 |
| dplyr_1.0.7 | tidyr_1.1.3 | readr_2.0.1 | stringr_1.4.0 | DT_0.18 |
| ggplot2_3.3.5 | plotly_4.9.4.1 | enrichplot_1.12.2 | ggVennDiagram_1.2.0 | ggrepel_0.9.1 |
| rhoR_1.3.0.3 | limma_3.48.3 | clusterProfiler_4.0.5 | limma_3.48.3 | GSVA_1.40.1 |
| BiocManager_1.30.16 | reshape2_1.4.4 | ggpubr_0.4.0 |  |  |

# Required Files

* **MSigDB Gene Set Names:**
  * These [gene set files](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration/tree/main/GeneSets) were gathered from the [Molecular Signatures Database (MSigDB)](http://www.gsea-msigdb.org/gsea/msigdb/index.jsp) as separate collections and processed through R to generate a master gene set file with catagorical labels to use for GSEA and ssGSEA analysis.
  * This is used mainly for the UI for gmt category selection.
* **MSigDB Gene Set RData List:**
  * The RData gene set list is a more refined format of the gene set table.
  * This is a named list with over 32,000 gene sets from MSigDB paired with the genes they consist of.
  * This list is used for the back end analysis.

# App Features

## Scatter Plot Comparison

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration/blob/main/App_Demo_Pictures/EASY_INT_Scatter.png?raw=true)

1. The user may upload teo files gathered from [ssGSEA analysis](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY#singer-sample-gsea-box-plot) in the main DRPPM-EASY app to compare 
   * In theory, the plot gathers values from the third column in both of the files uploaded as long as the first two columns are the same. So as long as the data follows this format a plot should be generated.
     * More adjustments will be made in future edits to ensure compatibility and freedom to use other data sets
2. The user may choose to log transform either axis
3. The axis label names may also be adjusted
4. The table shown below the figure may be downloaded for future use and analysis and even used in the following tab for the [correlation rank plot]()



## Ranked Feature Correlation



## Matrix Comparison with Reciprical GSEA



# Future Enhancments



# Questions and Comments

Please email Alyssa Obermayer at alyssa.obermayer@moffitt.org if you have any further comments or questions in regards to the R Shiny Application.


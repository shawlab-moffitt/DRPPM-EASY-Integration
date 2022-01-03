# DRPPM - EASY - Integration

# Introduction

This is an extention of the [DRPPM Expression Analysis ShinY (EASY) App](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY) for further integration of result files and matrices comparison with differential gene expression and gene set enrichment analysis. Here we allow for the comparison of data sets such as from single sample GSEA and expression matrices. The comparison of expression data is one of the main features we elaborate on in this app. The user may upload two expression matrices, in our example we uploaded mRNA and Protein expression data, with their corresponding meta data. Through the functionality of the DRPPM-EASY-Integration app we are abloe to compare log fold change of the expression and generate gene sets of the significantly up regulated and down regulated genes. These gene set were further used for GSEA and ssGSEA on their opposing matrix to illustrate any similarity in gene regulation. Additionally, these gene set are further compared with GMT files of published gene sets from the Molecular Signatures Database (MSigDB) with corresponding statiistal calculations to rank their similarity, such as Fishers Exact Test, Cohen's Kappa, and the Jaccard Index. Below you may see a flow chart of the DRPPM-EASY pipeline, where this Integration app represents segment C.

<img src="https://github.com/shawlab-moffitt/DRPPM-EASY-Integration/blob/main/App_Demo_Pictures/EASY_FlowChart_INT.png" width="900">


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

1. The user may upload two files gathered from [ssGSEA analysis](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY#singer-sample-gsea-box-plot) in the main DRPPM-EASY app to compare 
   * In theory, the plot gathers values from the third column in both of the files uploaded as long as the first two columns are the same. So as long as the data follows this format a plot should be generated.
     * More adjustments will be made in future edits to ensure compatibility and freedom to use other data sets
2. The user may choose to log transform either axis
3. The axis label names may also be adjusted
4. The table shown below the figure may be downloaded for future use and analysis and even used in the following tab for the [correlation rank plot](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration/blob/main/README.md#ranked-feature-correlation).

## Ranked Feature Correlation

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration/blob/main/App_Demo_Pictures/EASY_INT_Corr.png?raw=true)

1. The user must upload a primary feature file which consists of 3 columns, sample name, sample type, and a value and an expression matrix.
   * Typically, and as shown in the example data, we use a ssGSEA score file and the primary feature file.
   * In future edits more compatibility will be added
2. When using the ssGSEA score file as the primary feature file the name of the gene sets within the file will show up to choose which one to correlate the expression matrix with.
   * The user may use the merge ssGSEA table from the previous [scatter plot tab](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration/blob/main/README.md#scatter-plot-comparison) as input.
3. The user may chose the correlation method which options of Spearman, Pearson, or Kendall, as well as the correlation cutoff values used in coloring the genes on the figure
4. The table of correlation values below the figure may be download with a file name of your choice.

## Matrix Comparison with Reciprical GSEA

### Matrix Upload

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration/blob/main/App_Demo_Pictures/EASE_INT_MatUpload.png?raw=true)

1. Assign names to the two matrices the user will be comparing for easier identification
2. Upload expression matrices
4. Upload meta data

### Fold Change Scatter Plot

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration/blob/main/App_Demo_Pictures/EASE_INT_logFCscatter.png?raw=true)

1. The user may choose the sample type to compare the expression fold change of
2. Select genes to annotate on the plot
3. When hovering of a point this text box will appear giving information on that gene
4. The table below the figure may be downloaded for future use

### Reciprocal GSEA

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration/blob/main/App_Demo_Pictures/EASE_INT_GSEA.png?raw=true)

1. The user may select the number of top genes to include in the gene set generation of the matrices and the P-value cutoff for the GSEA calculation
   * For example, this figure shows that after performing differential gene expression on both matrices, we are selecting the top 100 up and down regulated genes of each matrices
   * These gene sets will be used to perform a GSEA calculation on their opposing gene set with the P-value cutoff designated
2. Above the enrichment plots displays the normalizxed enrichment scores and P-value for each gene set

### Reciprocal Single Sample GSEA

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration/blob/main/App_Demo_Pictures/EASY_INT_ssGSEA.png?raw=true)

1. The user may choose the ssGSEA scoring method based off the choices of ssGSEA, GSVA, z-score, or plage

### Venn Diagram and Statistical Analyses

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration/blob/main/App_Demo_Pictures/EASY_INT_Venn.png?raw=true)

1. The user may select significance cutoffs values for the fold change and P-value when determining the differentially expressed genes in each matrix which is going to be compared to the gene sets
2. The user may select a category of gene sets based on the MSigDB categories to compare to the differntially expressed genes
3. The user may also upload their own .gmt file to analyze with their matrices.
   * Proper formatting of .gmt files is provided by the Broad Institute and can be found [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)
4. The statistics table may be downloaded for future use


# Future Enhancments

* More compatibility with a wider range of files for intergrative analysis comparison
* Mouse-to-human gene name conversion for comparison analysis


# Questions and Comments

Please email Alyssa Obermayer at alyssa.obermayer@moffitt.org if you have any further comments or questions in regards to the R Shiny Application.


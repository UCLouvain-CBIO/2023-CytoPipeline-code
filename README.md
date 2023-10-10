# 2023-CytoPipeline-code

[![DOI](https://zenodo.org/badge/664524730.svg)](https://zenodo.org/badge/latestdoi/664524730)
[![license](https://img.shields.io/badge/license-GPL3.0-blue)](https://opensource.org/licenses/GPL-3.0)

Code and raw data to reproduce the results in the CytoPipeline and CytoPipelineGUI article.

User's instruction:

## Set-up
- do a `git clone` of this [github repo](https://github.com/UCLouvain-CBIO/2023-CytoPipeline-code)

- create a R project pointing to the root directory

- if not done yet, install the needed packages: 
*CytoPipeline*, *CytoPipelineGUI*, *pipeComp*, 
and *CytoPipelineUtils*: 

```
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("CytoPipeline")
# BiocManager::install("CytoPipelineGUI")
# BiocManager::install("pipeComp")
#
# if (!require("devtools", quietly = TRUE))
#     install.packages("devtools")
# devtools::install_github(
#     "https://github.com/UCLouvain-CBIO/CytoPipelineUtils/")
```

## Generating the CytoPipeline caches 

- Execute the R script : `source("./R/1_build_and_run_pipelines.R")`  

This is likely to take several minutes. Also note that during the run 
of the first pipeline (*HBVMouse_PQC_error* experiment), 
some error messages are encountered, these are **expected** errors:  

> Error : BiocParallel errors
>
>  55 remote errors, element index: 1, 2, 3, 4, 5, 6, ...
>
>  0 unevaluated and other errors
>
>  first remote error:
>
> Error in removeDebrisFlowClustTmix(new("flowFrame", exprs = structure(c(401.926147460938, : argument "nClust" is missing, with no default

## Generate the pipeComp results

- Execute the R script : `source("./R/2_build_and_run_benchmark.R")`  

This is also likely to take several minutes. Also note that during the run, 
*pipeComp* might generate warnings as below, this is not to worry about.  

> 'as.is' should be specified by the caller; using TRUE

## Produce the plots in article section : _Results/Visual Assessments_

- Execute the lines (one by one) of the following *R* script:   
`./R/3_Results_VisualAssessment.R`

## Produce the plots in article section : _Results/Benchmark results_

- Execute the lines (one by one) of the following *R* script:   
`./R/4_Results_VisualAssessment.R`



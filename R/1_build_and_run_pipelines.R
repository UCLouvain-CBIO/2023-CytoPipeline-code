require("CytoPipeline")
require("CytoPipelineUtils")

# preliminaries: path to raw data files, json inputs, and outputs
rawDataDir <- "./rawData"
jsonDir <- "./json"
#outputDir <- tempDir()
outputDir <- "."

# build CytoPipeline objects
###############################################################################
# here we define 4 pre-processing pipelines
# pipeline 1: with PeacoQC as QC step, but one of the steps has an error in it, 
# i.e. the 'remove_debris' is defined with 'nClusters' as a parameter 
# instead of 'nClust'
# pipeline 2: with PeacoQC as QC step, 'nClust' set to 3 in the 'remove_debris'
# step
# pipeline 3: with PeacoQC as QC step, 'nClust' in 'remove_debris' set to 2
# pipeline 4: with flowAI as QC step

expNames <- c("HBVMouse_PQC_error",
              "HBVMouse_PQC_nC3",
              "HBVMouse_PQC",
              "HBVMouse_flowAI")

jsonFiles <- file.path(jsonDir,
                       paste0(expNames, ".json"))

(sampleFiles <- list.files(path = rawDataDir,
                           pattern = ".fcs",
                           full.names = TRUE))

# construct phenoData data.frame
day <- substr(basename(sampleFiles), 1, 3)
well <- substr(basename(sampleFiles), 5, 7)

pData <- data.frame(day = day, well = well, 
                    row.names = basename(sampleFiles))

pipelineList <- list()

p <- 3

for (p in seq_along(expNames)) {
    pipelineList[[p]] <- 
        CytoPipeline(jsonFiles[p],
                     experimentName = expNames[p],
                     sampleFiles = sampleFiles,
                     pData = pData)
}

# run CytoPipeline objects
###############################################################################

for (p in seq_along(expNames)) {
    try(execute(pipelineList[[p]], 
                path = outputDir,
                rmCache = TRUE))
}
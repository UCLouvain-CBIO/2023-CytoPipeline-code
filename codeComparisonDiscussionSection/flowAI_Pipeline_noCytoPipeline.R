library(flowCore)
library(CytoPipeline)
library(CytoPipelineUtils)

rawDataDir <- "/rawData"           # input raw data directory here
compensationMatrixPath <- 
    "/compensation/compMatrix.csv" #input path to compensation matrix here 
resultsDir <- "."                  # input target results diretory here


# locates the 55 raw fcs files
sampleFiles <- list.files(path = rawDataDir, pattern = "*.fcs") 

res <- flowCore::read.flowSet(sampleFiles,
                              ...)

# reading the sample files from disk
fs <- flowCore::read.flowSet(sampleFiles = sampleFiles,
                             truncate_max_range = FALSE,
                             min.limit = NULL)

# reading compensation matrix
compensationMatrix <- read.csv(compensationMatrixPath,
                               check.names = FALSE,
                               row.names = 1)

## SCALE TRANSFORM steps

# remove margins (outliers) from all samples
channelSpecs <- list()
channelSpecs[["AllFluoChannels"]] <- c(-300, 262144)
fs <- CytoPipeline::removeMarginsPeacoQC(
    fs,
    channelSpecifications = channelSpecs
)

# randomly select a few samples to calculate scale transformation
# and create an aggregated sample of 10,000 events
nRandomSamples <- 4
randomSampleFiles <- NULL
aggSampleSize <- 10000
with_seed(seed = 0,
          {
            whichSamples <- sample(seq_along(sampleFiles), nRandomSamples)
            randomSampleFiles <- sampleFiles[whichSamples]
          })

aggFF <- CytoPipeline::aggregateAndSample(
    fs = fs[randomSampleFiles], 
    nTotalEvents = aggSampleSize,
    seed = 0)


# run compensation for the aggregated file
aggFF <- flowCore::compensate(aggFF, compensationMatrix)              

# estimate scale transformations
scaleTransfos <- CytoPipeline::estimateScaleTransforms(
    aggFF,
    fluoMethod = "estimateLogicle",
    scatterMethod = "linearQuantile",
    scatterRefMarker = "CD4"
)

## pre-processing steps (file by file)

for (i in seq_len(flowCore::length(fs))){
    ff <- fs[[i]]
    
    # QC in time using flowAI
    
    channel2Exclude <-
        flowCore::colnames(ff)[!CytoPipeline::areSignalCols(ff)]
    ff <- flowAI::flow_auto_qc(
        fcsfiles = ff,
        ChExcludeFS = channel2Exclude,
        chExcludeFM = channel2Exclude,
        html_report = FALSE,
        mini_report = FALSE,
        folder_results = FALSE,
        remove_from = "all",
        second_fractionFR = 0.1,
        deviationFR = "MAD",
        alphaFR = 0.01,
        decompFR = TRUE,
        outlier_binsFS = FALSE,
        pen_valueFS = 500,
        max_cptFS = 3,
        sideFM = "both",
        neg_valuesFM = 1
    )
    
    # run compensation
    ff <- flowCore::compensate(ff, compensationMatrix)
    
    # transform using scale transformation
    ff <- flowCore::transform(ff, transList)
    
    # doublets removal
    ff <- CytoPipeline::removeDoubletsCytoPipeline(
        ff,
        areaChannels = "FSC-A",
        heightChannels = "FSC-H",
        nmads = 3
    )
    
    # debris removal
    ff <- CytoPipelineUtils::removeDebrisFlowClustTmix(
        ff,
        FSCChannel = "FSC-A",
        SSCChannel = "SSC-A",
        nClust = 2,
        level = 0.97,
        B =100,
        verbose = TRUE
    )
    
    # dead cells removal
    ff <- CytoPipelineUtils::removeDeadCellsDeGate(
        ff,
        LDMarker = "Live & Dead"
    )
    
    # export results after pre-processing (including scale transformations)
    CytoPipeline::writeFlowFrame(
        ff,
        dir = resultsDir,
        suffix = "°preprocessed",
        format = "fcs")
    
}


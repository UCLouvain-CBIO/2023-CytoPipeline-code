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

# run compensation for each file
fs <- flowCore::compensate(fs, compensationMatrix)

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
              

# estimate scale transformations
scaleTransfos <- CytoPipeline::estimateScaleTransforms(
    aggFF,
    fluoMethod = "estimateLogicle",
    scatterMethod = "linearQuantile",
    scatterRefMarker = "CD4"
)

## pre-processing steps (file by file)

# transform all samples using scale transformations
fs <- flowCore::transform(fs, transList)

for (i in seq_len(flowCore::length(fs))){
    ff <- fs[[i]]
    
    # QC in time using PeacoQC
    channel4QualityControl <-
        flowCore::colnames(ff)[CytoPipeline::areSignalCols(ff)]
    ff <- PeacoQC::PeacoQC(
        ff = ff,
        channels = channel4QualityControl,
        report = FALSE,
        plot = FALSE,
        save_fcs = FALSE,
        output_directory = NULL,
        min_cells = 150,
        max_bins= 500,
        step = 500,
        MAD = 6,
        IT_limit = 0.6,
        force_IT = 150,
        peak_removal = 1/3,
        min_nr_bins_peakdetection = 10)
    
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


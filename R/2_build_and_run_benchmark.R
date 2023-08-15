require(pipeComp)
require(CytoPipeline)
require(CytoPipelineUtils)

#pipeCompOutputDir <- file.path(tempdir(),"pipeComp_output")
pipeCompOutputDir <- "./pipeComp_output"
if(!dir.exists(pipeCompOutputDir)) {
    dir.create(pipeCompOutputDir)
}

rawDataDir <- "./rawData"
compensationDir <- "./compensation"

resultsDir <- "."

executePipeComp <- TRUE

# create appropriate BiocParallel back-end
logDir <- "./BPlog2"
if(!dir.exists(logDir)) {
    dir.create(logDir)
}

bp <- BiocParallel::SnowParam(workers = 6,
                              progressbar = TRUE,
                              log = TRUE,
                              logdir = logDir,
                              stop.on.error = FALSE)

runPipeComp <- function(rawDataDir, 
                        gatingResultsDir, 
                        compensationDir, 
                        outputDir,
                        backend,
                        sampleFileIndicesToRemove = integer()){
    
    computeConcordanceIndicators <- function(ff, gtFF, originalIDs){
        if (!inherits(ff, "flowFrame")) {
            stop("ff should be a flowFrame")
        }       
        if (!inherits(gtFF, "flowFrame")) {
            stop("gtFF should be a flowFrame")
        }
        if (!"Original_ID" %in% flowCore::colnames(ff)) {
            stop("ff should contain `Original_ID` channel")
        }
        if (!"Original_ID" %in% flowCore::colnames(gtFF)) {
            stop("gtFF should contain `Original_ID` channel")
        }
        
        flaggedGood <- flowCore::exprs(ff)[, "Original_ID"]
        if(!all(flaggedGood %in% originalIDs)) {
            stop("Found original ids in ff which are not ",
                 "in provided originalIDs vector")
        }
        flaggedBad <- setdiff(originalIDs, flaggedGood)
        gtGood <- flowCore::exprs(gtFF)[, "Original_ID"]
        if(!all(gtGood %in% originalIDs)) {
            stop("Found original ids in gtFF which are not ",
                 "in provided originalIDs vector")
        }
        gtBad <- setdiff(originalIDs, gtGood)
        
        truePositives <- intersect(flaggedBad, gtBad)
        falsePositives <- setdiff(flaggedBad, gtBad)
        trueNegatives <- intersect(flaggedGood, gtGood)
        falseNegatives <- setdiff(flaggedGood, gtGood)
        
        nP <- length(gtBad)
        nTP <- length(truePositives)
        nFP <- length(falsePositives)
        nN <- length(gtGood)
        nTN <- length(trueNegatives)
        nFN <- length(falseNegatives)
        
        sensi <- ifelse(nP > 0, nTP/nP, 1.)
        speci <- ifelse(nN > 0, nTN/nN, 1.)
        preci <- ifelse((nTP+nFP)>0, nTP/(nTP+nFP), 1.)
        recall <- ifelse((nTN+nFN)>0, nTN/(nTN+nFN), 1.)
        
        return(list(flaggedGood = flaggedGood,
                    flaggedBad = flaggedBad,
                    gtGood = gtGood,
                    gtBad = gtBad,
                    truePositives = truePositives,
                    falsePositives = falsePositives,
                    trueNegatives = trueNegatives,
                    falseNegatives = falseNegatives,
                    nP = nP, 
                    nTP = nTP, 
                    nFP = nFP, 
                    nN = nN, 
                    nTN = nTN, 
                    nFN = nFN,
                    sensi = sensi, 
                    speci = speci, 
                    preci = preci, 
                    recall = recall))
    }
    
    

    
    # recreate CytoPipeline objects from cache
    pipL_FJ <- CytoPipeline::buildCytoPipelineFromCache(
        experimentName = "HBVMouse_FJ", path = gatingResultsDir
    )
    
    sampleFiles <- file.path(rawDataDir, sampleFiles(pipL_FJ))

    pData <- CytoPipeline::pData(pipL_FJ)
    
    # names will be used by pipeComp to designate the datasets
    names(sampleFiles) <- substr(row.names(pData), 1, nchar(row.names(pData))-4)
    
    compensationMatrixPathD91 <- file.path(
        compensationDir,
        "Compensations Liver D91.csv")
    compensationMatrixPathD93 <- file.path(
        compensationDir,
        "Compensations Liver D93.csv")
    
    transFormListPath <- file.path(
        gatingResultsDir, "HBVMouse_FJ/output/RDS",
        "scaleTransformList.rds")
    
    transformList <- readRDS(file = transFormListPath)
    
    
    doNothing <- function(x, ...) {
        # does nothing
        x
    }
    
    pip <- 
        PipelineDefinition(
            functions = list(
                flowframe_read = function(x){
                    #message("in flowframe_read()")
                    library(CytoPipeline)
                    res <- list()
                    res$ff <- readSampleFiles(x, 
                                              truncate_max_range = FALSE,
                                              min.limit = NULL,
                                              pData = pData
                    )
                    res$GT <- CytoPipeline::getCytoPipelineFlowFrame(
                        pipL_FJ, path = gatingResultsDir, 
                        sampleFile = basename(x), 
                        whichQueue = "pre-processing",
                        objectName = "lymphocytes_gate_obj"
                    )
                    res$nOriginalEvents <- flowCore::nrow(res$ff)
                    #message("end flowframe_read()")
                    res
                },
                remove_margins = function(x){
                    #message("in remove_margins()")
                    library(CytoPipeline)
                    x$ff <- 
                        removeMarginsPeacoQC(x$ff,
                                             channelSpecifications = 
                                                 list(
                                                     AllFluoChannels = 
                                                         c(-300, 262144)
                                                 ))
                    #message("end remove_margins()")
                    x
                },
                compensate1 = function(x, compensate1_meth) {
                    #message("compensate1()")
                    library(CytoPipeline)
                    x$ff <- 
                        get(compensate1_meth)(
                            x$ff, 
                            matrixSource = "pData",
                            pDataVar = "day",
                            pDataPathMapping = list(
                                D91 = compensationMatrixPathD91,
                                D93 = compensationMatrixPathD93))
                    #message("end compensate1()")
                    x
                },
                qualityControl = function(x, qc_meth, preTransformForQC) {
                    #message("in qualityControl()")
                    library(CytoPipeline)
                    x$ff <-
                        get(qc_meth)(
                            x$ff,
                            preTransform = preTransformForQC,
                            transList = transformList)
                    #message("end qualityControl()")
                    x
                },
                compensate2 = function(x, compensate2_meth) {
                    #message("in compensate2()")
                    library(CytoPipeline)
                    x$ff <-
                        get(compensate2_meth)(
                            x$ff,
                            matrixSource = "pData",
                            pDataVar = "day",
                            pDataPathMapping = list(
                                D91 = compensationMatrixPathD91,
                                D93 = compensationMatrixPathD93))
                    #message("end compensate2()")
                    x
                },
                remove_doublets = function(x) {
                    #message("in remove_doublets()")
                    library(CytoPipeline)
                    x$ff <- removeDoubletsCytoPipeline(
                        x$ff,
                        areaChannels = "FSC-A",
                        heightChannels = "FSC-H",
                        nmads = 3)
                    #message("end remove_doublets()")
                    x
                },
                remove_debris = function(x) {
                    #message("in remove_debris()")
                    library(CytoPipelineUtils)
                    x$ff <- removeDebrisFlowClustTmix(
                        x$ff,
                        FSCChannel = "FSC-A",
                        SSCChannel = "SSC-A",
                        nClust = 2,
                        level = 0.97,
                        B = 100)
                    #message("end remove_debris()")
                    x
                },
                remove_dead_cells = function(x) {
                    #message("in remove_dead_cells()")
                    library(CytoPipelineUtils)
                    x$ff <- removeDeadCellsDeGate(
                        x$ff,
                        preTransform = TRUE,
                        transList = transformList,
                        LDMarker = "Live & Dead")
                    #message("end remove_dead_cells()")
                    x
                }
            ),
            descriptions = list(
                flowframe_read = "reading flowFrame from fcs file on disk",
                remove_margins = "removing margin events",
                compensate1 = "optional compensation before QC",
                qualityControl = 
                    "QC, i.e. removing chunks of events with non stable distribution in time",
                compensate2 = "optional compensation after QC",
                remove_doublets = "removing non singlet events",
                remove_debris = "removing debris events",
                remove_dead_cells = "removing dead cells events"
            ),
            evaluation = list(
                flowframe_read = function(x) {
                    c(nEvents = flowCore::nrow(x$ff))
                },
                remove_margins = function(x) {
                    c(nEvents = flowCore::nrow(x$ff))
                },
                qualityControl = function(x) {
                    c(nEvents = flowCore::nrow(x$ff))
                },
                remove_doublets = function(x) {
                    c(nEvents = flowCore::nrow(x$ff))
                },
                remove_debris = function(x) {
                    c(nEvents = flowCore::nrow(x$ff))
                },
                remove_dead_cells = function(x) {
                    concordanceIndicators <-
                        computeConcordanceIndicators(
                            ff = x$ff,
                            gtFF = x$GT,
                            originalIDs = 1:(x$nOriginalEvents)
                        ) 
                    
                    c(nEvents = concordanceIndicators$nTN +
                          concordanceIndicators$nFN,
                      nGTEvents = concordanceIndicators$nN,
                      sensitivity = concordanceIndicators$sensi,
                      specificity = concordanceIndicators$speci,
                      precision = concordanceIndicators$preci,
                      recall = concordanceIndicators$recall)
                }
            )
        )
    
    pip
    
    alternatives <- list(
        compensate1_meth=c("compensateFromMatrix", "doNothing"),
        compensate2_meth=c("compensateFromMatrix", "doNothing"),
        qc_meth=c("qualityControlPeacoQC", "qualityControlFlowAI"),
        preTransformForQC = c(TRUE, FALSE))
    
    comb <- buildCombMatrix(alternatives)
    comb <- comb[ comb$compensate1_meth != comb$compensate2_meth ,]
    comb <- comb[ comb$compensate1_meth == "compensateFromMatrix" & 
                      comb$qc_meth == "qualityControlPeacoQC" &
                      comb$preTransformForQC == TRUE |
                      comb$compensate1_meth == "doNothing" & 
                      comb$qc_meth == "qualityControlFlowAI" &
                      comb$preTransformForQC == FALSE, ]
    #(comb)
    
    res <- runPipeline(datasets = sampleFiles[-sampleFileIndicesToRemove], 
                       alternatives = alternatives, 
                       pipelineDef = pip,
                       comb = comb,
                       nthreads = backend,
                       output.prefix = paste0(outputDir, "/"),
                       debug = FALSE,
                       skipErrors = TRUE)
    
    res
}

if (executePipeComp) {
    res <- runPipeComp(rawDataDir = rawDataDir, 
                       gatingResultsDir = resultsDir,
                       compensationDir = compensationDir,
                       outputDir = pipeCompOutputDir,
                       backend = bp,
                       sampleFileIndicesToRemove = 3) # issue with flowAI on #3 
} 
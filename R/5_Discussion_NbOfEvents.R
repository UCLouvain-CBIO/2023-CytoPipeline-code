require(CytoPipeline)
require(reshape2)
require(ggplot2)
require(patchwork)

# Figure 14

collectNbOfEvents <- function( 
        experimentName,
        path = ".",
        whichSampleFiles) {
    
    pipL <- CytoPipeline::buildCytoPipelineFromCache(
        experimentName = experimentName,
        path = path
    )
    
    if (missing(whichSampleFiles)) {
        whichSampleFiles <- CytoPipeline::sampleFiles(pipL)
    }
    
    nEventPerSampleList <- list()
    allStepNames <- c()
    for (s in seq_along(whichSampleFiles)) {
        message("Collecting nb of events for sample file ", whichSampleFiles[s], 
                "...")
        objInfos <- CytoPipeline::getCytoPipelineObjectInfos(
            pipL,
            whichQueue = "pre-processing",
            sampleFile = whichSampleFiles[s])
        objInfos <- objInfos[objInfos[,"ObjectClass"] == "flowFrame",]
        
        nEventPerSampleList[[s]] <- lapply(
            objInfos[,"ObjectName"],
            FUN = function(objName) {
                message("Reading object ", objName, "...")
               ff <- CytoPipeline::getCytoPipelineFlowFrame(
                   pipL,
                   path = path,
                   whichQueue = "pre-processing",
                   sampleFile = whichSampleFiles[s],
                   objectName = objName)
               flowCore::nrow(ff)})
        
        stepNames <- vapply(objInfos[,"ObjectName"],
                            FUN = function(str){
                                gsub(x = str,
                                     pattern = "_obj",
                                     replacement = "")
                            },
                            FUN.VALUE = character(length = 1))   
        
        names(nEventPerSampleList[[s]]) <- stepNames
        
        allStepNames <- union(allStepNames, stepNames)
    }
    
    nSampleFiles <- length(whichSampleFiles)
    nAllSteps <- length(allStepNames)
    eventNbs <- matrix(rep(NA, nSampleFiles * nAllSteps),
                       nrow = nSampleFiles)
    rownames(eventNbs) <- as.character(whichSampleFiles)
    colnames(eventNbs) <- allStepNames
    
    for (s in seq_along(whichSampleFiles)) {
        stepNames <- names(nEventPerSampleList[[s]])
        for (st in seq_along(stepNames)) {
            eventNbs[as.character(whichSampleFiles)[s],
                     stepNames[st]] <- 
                nEventPerSampleList[[s]][[st]]
        }
    }
    
    eventNbs
    
}

selectedExpName <- "HBVMouse_PQC"
selectedSamples <- c("D91_C07.fcs", "D93_A05.fcs", "D91_D03.fcs")

eventNbPeacoQC <- collectNbOfEvents( 
    experimentName = selectedExpName,
    whichSampleFiles = selectedSamples
)

eventFracPeacoQC <- t(apply(
    eventNbPeacoQC,
    MARGIN = 1,
    FUN = function(line) {
        if (length(line) == 0 || is.na(line[1])) {
            as.numeric(rep(NA, length(line)))
        } else {
            line/line[1]
        }
    }
))

selectedExpName <- "HBVMouse_flowAI"

eventNbFlowAI <- collectNbOfEvents( 
    experimentName = selectedExpName,
    whichSampleFiles = selectedSamples
)

eventFracFlowAI <- t(apply(
    eventNbFlowAI,
    MARGIN = 1,
    FUN = function(line) {
        if (length(line) == 0 || is.na(line[1])) {
            as.numeric(rep(NA, length(line)))
        } else {
            line/line[1]
        }
    }
))

stepNames <- colnames(eventFracPeacoQC)
DFPeacoQC <- as.data.frame(eventFracPeacoQC)
DFPeacoQC$sampleFile <- rownames(DFPeacoQC)
DFPeacoQC2 <- reshape(data = DFPeacoQC,
               direction = "long",
               v.names = "eventFrac",
               varying = stepNames, 
               timevar = "step",
               times = stepNames)

DFPeacoQC2$step <- factor(DFPeacoQC2$step,
                          levels = stepNames)

p1 <- ggplot(DFPeacoQC2,
             mapping = aes(x = step,
                           y = eventFrac,
                           group = sampleFile,
                           col = sampleFile)) +
    geom_line() +
    scale_y_continuous(limits=c(0,1)) + 
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none") +
    labs(y = "fraction of events kept", 
         title = "PeacoQC-based pipeline")
        
stepNames <- colnames(eventFracFlowAI)
DFFlowAI <- as.data.frame(eventFracFlowAI)
DFFlowAI$sampleFile <- rownames(DFFlowAI)
DFFlowAI2 <- reshape(data = DFFlowAI,
                      direction = "long",
                      v.names = "eventFrac",
                      varying = stepNames, 
                      timevar = "step",
                      times = stepNames)

DFFlowAI2$step <- factor(DFFlowAI2$step,
                          levels = stepNames)

p2 <- ggplot(DFFlowAI2,
             mapping = aes(x = step,
                           y = eventFrac,
                           group = sampleFile,
                           col = sampleFile)) +
    geom_line() +
    scale_y_continuous(limits=c(0,1)) + 
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y = "fraction of events kept", 
         title = "flowAI-based pipeline")


p1+p2


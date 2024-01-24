# plot some benchmark results

require("ggplot2")
require("ggrepel")
require("patchwork")
require("reshape2")
require("CytoPipelineGUI")

pipeCompOutputDir <- "./pipeComp_output"
resultsDir <- "."
rawDataDir <- "./rawData"

res <- readRDS(file = file.path(pipeCompOutputDir, "aggregated.rds"))

#res$evaluation

# figure 10

benchmarkDFLongFmt <- 
    res$evaluation$remove_dead_cells[,c("dataset",
                                        "qc_meth",
                                        "nEvents",
                                        "nGTEvents",
                                        "sensitivity",
                                        "specificity",
                                        "precision",
                                        "recall")]

benchmarkDFLongFmt$qc_meth <- 
    substr(benchmarkDFLongFmt$qc_meth,
           15, nchar(benchmarkDFLongFmt$qc_meth))

benchmarkDFPipeComp <- 
    reshape(data = benchmarkDFLongFmt,
            idvar = "dataset",
            timevar = "qc_meth",
            direction = "wide",
            sep = "")



summaryDF <- melt(benchmarkDFPipeComp,
                  id.vars = c("dataset"),
                  variable.name = "indicator_pipe",
                  measure.vars = c("nEventsPeacoQC", "nEventsFlowAI",
                                   "sensitivityPeacoQC", "sensitivityFlowAI",
                                   "specificityPeacoQC", "specificityFlowAI",
                                   "precisionPeacoQC", "precisionFlowAI",
                                   "recallPeacoQC", "recallFlowAI"))

summaryDF$pipeline <- ifelse(grepl("PeacoQC", summaryDF$indicator_pipe),
                             "PeacoQC",
                             "flowAI")
summaryDF$pipeline <- factor(summaryDF$pipeline, levels = c("PeacoQC", "flowAI"))

summaryDF$indicator <- ifelse(summaryDF$pipeline == "PeacoQC",
                              gsub("PeacoQC", "", summaryDF$indicator_pipe),
                              gsub("FlowAI", "", summaryDF$indicator_pipe))

summaryDF$indicator <- factor(summaryDF$indicator, levels = c("sensitivity",
                                                              "specificity",
                                                              "precision",
                                                              "recall",
                                                              "nEvents"))

summaryDF$indicator_pipe <- NULL

summaryDF


ggplot(summaryDF[summaryDF$indicator != "nEvents",], aes(x = pipeline, y = value)) + 
    ggplot2::facet_wrap(~ indicator) + 
    geom_boxplot()


# figure 11

summaryDF2 <- reshape(data = summaryDF,
                      idvar = c("dataset","indicator"),
                      timevar = "pipeline",
                      direction = "wide",
                      sep = "")

colnames(summaryDF2) <- c("dataset", "indicator", "PeacoQC", "flowAI")

useRepel <- TRUE

geom_text_repel_high_overlap <- function(...){
    geom_text_repel(..., max.overlaps = 100)
}

toDisplay <- c("D93_A05", "D91_C07", "D91_D03")

summaryDF2$label <- summaryDF2$dataset
summaryDF2$label[!summaryDF2$dataset %in% toDisplay] <- ""
summaryDF2$highlight <- "false"
summaryDF2$highlight[summaryDF2$dataset %in% toDisplay] <- "true"
summaryDF2 <- 
    summaryDF2[order(summaryDF2$highlight),]


if (useRepel){
    geom_text_FUN <- geom_text_repel_high_overlap
} else {
    geom_text_FUN <- geom_text
}

highlightCol <- "red"
normalCol <- "black"

p <- 
    ggplot(data = summaryDF2[summaryDF2$indicator != "nEvents",], 
           mapping = aes(x = PeacoQC,
                         y = flowAI)) + 
    geom_point(mapping = aes(colour = highlight)) +
    geom_text_FUN(mapping = aes(label = label), hjust=0.5, vjust=1,
                  size = 2.5) +
    geom_abline(intercept = 0., slope = 1., colour = "blue", linetype = 2) +
    coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) + 
    scale_colour_manual(values = c("false" = normalCol, 
                                   "true" = highlightCol)) +
    xlab("PeacoQC pipeline") + ylab("flowAI pipeline") + 
    theme(legend.position="none") + 
    ggplot2::facet_wrap(~ indicator)
p   

# figure 12.A
selectedExpName <- "HBVMouse_flowAI"
selectedSampleFile <- file.path(
    rawDataDir,
    "D91_C07.fcs")
p1 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile,
    flowFrameName = "remove_debris_obj",
    path = resultsDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "Comp-APC-Cy7-A : Live & Dead",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p2 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile,
    flowFrameName = "remove_dead_cells_obj",
    path = resultsDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "Comp-APC-Cy7-A : Live & Dead",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p1+p2


# figure 12.B
selectedExpName <- "HBVMouse_PQC"
selectedSampleFile <- file.path(
    rawDataDir,
    "D91_C07.fcs")
p1 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile,
    flowFrameName = "remove_debris_obj",
    path = resultsDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "Comp-APC-Cy7-A : Live & Dead",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p2 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile,
    flowFrameName = "remove_dead_cells_obj",
    path = resultsDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "Comp-APC-Cy7-A : Live & Dead",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p1+p2


# Figure S5

selectedExpName <- "HBVMouse_PQC"
selectedSampleFile <- file.path(
    rawDataDir,
    "D93_A05.fcs")

p1 <- CytoPipelineGUI:::plotDiffFlowFrame(
    experimentNameFrom = selectedExpName,
    whichQueueFrom = "pre-processing",
    sampleFileFrom = selectedSampleFile,
    flowFrameNameFrom = "perform_QC_obj",
    xChannelLabelFrom = "FSC-A : NA",
    yChannelLabelFrom = "SSC-A : NA",
    path = resultsDir,
    experimentNameTo = selectedExpName,
    whichQueueTo = "pre-processing",
    sampleFileTo = selectedSampleFile,
    flowFrameNameTo = "remove_dead_cells_obj",
    xChannelLabelTo = "FSC-A : NA",
    yChannelLabelTo = "SSC-A : NA",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))


selectedExpName <- "HBVMouse_flowAI"
p2 <- CytoPipelineGUI:::plotDiffFlowFrame(
    experimentNameFrom = selectedExpName,
    whichQueueFrom = "pre-processing",
    sampleFileFrom = selectedSampleFile,
    flowFrameNameFrom = "compensate_obj",
    xChannelLabelFrom = "FSC-A : NA",
    yChannelLabelFrom = "SSC-A : NA",
    path = resultsDir,
    experimentNameTo = selectedExpName,
    whichQueueTo = "pre-processing",
    sampleFileTo = selectedSampleFile,
    flowFrameNameTo = "remove_dead_cells_obj",
    xChannelLabelTo = "FSC-A : NA",
    yChannelLabelTo = "SSC-A : NA",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

selectedExpName <- "HBVMouse_FJ"
p3 <- CytoPipelineGUI:::plotDiffFlowFrame(
    experimentNameFrom = selectedExpName,
    whichQueueFrom = "pre-processing",
    sampleFileFrom = selectedSampleFile,
    flowFrameNameFrom = "time_gate_obj",
    xChannelLabelFrom = "FSC-A : NA",
    yChannelLabelFrom = "SSC-A : NA",
    path = resultsDir,
    experimentNameTo = selectedExpName,
    whichQueueTo = "pre-processing",
    sampleFileTo = selectedSampleFile,
    flowFrameNameTo = "lymphocytes_gate_obj",
    xChannelLabelTo = "FSC-A : NA",
    yChannelLabelTo = "SSC-A : NA",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p1+p2+p3


# Figure S6

selectedExpName <- "HBVMouse_PQC"
selectedSampleFile <- file.path(
    rawDataDir,
    "D91_D03.fcs")

p1 <- CytoPipelineGUI:::plotDiffFlowFrame(
    experimentNameFrom = selectedExpName,
    whichQueueFrom = "pre-processing",
    sampleFileFrom = selectedSampleFile,
    flowFrameNameFrom = "compensate_obj",
    xChannelLabelFrom = "Time : NA",
    yChannelLabelFrom = "FSC-A : NA",
    path = resultsDir,
    experimentNameTo = selectedExpName,
    whichQueueTo = "pre-processing",
    sampleFileTo = selectedSampleFile,
    flowFrameNameTo = "perform_QC_obj",
    xChannelLabelTo = "Time : NA",
    yChannelLabelTo = "FSC-A : NA",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = FALSE,
    linearRange = NULL)


selectedExpName <- "HBVMouse_flowAI"
p2 <- CytoPipelineGUI:::plotDiffFlowFrame(
    experimentNameFrom = selectedExpName,
    whichQueueFrom = "pre-processing",
    sampleFileFrom = selectedSampleFile,
    flowFrameNameFrom = "remove_margins_obj",
    xChannelLabelFrom = "Time : NA",
    yChannelLabelFrom = "FSC-A : NA",
    path = resultsDir,
    experimentNameTo = selectedExpName,
    whichQueueTo = "pre-processing",
    sampleFileTo = selectedSampleFile,
    flowFrameNameTo = "perform_QC_obj",
    xChannelLabelTo = "Time : NA",
    yChannelLabelTo = "FSC-A : NA",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = FALSE,
    linearRange = NULL)

selectedExpName <- "HBVMouse_FJ"
p3 <- CytoPipelineGUI:::plotDiffFlowFrame(
    experimentNameFrom = selectedExpName,
    whichQueueFrom = "pre-processing",
    sampleFileFrom = selectedSampleFile,
    flowFrameNameFrom = "flowframe_read_obj",
    xChannelLabelFrom = "Time : NA",
    yChannelLabelFrom = "FSC-A : NA",
    path = resultsDir,
    experimentNameTo = selectedExpName,
    whichQueueTo = "pre-processing",
    sampleFileTo = selectedSampleFile,
    flowFrameNameTo = "time_gate_obj",
    xChannelLabelTo = "Time : NA",
    yChannelLabelTo = "FSC-A : NA",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = FALSE,
    linearRange = NULL)

p1+p2+p3




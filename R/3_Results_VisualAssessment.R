require(CytoPipeline)
require(CytoPipelineGUI)
require(patchwork)


if (!exists("pipelineList"))
    stop("PipelineList object not found => execute script #1 first!")


# USE CASE #1


CytoPipeline::plotCytoPipelineProcessingQueue(
    x <- pipelineList[[1]],
    whichQueue = "pre-processing",
    sampleFile = 1,
    path = outputDir,
    title = ""
)
    

# USE CASE #2

selectedExpName <- "HBVMouse_PQC"
selectedSampleFile <- file.path(
    rawDataDir,
    "D91_G01.fcs")

p1 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile,
    flowFrameName = "remove_margins_obj",
    path = outputDir,
    xChannelLabel = "BUV805-A : CD8",
    yChannelLabel = "BB700-A : CD38",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p2 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile,
    flowFrameName = "compensate_obj",
    path = outputDir,
    xChannelLabel = "Comp-BUV805-A : CD8",
    yChannelLabel = "Comp-BB700-A : CD38",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p1+p2

# USE CASE #3

selectedExpName1 <- "HBVMouse_PQC_nC3"
selectedExpName2 <- "HBVMouse_PQC"
selectedSampleFile <- file.path(
    rawDataDir,
    "D91_G01.fcs")

p1 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName1,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile,
    flowFrameName = "remove_debris_obj",
    path = outputDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "SSC-A : NA",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p2 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName2,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile,
    flowFrameName = "remove_debris_obj",
    path = outputDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "SSC-A : NA",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p3 <- CytoPipelineGUI:::plotDiffFlowFrame(
    experimentNameFrom = selectedExpName1,
    whichQueueFrom = "pre-processing",
    sampleFileFrom = selectedSampleFile,
    flowFrameNameFrom = "remove_debris_obj",
    xChannelLabelFrom = "FSC-A : NA",
    yChannelLabelFrom = "SSC-A : NA",
    path = outputDir,
    experimentNameTo = selectedExpName2,
    whichQueueTo = "pre-processing",
    sampleFileTo = selectedSampleFile,
    flowFrameNameTo = "remove_debris_obj",
    xChannelLabelTo = "FSC-A : NA",
    yChannelLabelTo = "SSC-A : NA",
    useAllCells = FALSE,
    nDisplayCells = 10000,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p1+p2+p3

# USE CASE #4

selectedExpName1 <- "HBVMouse_PQC"
selectedExpName2 <- "HBVMouse_flowAI"
selectedSampleFile <- file.path(
    rawDataDir,
    "D91_A01.fcs")

p1 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName1,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile,
    flowFrameName = "remove_dead_cells_obj",
    path = outputDir,
    xChannelLabel = "Time : NA",
    yChannelLabel = "FSC-A : NA",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = FALSE,
    linearRange = NULL)

p2 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName2,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile,
    flowFrameName = "remove_dead_cells_obj",
    path = outputDir,
    xChannelLabel = "Time : NA",
    yChannelLabel = "FSC-A : NA",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = FALSE,
    linearRange = NULL)

p3 <- CytoPipelineGUI:::plotDiffFlowFrame(
    experimentNameFrom = selectedExpName1,
    whichQueueFrom = "pre-processing",
    sampleFileFrom = selectedSampleFile,
    flowFrameNameFrom = "remove_dead_cells_obj",
    xChannelLabelFrom = "Time : NA",
    yChannelLabelFrom = "FSC-A : NA",
    path = outputDir,
    experimentNameTo = selectedExpName2,
    whichQueueTo = "pre-processing",
    sampleFileTo = selectedSampleFile,
    flowFrameNameTo = "remove_dead_cells_obj",
    xChannelLabelTo = "Time : NA",
    yChannelLabelTo = "FSC-A : NA",
    useAllCells = FALSE,
    nDisplayCells = 10000,
    useFixedLinearRange = FALSE,
    linearRange = NULL)

p1+p2+p3

# USE CASE #5

selectedExpName <- "HBVMouse_PQC"
selectedSampleFile1 <- file.path(
    rawDataDir,
    "D91_A01.fcs")
selectedSampleFile2 <- file.path(
    rawDataDir,
    "D93_B05.fcs")

p1 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile1,
    flowFrameName = "remove_dead_cells_obj",
    path = outputDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "Comp-APC-Cy7-A : Live & Dead",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p2 <- CytoPipelineGUI:::plotSelectedFlowFrame(
    experimentName = selectedExpName,
    whichQueue = "pre-processing",
    sampleFile = selectedSampleFile2,
    flowFrameName = "remove_dead_cells_obj",
    path = outputDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "Comp-APC-Cy7-A : Live & Dead",
    useAllCells = TRUE,
    nDisplayCells = 0,
    useFixedLinearRange = TRUE,
    linearRange = c(-100, 262144))

p1+p2

# USE CASE #6
CytoPipelineGUI::ScaleTransformApp()

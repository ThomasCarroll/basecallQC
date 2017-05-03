readPattern <- function(bcl2fastqparams){
  c(readlengths(bcl2fastqparams)[1],indexlengths(bcl2fastqparams)[1],
    indexlengths(bcl2fastqparams)[2],readlengths(bcl2fastqparams)[2])
}
readBCLStatFile <- function(read.filename){
  read.filename <- file(read.filename, "rb")
  clusterNumber <- readBin(read.filename,"int",size=4,n = 1)
  ACI <- readBin(read.filename,"double",size=8,n = 1)
  AIAOverAIAA <- readBin(read.filename,"double",size=8,n = 1)
  AICOverAICA <- readBin(read.filename,"double",size=8,n = 1)
  AIGOverAIGA <- readBin(read.filename,"double",size=8,n = 1)
  AITOverAITA <- readBin(read.filename,"double",size=8,n = 1)
  AIAOverA <- readBin(read.filename,"double",size=8,n = 1)
  AICOverC <- readBin(read.filename,"double",size=8,n = 1)
  AIGOverG <- readBin(read.filename,"double",size=8,n = 1)
  AITOverT <- readBin(read.filename,"double",size=8,n = 1)
  NCA <- readBin(read.filename,"int",size=4,n = 1)
  NCC <- readBin(read.filename,"int",size=4,n = 1)
  NCG <- readBin(read.filename,"int",size=4,n = 1)
  NCT <- readBin(read.filename,"int",size=4,n = 1)
  NCX <- readBin(read.filename,"int",size=4,n = 1)
  NCIA <- readBin(read.filename,"int",size=4,n = 1)
  NCIC <- readBin(read.filename,"int",size=4,n = 1)
  NCIG <- readBin(read.filename,"int",size=4,n = 1)
  NCIT <- readBin(read.filename,"int",size=4,n = 1)
  
  
  bclClusterTileStats <- c(clusterNumber,ACI,AIAOverAIAA,AICOverAICA,AIGOverAIGA,
                           AITOverAITA,AIAOverA,AICOverC,AIGOverG,
                           AITOverT,NCA,NCC,NCG,NCT,NCX,NCIA,NCIC,NCIG,NCIT)
  close(read.filename)
  return(bclClusterTileStats)
}
readBCLStatFiles <- function(read.filenames){
  bclStatList <- lapply(read.filenames,readBCLStatFile)
  bclStatMat <- do.call(rbind,bclStatList)
}
readExtractionMetrics <- function(extractionMetricsBin){
  bytesOfFile <- file.info(extractionMetricsBin)$size
  read.filename <- file(extractionMetricsBin, "rb")
  bytesRead <- 0
  k <- 1
  versionNumber <- readBin(read.filename,"int",size=1,n = 1)
  recordSize <- readBin(read.filename,"int",size=1,n = 1,signed = FALSE)
  bytesRead <- bytesRead+1+1
  extractionMetricsOut <- list()
  while(bytesRead < bytesOfFile){
    laneNumber <- readBin(read.filename,"int",size=2,n = 1)
    tileNumber <- readBin(read.filename,"int",size=2,n = 1)
    cycleNumber <- readBin(read.filename,"int",size=2,n = 1)
    focusForChannelA <- readBin(read.filename,"double",size=4,n = 1)
    focusForChannelC <- readBin(read.filename,"double",size=4,n = 1)
    focusForChannelG <- readBin(read.filename,"double",size=4,n = 1)
    focusForChannelT <- readBin(read.filename,"double",size=4,n = 1)
    maxIntensityForChannelA <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    maxIntensityForChannelC <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    maxIntensityForChannelG <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    maxIntensityForChannelT <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    timeStamp <- readBin(read.filename,"int",size=8,n = 1)
    bytesRead <- bytesRead+(2*3)+(4*4)+(2*4)+8
    extractionMetricsOut[[k]] <- c(
      laneNumber,tileNumber,cycleNumber,
      focusForChannelA,focusForChannelC,
      focusForChannelG,focusForChannelT,
      maxIntensityForChannelA,maxIntensityForChannelC,
      maxIntensityForChannelG,maxIntensityForChannelT,
      timeStamp
    )
    k <- k+1
  }
  extractionMetricsOutMat <- do.call(rbind,extractionMetricsOut)
  colnames(extractionMetricsOutMat) <-  c(
    "laneNumber","tileNumber","cycleNumber",
    "focusForChannelA","focusForChannelC",
    "focusForChannelG","focusForChannelT",
    "maxIntensityForChannelA","maxIntensityForChannelC",
    "maxIntensityForChannelG","maxIntensityForChannelT",
    "timeStamp"
  )
  close(read.filename)
  extractionMetricsOutMat <- tbl_df(extractionMetricsOutMat)
}
readIndexMetrics <- function(indexMetricsBin){
  bytesOfFile <- file.info(indexMetricsBin)$size
  indexMetricsBin <- file(indexMetricsBin, "rb")
  
  bytesRead <- 0
  versionNumber <- readBin(indexMetricsBin,"int",size=1,n = 1)
  bytesRead <- bytesRead+1
  k <- 1
  IndexMetricsOut <- list()
  while(bytesRead < bytesOfFile){
    laneNumber <- readBin(indexMetricsBin,"int",size=2,n = 1)
    tileNumber <- readBin(indexMetricsBin,"int",size=2,n = 1)
    readNumber <- readBin(indexMetricsBin,"int",size=2,n = 1)
    indexNameLength <- readBin(indexMetricsBin,"int",size=2,n = 1)
    indexName <- rawToChar(readBin(indexMetricsBin,"raw",n=indexNameLength))
    identifiedAsIndex <- readBin(indexMetricsBin,"int",size=4,n = 1)
    sampleNameLength <- readBin(indexMetricsBin,"int",size=2,n = 1)
    sampleName <- rawToChar(readBin(indexMetricsBin,"raw",n=sampleNameLength))
    projectNameLength <- readBin(indexMetricsBin,"int",size=2,n = 1)
    projectName <- rawToChar(readBin(indexMetricsBin,"raw",n=projectNameLength))
    bytesRead <- bytesRead + 2+2+2+2+indexNameLength+4+2+sampleNameLength+2+projectNameLength
    IndexMetricsOut[[k]] <- c(
      laneNumber,tileNumber,readNumber,indexName,identifiedAsIndex,sampleName, projectName
    )
    k  <- k+1
  }
  close(indexMetricsBin)
  IndexMetricsFrame <- do.call(rbind,IndexMetricsOut)
  colnames(IndexMetricsFrame) <- c(
    "laneNumber","tileNumber","readNumber","indexName","identifiedAsIndex","sampleName","projectName"
  )
  tbl_df(IndexMetricsFrame)
}
readQMetrics <- function(qMetricsBin){
  bytesOfFile <- file.info(qMetricsBin)$size
  read.filename <- file(qMetricsBin, "rb")
  bytesRead <- 0
  versionNumber <- readBin(read.filename,"int",size=1,n = 1)
  recordSize <- readBin(read.filename,"int",size=1,n = 1,signed = FALSE)
  hasBins <- readBin(read.filename,"logical",size=1,n = 1)
  numberOfBins <- readBin(read.filename,"int",size=1,n = 1,signed = FALSE)
  lowEndsOfBins <- readBin(read.filename,"int",size=1,n = numberOfBins,signed = FALSE)
  highEndsOfBins <- readBin(read.filename,"int",size=1,n = numberOfBins,signed = FALSE)
  valueOfBins <- readBin(read.filename,"int",size=1,n = numberOfBins,signed = FALSE)
  bytesRead <- bytesRead+1+1+1+1+numberOfBins*3
  k <- 1
  qMetricsOut <- list()
  while(bytesRead < bytesOfFile){
    Lane <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    Tile <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    Cycle <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    qhist <- readBin(read.filename,"integer",size=4,n = 50)
    qMetricsOut[[k]] <-
      c(Lane,Tile,Cycle,qhist)
    bytesRead <- bytesRead+2+2+2+4*50
    k <- k+1
  }
  
  close(read.filename)
  qMetricsFrame <- do.call(rbind,qMetricsOut)
  colnames(qMetricsFrame) <- c("Lane","Tile","Cycle",paste0("Q",1:50))
  tbl_df(qMetricsFrame)
}
readCorrectedIntMetrics <- function(correctedIntMetricsBin){
  bytesOfFile <- file.info(correctedIntMetricsBin)$size
  read.filename <- file(correctedIntMetricsBin, "rb")
  bytesRead <- 0
  k <- 1
  versionNumber <- readBin(read.filename,"int",size=1,n = 1)
  recordSize <- readBin(read.filename,"int",size=1,n = 1,signed = FALSE)
  bytesRead <- bytesRead+1+1
  correctedIntMetricsOut <- list()
  while(bytesRead < bytesOfFile){
    laneNumber <- readBin(read.filename,"int",size=2,n = 1)
    tileNumber <- readBin(read.filename,"int",size=2,n = 1)
    cycleNumber <- readBin(read.filename,"int",size=2,n = 1)
    averageCycleIntensity <- readBin(read.filename,"int",size=2,n = 1)
    averageCorrectedIntensityForChannelA <- readBin(read.filename,"int",size=2,n = 1)
    averageCorrectedIntensityForChannelC <- readBin(read.filename,"int",size=2,n = 1)
    averageCorrectedIntensityForChannelG <- readBin(read.filename,"int",size=2,n = 1)
    averageCorrectedIntensityForChannelT <- readBin(read.filename,"int",size=2,n = 1)
    averageCorrectedIntForCalledClustersForBaseA <- readBin(read.filename,"int",size=2,n = 1)
    averageCorrectedIntForCalledClustersForBaseC <- readBin(read.filename,"int",size=2,n = 1)
    averageCorrectedIntForCalledClustersForBaseG <- readBin(read.filename,"int",size=2,n = 1)
    averageCorrectedIntForCalledClustersForBaseT <- readBin(read.filename,"int",size=2,n = 1)
    numberOfBaseCallsForNoCall  <- readBin(read.filename,"int",size=4,n = 1)
    numberOfBaseCallsForBaseA  <- readBin(read.filename,"int",size=4,n = 1)
    numberOfBaseCallsForBaseC  <- readBin(read.filename,"int",size=4,n = 1)
    numberOfBaseCallsForBaseG  <- readBin(read.filename,"int",size=4,n = 1)
    numberOfBaseCallsForBaseT  <- readBin(read.filename,"int",size=4,n = 1)
    signalToNoiseRatio   <- readBin(read.filename,"double",size=4,n = 1)
    bytesRead <- bytesRead+(2*3)+(2*9)+(4*5)+4
    correctedIntMetricsOut[[k]] <- c(laneNumber,tileNumber,cycleNumber,averageCycleIntensity,
                                     averageCorrectedIntensityForChannelA,averageCorrectedIntensityForChannelC,
                                     averageCorrectedIntensityForChannelG,averageCorrectedIntensityForChannelT,
                                     averageCorrectedIntForCalledClustersForBaseA,averageCorrectedIntForCalledClustersForBaseC,
                                     averageCorrectedIntForCalledClustersForBaseG,averageCorrectedIntForCalledClustersForBaseT,
                                     numberOfBaseCallsForNoCall,
                                     numberOfBaseCallsForBaseA,numberOfBaseCallsForBaseC,
                                     numberOfBaseCallsForBaseG,numberOfBaseCallsForBaseT)
    k <- k+1
  }
  close(read.filename)
  correctedIntMetricsMat <- do.call(rbind,correctedIntMetricsOut)
}
readImageMetrics <- function(imageMetricsBin){
  bytesOfFile <- file.info(imageMetricsBin)$size
  read.filename <- file(imageMetricsBin, "rb")
  bytesRead <- 0
  k <- 1
  versionNumber <- readBin(read.filename,"int",size=1,n = 1)
  recordSize <- readBin(read.filename,"int",size=1,n = 1,signed = FALSE)
  bytesRead <- bytesRead+1+1
  imageMetricsOut <- list()
  while(bytesRead < bytesOfFile){
    laneNumber <- readBin(read.filename,"int",size=2,n = 1)
    tileNumber <- readBin(read.filename,"int",size=2,n = 1)
    cycleNumber <- readBin(read.filename,"int",size=2,n = 1)
    channelNumber <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    minimumContrast <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    maximumContrast <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    bytesRead <- bytesRead+(2*3)+(2*3)
    imageMetricsOut[[k]] <- c(laneNumber,tileNumber,cycleNumber,
                              channelNumber,minimumContrast,maximumContrast)
    k <- k+1
  }
  close(read.filename)
  imageMetricsMat <- do.call(rbind,imageMetricsOut)
}
readTileMetrics <- function(tileMetricsBin){
  bytesOfFile <- file.info(tileMetricsBin)$size
  read.filename <- file(tileMetricsBin, "rb")
  bytesRead <- 0
  k <- 1
  versionNumber <- readBin(read.filename,"int",size=1,n = 1)
  recordSize <- readBin(read.filename,"int",size=1,n = 1,signed = FALSE)
  bytesRead <- bytesRead+1+1
  tileMetricsOut <- list()
  while(bytesRead < bytesOfFile){
    laneNumber <- readBin(read.filename,"int",size=2,n = 1)
    tileNumber <- readBin(read.filename,"int",size=2,n = 1)
    codeTile <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    valueTile <- readBin(read.filename,"double",size=4,n = 1)
    bytesRead <- bytesRead+(2*2)+2+4
    tileMetricsOut[[k]] <- c(laneNumber,tileNumber,
                             codeTile,valueTile)
    k <- k+1
  }
  close(read.filename)
  tileMetricsMat <- do.call(rbind,tileMetricsOut)
  colnames(tileMetricsMat) <-  c(
    "laneNumber","tileNumber",
    "code","value"
  )
  tbl_df(tileMetricsMat) %>%
    mutate(code = str_c("c",code)) %>%
    group_by(laneNumber,tileNumber,code) %>%
    summarise_all(max) %>%
    ungroup() %>%
    mutate(Lane=laneNumber,Tile=tileNumber) %>% 
    mutate(code=stringr::str_replace_all(code,"c100","ClusterDensity")) %>% 
    mutate(code=stringr::str_replace_all(code,"c101","ClusterDensityPF")) %>% 
    mutate(code=stringr::str_replace_all(code,"c102","NumberOfClusters")) %>%        
    mutate(code=stringr::str_replace_all(code,"c103","NumberOfClustersPF")) %>%                       
    mutate(code=stringr::str_replace_all(code,"c200",
                                         paste0("PhasingFor",
                                                names(readPattern(bcl2fastqparams))[readPattern(bcl2fastqparams) != 0][1]))) %>% 
    mutate(code=stringr::str_replace_all(code,"c201",
                                         paste0("PrePhasingFor",
                                                names(readPattern(bcl2fastqparams))[readPattern(bcl2fastqparams) != 0][1]))) %>%    
    mutate(code=stringr::str_replace_all(code,"c202",
                                         paste0("PhasingFor",
                                                names(readPattern(bcl2fastqparams))[readPattern(bcl2fastqparams) != 0][2]))) %>% 
    mutate(code=stringr::str_replace_all(code,"c203",
                                         paste0("PrePhasingFor",
                                                names(readPattern(bcl2fastqparams))[readPattern(bcl2fastqparams) != 0][2]))) %>%
    mutate(code=stringr::str_replace_all(code,"c204",
                                         paste0("PhasingFor",
                                                names(readPattern(bcl2fastqparams))[readPattern(bcl2fastqparams) != 0][3]))) %>% 
    mutate(code=stringr::str_replace_all(code,"c205",
                                         paste0("PrePhasingFor",
                                                names(readPattern(bcl2fastqparams))[readPattern(bcl2fastqparams) != 0][3]))) %>% 
    mutate(code=stringr::str_replace_all(code,"c206",
                                         paste0("PhasingFor",
                                                names(readPattern(bcl2fastqparams))[readPattern(bcl2fastqparams) != 0][4]))) %>% 
    mutate(code=stringr::str_replace_all(code,"c207",
                                         paste0("PrePhasingFor",
                                                names(readPattern(bcl2fastqparams))[readPattern(bcl2fastqparams) != 0][4]))) %>%
    filter(code!="PrePhasingForNA" & code!="PhasingForNA" ) %>% 
    spread(code,value) %>%
    dplyr::select(-laneNumber,-tileNumber,
                          -starts_with("c",ignore.case = FALSE))
}
readPhasingMetrics <- function(phasingMetricsBin){
  bytesOfFile <- file.info(phasingMetricsBin)$size
  read.filename <- file(phasingMetricsBin, "rb")
  bytesRead <- 0
  k <- 1
  versionNumber <- readBin(read.filename,"int",size=1,n = 1)
  recordSize <- readBin(read.filename,"int",size=1,n = 1,signed = FALSE)
  bytesRead <- bytesRead+1+1
  phasingMetricsOut <- list()
  while(bytesRead < bytesOfFile){
    laneNumber <- readBin(read.filename,"int",size=2,n = 1)
    tileNumber <- readBin(read.filename,"int",size=2,n = 1)
    cycleNumber <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    phasingWeight  <- readBin(read.filename,"double",size=4,n = 1)
    prephasingWeight  <- readBin(read.filename,"double",size=4,n = 1)
    bytesRead <- bytesRead+(2*3)+(4*2)
    phasingMetricsOut[[k]] <- c(laneNumber,tileNumber,cycleNumber,
                                phasingWeight,prephasingWeight)
    k <- k+1
  }
  close(read.filename)
  phasingMetricsMat <- do.call(rbind,phasingMetricsOut)
}
readErrorMetrics <- function(errorMetricsBin){
  bytesOfFile <- file.info(errorMetricsBin)$size
  read.filename <- file(errorMetricsBin, "rb")
  bytesRead <- 0
  k <- 1
  versionNumber <- readBin(read.filename,"int",size=1,n = 1)
  recordSize <- readBin(read.filename,"int",size=1,n = 1,signed = FALSE)
  bytesRead <- bytesRead+1+1
  tileMetricsOut <- list()
  while(bytesRead < bytesOfFile){
    laneNumber <- readBin(read.filename,"int",size=2,n = 1)
    tileNumber <- readBin(read.filename,"int",size=2,n = 1)
    cycleNumber <- readBin(read.filename,"int",size=2,n = 1,signed = FALSE)
    errorRate <- readBin(read.filename,"double",size=4,n = 1)
    numberOfPerfectReads <- readBin(read.filename,"int",size=4,n = 1)
    numberOfReadsWithOneError <- readBin(read.filename,"int",size=4,n = 1)
    numberOfReadsWithTwoError <- readBin(read.filename,"int",size=4,n = 1)
    numberOfReadsWithThreeError <- readBin(read.filename,"int",size=4,n = 1)
    numberOfReadsWithFourError <- readBin(read.filename,"int",size=4,n = 1)
    bytesRead <- bytesRead+(2*3)+(4*6)
    tileMetricsOut[[k]] <- c(laneNumber,tileNumber,cycleNumber,
                             errorRate,numberOfPerfectReads,
                             numberOfReadsWithOneError,numberOfReadsWithTwoError,
                             numberOfReadsWithThreeError,numberOfReadsWithFourError)
    k <- k+1
  }
  close(read.filename)
  tileMetricsMat <- do.call(rbind,tileMetricsOut)
}
readInterOpsMetrics <- function(bcl2fastqparams,verbose=TRUE,
                                interopsToParse=c("TileMetricsOut","QMetricsOut")){
  interOpsMetrics <- list()
  if(any(interopsToParse %in% "TileMetricsOut")){
    if(verbose){message("Parsing InterOps Tile binary file..",appendLF = FALSE)}
    tileMetricsFile <- file.path(bcl2fastqparams@RunDir,"InterOp","TileMetricsOut.bin")
    if(file.exists(tileMetricsFile)){
      interOpsMetrics[["tileMet"]] <- readTileMetrics(tileMetricsFile)
    }
    if(verbose){message("done",appendLF = TRUE)}
  }
  if(any(interopsToParse %in% "QMetricsOut")){
    if(verbose){message("Parsing InterOps QMetrics binary file..",appendLF = FALSE)}
    qMetricsFile <- file.path(bcl2fastqparams@RunDir,"InterOp","QMetricsOut.bin")
    if(file.exists(qMetricsFile)){
      interOpsMetrics[["qMet"]] <- readQMetrics(qMetricsFile)
    }
    if(verbose){message("done",appendLF = TRUE)}
  }
  return(interOpsMetrics)
}


#' Function to parse InterOps files and generate summary reports
#'
#' Parses the InterOps binary files produced by Illumina's Real Time Analysis sofware and
#' used by Illumina's SAV sofware.
#' InterOp binary files contain information on phasing/prephsing, yield,read numbers
#' and basecalling quality score distributions per cycle.
#' This interOpsReport functions parses and summarises the InterOps files, TileMetrics.bin and QMetrics.bin, and the 
#' Stats directory XML files, ConversionStats.xml and DemultiplexingStats.xml. 
#'
#'
#' @docType methods
#' @name interOpsReport
#' @rdname interOpsReport
#'
#' @author Thomas Carroll.
#' @param bcl2fastqparams A BCL2FastQparams object.
#' @param verbose TRUE or FALSE . TRUE reports progress through file parsing.
#' @return A list of length 3 containing machine and run information, 
#' basecalling quality information and demultiplexing information.
#' 
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#'

#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#' 
#' # myRes_BCAGJ8ANXX <- interOpsReport(bcl2fastqparams,verbose=TRUE) 
#' @export
interOpsReport <- function(bcl2fastqparams,verbose=TRUE){
  if(verbose){message("Parsing XML files..",appendLF = FALSE)}
  basecallmetrics <- baseCallMetrics(bcl2fastqparams)
  demultiplexmetrics <- demultiplexMetrics(bcl2fastqparams)
  dmx <- demultiplexmetrics$demuxStatsProcessed
  conv <- basecallmetrics$convStatsProcessed
  if(verbose){message("done",appendLF = TRUE)}
  
  dmxSum <- dmx %>% 
    tbl_df %>% 
    dplyr::select(Sample,Barcode,Lane) %>%
    distinct(Sample,Barcode,Lane) %>%
    mutate(Lane=factor(as.numeric(str_replace_all(Lane,"Lane",""))))
  
  convSum <- conv %>%
    tbl_df %>%
    group_by(Sample,Lane,Filter) %>%
    summarise(YieldSum=sum(as.numeric(Yield)),
              Yield30Sum=sum(as.numeric(Yield30)),
              Reads=sum(ClusterCount)) %>%
    filter(Filter=="Pf") %>%
    group_by(Lane) %>%
    mutate(PercentageOfGreaterQ30_Bases_PF=(Yield30Sum/YieldSum)*100,
           PercentOfLane=(Reads/sum(Reads[Sample != "all"]))*100) %>%
    dplyr::select(-Filter,-Yield30Sum)
  
  dmxReport <- full_join(convSum,tbl_df(dmxSum),by = c("Sample", "Lane"))
  
  if(verbose){message("Parsing machine information..",appendLF = FALSE)}
  
  runPReport <- bcl2fastqparams@RunParameters$runParams %>% 
    dplyr::select(ComputerName,RunID,ExperimentName,
                  Flowcell,
                  ChemistryVersion,RunMode,
                  ApplicationName,ApplicationVersion,
                  RTAVersion
    ) %>%  t
  
  if(verbose){message("done",appendLF = TRUE)}
  if(verbose){message("Parsing InterOps binary files..",appendLF = FALSE)}
  interOpsMetrics <- readInterOpsMetrics(bcl2fastqparams,verbose = FALSE)
  tileMet <- interOpsMetrics$tileMet  %>%
    group_by(Lane) %>%
    summarise(meanClusterDensity=mean(ClusterDensity),
              sdClusterDensity=sd(ClusterDensity),
              meanPhasingForRead1=mean(PhasingForRead1,na.rm=TRUE),
              meanPrePhasingForRead1=mean(PrePhasingForRead1,na.rm=TRUE),
              NumberOfReads=sum(NumberOfClusters),
              NumberOfReadsPF=sum(NumberOfClustersPF),
              percent_PF_Clusters=mean(NumberOfClustersPF/NumberOfClusters)*100,
              percent_PF_ClustersSD=sd(NumberOfClustersPF/NumberOfClusters)*100,
              meanPhasingForRead2 = if(exists('PhasingForRead2', where = .)) mean(PhasingForRead2,na.rm=TRUE) else 0,
              meanPrePhasingForRead2 = if(exists('PrePhasingForRead2', where = .)) mean(PrePhasingForRead2,na.rm=TRUE) else 0
    ) %>% 
    mutate(meanClusterDensity=stringr::str_c(signif(meanClusterDensity/1000,4),"+/-",signif(sdClusterDensity/1000,4)),
           percent_PF_Clusters=stringr::str_c(signif(percent_PF_Clusters,3),"+/-",signif(percent_PF_ClustersSD,3)),
           Phasing_PrephasingRead1=stringr::str_c(signif(meanPhasingForRead1*100,3)," / ",signif(meanPrePhasingForRead1*100,3))) %>%
           {if(exists('meanPhasingForRead2', where = .)) mutate(.,Phasing_PrephasingRead2=stringr::str_c(signif(meanPhasingForRead2*100,3)," / ",signif(meanPrePhasingForRead2*100,3))) else .} %>%  
    dplyr::select(Lane,meanClusterDensity,percent_PF_Clusters,
                          starts_with("Phasing_Prephasing"),
                          NumberOfReads,NumberOfReadsPF)
  
  
  qMet <- interOpsMetrics$qMet %>% 
    group_by(Lane,Cycle) %>% 
    dplyr::select(Lane,Cycle,starts_with("Q")) %>% 
    summarise_all(sum) %>% 
    ungroup() %>% 
    mutate(q30=rowSums(.[32:52])/rowSums(.[3:52])) %>% 
    mutate(Read = cut(Cycle,
                      breaks=unique(cumsum(c(0,
                                             readPattern(bcl2fastqparams)                             
                      ))),
                      labels= names(readPattern(bcl2fastqparams))[readPattern(bcl2fastqparams) !=0]
    )) %>% 
    group_by(Lane,Read) %>% 
    summarise_all(mean) %>% 
    dplyr::select(Lane,Read,q30) %>% 
    mutate(q30=signif(q30*100,3)) %>%
    mutate(Read=paste0("Q30-",Read)) %>% 
    spread(Read,q30)
  
  fullSummaryDual <- full_join(tileMet,qMet,by = "Lane")
  if(verbose){message("done",appendLF = TRUE)}
  toReport <- list(machineReport=as.data.frame(runPReport),
                   sequencingReport=as.data.frame(fullSummaryDual),
                   demuxReport=as.data.frame(dmxReport))
  return(toReport)
}

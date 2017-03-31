#' The Parameters for BCL2FastQparams object.
#'
#' Parameter class and accessors for use with basecallQC
#'
#' @aliases BCL2FastQparams BCL2FastQparams-BCL2FastQparams
#'
#' @rdname BCL2FastQparams
#' @docType class
#' @export
#'
setClass("BCL2FastQparams", representation(RunDir="character",OutDir="character",RunParameters = "list"))

#' The basecallQC object and constructor.
#'
#' Object and method to handle Illumina basecalling/demultiplexing
#' inputs and output files.
#' Provides sample sheet cleanup, basecall command and
#' summary QC statistics for basecalling/demultiplexing.
#'
#' @aliases basecallQC basecallQC-basecallQC
#'
#' @rdname basecallQC
#' @docType class
#' @export

setClass("basecallQC", representation(BCL2FastQparams="BCL2FastQparams",RunMetadata = "data.frame",
                                      cleanedSampleSheet="data.frame",BCLCommand="character",baseMasks= "data.frame",
                                      baseCallMetrics="list",demultiplexMetrics="list",fqQCmetrics="list"))


#' Set Parameters for BCL2FastQparameters object.
#'
#' Parameter class and accessors
#'
#' @aliases BCL2FastQparams BCL2FastQparams-BCL2FastQparams
#'
#' @rdname BCL2FastQparams
#' @docType methods
#' @param runXML file path to runParameters.xml ,if not specified
#' looks in top level of run directory.
#' @param config file path to config.ini ,if not specified
#' looks in top level of run directory.
#' @param runDir file path to run directory.
#' @param outDir file path to out directory.
#' @param verbose TRUE or FALSE. Messages on or off. Warnings/errors persist
#' @details The BCL2FastQparams object contains slots RunDir, OutDir and RunParameters
#' \itemize{
#' \item{"RunDir"}{ Character string specifying the top level Run directory}
#' \item{"OutDir"}{ Character string specifying the output directory}
#' \item{"RunParameters"}{ A data.frame containing the information from runParameters.xml (See vignette for more details).
#' }
#' }
#' @return A BCL2FastQparams object (See details).
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' BCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#' @export
BCL2FastQparams <- function(runXML=NULL,config=NULL,runDir=NULL,outDir=NULL,verbose=TRUE){
  if(is.null(runDir)) runDir <- getwd(); if(verbose) message("No runDir specified, run directory set to working directory");
  if(is.null(runXML)){
    if(verbose) message("No location for runParameters.xml specified")
    runParameters <- file.path(runDir,"runParameters.xml")
    if(!file.exists(runParameters)) stop("No runParameters.xml found in run directory")
  }
  if(is.null(config)){
    if(verbose) message("No location for config.ini specified")
    config <- file.path(runDir,"config.ini")
    if(!file.exists(config)) stop("No config.ini found in run directory")
  }

  runParameters = runParams(runXML,config)

  if(is.null(outDir)){
    if(verbose) message("No location for outDir specified")
    outDir <- file.path(runDir,runParameters$runParams$Barcode)
    message("outDir set to",outDir)
  }
  new("BCL2FastQparams",
      RunDir=runDir,
      OutDir=outDir,
      RunParameters=runParameters
  )

}
#' The basecallQC object and constructor.
#'
#'
#' @name basecallQC
#' @rdname basecallQC
#' @param bcl2fastqparams A BCL2FastQparams object as created by BCL2FastQparams() constructor.
#' @param RunMetaData Any run metadata to attach (data.frame)
#' @param sampleSheet A sample sheet for Illumina basecalling using bcl2Fastq (See vignette for more details).
#' @param doFQMetric TRUE or FALSE. Perform ShortRead FastQ quality assessment
#' using ShortRead's qa and report function
#' @return basecallQC a basecallQC object (See details for more information)
#' @details The basecallQC object contains slots BCL2FastQparams,
#'  cleanedSampleSheet, baseMasks, BCLCommand, baseCallMetrics, demultiplexMetrics and fqQCmetrics.
#' \itemize{
#' \item{"BCL2FastQparams"}{ A BCL2FastQparams object}
#' \item{"cleanedSampleSheet"}{ A data.frame containing the cleaned sample sheet for
#'  Illumina basecalling using bcl2Fastq versions >= 2.1.7
#' }
#' \item{"baseMasks"}{ A data.frame containing basecall masks per lane for use with bcl2Fastq versions >= 2.1.7. Basemasks in data.frame for reads and indexes as well as the total basemasks for each lane.
#' }
#' \item{"BCLCommand"}{ A character string containing the command to be used for basecalling using bcl2Fastq (versions >= 2.1.7).
#' }
#' \item{"baseCallMetrics"}{ A list containing the full basecalling metrics from ConversionStats.xml. Contains an unsummarised data.frame and basecalling metrics summarised to Sample, Lane, Sample by lane, and Sample by Lane and Tile
#' }
#' \item{"demultiplexMetrics"}{ A list containing the full demultiplexing metrics from DemultiplexingStats.xml. Contains an unsummarised data.frame and demultiplexing metrics filtered to per Sample metrics
#' }
#' \item{"fqQCmetrics"}{ A list containing a data.frame of read counts and links to ShortRead QA reports and a ShortRead QA object containing quality information for generated fastQs.
#' }
#' }
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' @export
basecallQC <- function(bcl2fastqparams,RunMetaData=NULL,sampleSheet=NULL,doFQMetric=FALSE){

  cleanedSampleSheet <- validateBCLSheet(sampleSheet,bcl2fastqparams)
  baseMasks <- createBasemasks(cleanedSampleSheet,bcl2fastqparams)
  toSubmit <- createBCLcommand(bcl2fastqparams,cleanedSampleSheet,baseMasks)
  basecallmetrics <- baseCallMetrics(bcl2fastqparams)
  demultiplexmetrics <- demultiplexMetrics(bcl2fastqparams)
  fastqs <- dir(bcl2fastqparams@OutDir,pattern="*.fastq.gz",full.names=TRUE)
  if(doFQMetric==TRUE & length(fastqs) > 0){
    fqQCmetrics <- qcShortRead(fastqs)
  }else{
    fqQCmetrics <- list(FQQC_Table = NULL,ShortReadQC=NULL)
  }
  basecallQC <- new("basecallQC",
                    BCL2FastQparams = bcl2fastqparams,
                    cleanedSampleSheet = cleanedSampleSheet,
                    baseMasks = baseMasks,
                    BCLCommand=toSubmit,
                    baseCallMetrics = basecallmetrics,
                    demultiplexMetrics = demultiplexmetrics,
                    fqQCmetrics=fqQCmetrics)
  return(basecallQC)
}

runParams <- function(runXML=NULL,config=NULL){
  runParams <- runParameters(runXML)
  configParams <- configParams(config)
  return(list(runParams=runParams,configParams=configParams))
}

#' Gather basecalling metrics from a Run (using Run's ConversionStats.xml file).
#'
#' @name baseCallMetrics
#' @rdname baseCallMetrics
#' @param bcl2fastqparams A BCL2FastQparams object as created by BCL2FastQparams() constructor.
#' @return A list of length two containing the full basecalling metrics from a Run (using Run's ConversionStats.xml file). Contains an unsummarised data.frame and basecalling metrics summarised to Sample, Lane, Sample by lane, and Sample by Lane and Tile.
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' convMetrics <- baseCallMetrics(bcl2fastqparams)
#' @export
baseCallMetrics <- function(bcl2fastqparams){
  convStatsXML <- file.path(bcl2fastqparams@OutDir,"Stats","ConversionStats.xml")
  if(!file.exists(convStatsXML)) return(list(convStatsProcessed=NULL,summarisedConvStats=NULL))
  convStatsProcessed <- processConvStats(convStatsXML)
  summarisedConvStats <- summariseConvStats(convStatsProcessed)
  return(list(convStatsProcessed=convStatsProcessed,
                       summarisedConvStats=summarisedConvStats))
}

#' Gather demultiplexing metrics from a Run (using Run's DemultiplexingStats.xml file).
#'
#' @name demultiplexMetrics
#' @rdname demultiplexMetrics
#' @param bcl2fastqparams A BCL2FastQparams object as created by BCL2FastQparams() constructor.
#' @return A list of length two containing the full demultiplexing metrics from a Run (using Run's DemultiplexingStats.xml file). Contains an unsummarised data.frame and demultiplexing metrics filtered to per Sample metrics 
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' demuxMetrics <- demultiplexMetrics(bcl2fastqparams)
#' @export
demultiplexMetrics <- function(bcl2fastqparams){
  demuxStatsXML <- file.path(bcl2fastqparams@OutDir,"Stats","DemultiplexingStats.xml")
  if(!file.exists(demuxStatsXML)) return(list(demuxStatsProcessed=NULL,summarisedDemuxStats=NULL))
  demuxStatsProcessed <- processDemultiplex(demuxStatsXML)
  summarisedDemuxStats <- summariseDemuxStats(demuxStatsProcessed)
  return(list(demuxStatsProcessed=demuxStatsProcessed,
                       summarisedDemuxStats=summarisedDemuxStats))
}







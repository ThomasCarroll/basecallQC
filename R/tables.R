#' Generate an HTML table of per sample summary demultiplexing statistics
#'
#' @usage
#' \S4method{summaryDemuxTable}{baseCallQC}(object)
#'
#' @docType methods
#' @name summaryDemuxTable
#' @rdname summaryDemuxTable
#' @aliases summaryDemuxTable summaryDemuxTable,baseCallQC-method
#' @author Thomas Carroll
#'
#' @param object A  basecallQC object or list from call to demultiplexMetrics()
#' @param output Whether the report contains frozen or sortable tables. Options are "static" and "html"
#' @return  An HTML table for reporting demultiplexing results.
#' @import stringr XML RColorBrewer methods raster
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' summaryDemuxTable(bclQC,output="static")
#' @export
summaryDemuxTable.basecallQC <- function(object,output="static"){

  toTable <- object@demultiplexMetrics$summarisedDemuxStats$Summary
  if(output=="static"){
    return(kable(toTable))
  }
  if(output=="html"){
    return(DT::datatable(toTable))
  }
}


setGeneric("summaryDemuxTable", function(object="basecallQC",output="character") standardGeneric("summaryDemuxTable"))

#' @rdname summaryDemuxTable
#' @export
setMethod("summaryDemuxTable", signature(object="basecallQC"), summaryDemuxTable.basecallQC)


summaryDemuxTable.list <- function(object,output="static"){
  toTable <- object$summarisedDemuxStats$Summary
  if(output=="static"){
    return(kable(toTable))
  }
  if(output=="html"){
    return(DT::datatable(toTable))
  }
}

#' @rdname summaryDemuxTable
#' @export
setMethod("summaryDemuxTable", signature(object="list"),summaryDemuxTable.list)





#' Creates an HTML table of per sample summary statistics from basecalling results
#'
#' @usage
#' \S4method{summaryConvStatsTable}{baseCallQC}(object)
#' 
#' @docType methods
#' @name summaryConvStatsTable
#' @rdname summaryConvStatsTable
#' @aliases summaryConvStatsTable summaryConvStatsTable,baseCallQC-method
#' @author Thomas Carroll
#'
#' @param object A  basecallQC object or list from call to baseCallMetrics()
#' @param output Whether the report contains frozen or sortable tables. Options are "static" and "html"
#' @return An HTML table for reporting basecalling results.
#' @import stringr XML RColorBrewer methods raster
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' summaryDemuxTable(bclQC,output="static")
#' @export
summaryConvStatsTable.basecallQC <- function(object,output="static"){

  toTable <- object@baseCallMetrics$summarisedConvStats$Sample_Stats
  if(output=="static"){
    return(kable(toTable))
  }
  if(output=="html"){
    return(DT::datatable(toTable))
  }
}

setGeneric("summaryConvStatsTable", function(object="basecallQC",output="character") standardGeneric("summaryConvStatsTable"))

#' @rdname summaryConvStatsTable
#' @export
setMethod("summaryConvStatsTable", signature(object="basecallQC"), summaryConvStatsTable.basecallQC)


summaryConvStatsTable.list <- function(object,output="static"){
  
  toTable <- object$summarisedConvStats$Sample_Stats
  if(output=="static"){
    return(kable(toTable))
  }
  if(output=="html"){
    return(DT::datatable(toTable))
  }
}

#' @rdname summaryConvStatsTable
#' @export
setMethod("summaryConvStatsTable", signature(object="list"),summaryConvStatsTable.list)


#' Generate an HTML table linking to per sample summary fastq QC statistics from ShortRead
#'
#' Creates an HTML table linking to per sample summary fastq QC statistics from ShortRead
#'
#'
#' @docType methods
#' @name makeFQTable
#' @rdname makeFQTable
#'
#' @author Thomas Carroll
#'
#' @param object A basecall QC object as returned from basecallQC function
#' @param output Whether the report contains frozen or sortable tables. Options are "static" and "html"
#' @return A HTML table for reporting fastq QC results from ShortRead. 
#' Table contains read counts and links to ShortRead QA reports per sample.
#' @import stringr XML RColorBrewer methods raster
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' #makeFQTable(bclQC,output="static")
#' @export

makeFQTable <- function(object,output="static"){
  fqQCTable <- object@fqQCmetrics$FQQC_Table
    if(!is.null(fqQCTable)){
      if(output=="static"){
        table <- kable(fqQCTable,escape = FALSE)
      }
      if(output=="html"){
        table <- DT::datatable(fqQCTable,escape=FALSE)
      }
      return(table)
    }else{
      return(NULL)
  }
}

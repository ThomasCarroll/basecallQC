
#' Barplot of Illumina basecalling statistics for reads passing filter.
#'
#' Produces a barplot of basecalling statistics for reads passing filter.
#'
#' @usage
#' \S4method{passFilterBar}{baseCallQC}(object,groupBy,metricToPlot)
#'
#' @docType methods
#' @name passFilterBar
#' @rdname passFilterBar
#' @aliases passFilterBar passFilterBar,baseCallQC-method
#' 
#' @author Thomas Carroll and Marian Dore
#' @param object baseCallQC A  basecallQC object
#' @param groupBy Character vector defining how plot will be grouped
#' @param metricToPlot Character vector defining which metric will be displayed in plot.
#' @return ggPlotObject A ggplot2 object.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' plot <- passFilterBar(bclQC)

#' @export
passFilterBar.basecallQC <- function(object,groupBy=c("Lane"),metricToPlot="Yield"){
  groupByS <- unique(c(groupBy,"Filter"))
  groupByG <- unique(c(groupBy))
  toPlot <- object@baseCallMetrics$convStatsProcessed %>%
    group_by_(.dots=as.list(groupByS)) %>%
    filter(Sample != "all") %>%
    summarise_(.dots = setNames(list(interp( ~sum(as.numeric(var)),
                                             var=as.name(metricToPlot))
    )
    ,metricToPlot)) %>%
    spread_("Filter",metricToPlot) %>%
    mutate(Ff=Raw-Pf) %>%
    dplyr::select(-Raw) %>%
    tbl_df %>%
    gather_(key_col="PassFilter",value_col=as.name(metricToPlot),c("Ff","Pf"))
  p <- ggplot(data=toPlot,aes_string(x=groupByG,y=metricToPlot,fill="PassFilter"))+geom_bar(stat = "identity")+ coord_flip()
  return(p)
}

setGeneric("passFilterBar", function(object="basecallQC",groupBy="character",metricToPlot="character") standardGeneric("passFilterBar"))

#' @rdname passFilterBar
#' @export
setMethod("passFilterBar", signature(object="basecallQC"), passFilterBar.basecallQC)


#' Boxplot of Illumina basecalling statistics for reads passing filter.
#'
#' Produces a boxplot of basecalling statistics for reads passing filter.
#'
#' @usage
#' \S4method{passFilterBoxplot}{baseCallQC}(object,groupBy,metricToPlot)
#' @docType methods
#' @name passFilterBoxplot
#' @rdname passFilterBoxplot
#' @aliases passFilterBoxplot passFilterBoxplot,baseCallQC-method
#' 
#' @author Thomas Carroll and Marian Dore
#' @param object basecallQC A  basecall QC object
#' @param groupBy Character vector of how data is grouped for plotting.
#' @param metricToPlot Character vector defining which metric will be displayed in plot.
#' @return ggPlotObject A ggplot2 object.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' plot <- passFilterBoxplot(bclQC,groupBy = "Sample")
#' @export
passFilterBoxplot.basecallQC <- function(object,groupBy=c("Lane"),metricToPlot="Yield"){
  groupByS <- unique(c("Lane","Sample","Tile","Filter"))
  groupByG <- unique(c(groupBy))
  toPlot <- object@baseCallMetrics$convStatsProcessed %>%
    group_by_(.dots=as.list(groupByS)) %>%
    filter(Sample != "all") %>%
    summarise_(.dots = setNames(list(interp( ~sum(as.numeric(var)),
                                             var=as.name(metricToPlot))
    )
    ,metricToPlot)) %>%
    spread_("Filter",metricToPlot) %>%
    mutate(Ff=Raw-Pf) %>%
    dplyr::select(-Raw) %>%
    tbl_df %>%
    gather_(key_col="PassFilter",value_col=as.name(metricToPlot),c("Ff","Pf"))
  p <- ggplot(data=toPlot,aes_string(x=groupByG,y="Yield",fill="PassFilter"))+geom_violin(scale="width")+ coord_flip()+facet_grid(PassFilter~.,scales = "free")+geom_jitter()
  return(p)
}

setGeneric("passFilterBoxplot", function(object="basecallQC",groupBy="character",metricToPlot="character") standardGeneric("passFilterBoxplot"))

#' @rdname passFilterBoxplot
#' @export
setMethod("passFilterBoxplot", signature(object="basecallQC"), passFilterBoxplot.basecallQC)


#' Tile plot of Illumina basecalling statistics for reads passing filter.
#'
#' Produces a plot of metrics per Tile for basecalling statistics of reads passing/failing filter.
#'
#' @usage
#' \S4method{passFilterTilePlot}{baseCallQC}(object,metricToPlot)
#' @docType methods
#' @name passFilterTilePlot
#' @rdname passFilterTilePlot
#' @aliases passFilterTilePlot passFilterTilePlot,baseCallQC-method
#' @author Thomas Carroll and Marian Dore
#' @param object baseCallQC A  basecall QC object
#' @param metricToPlot Character vector defining which metric will be displayed in plot.
#' @return ggPlotObject A ggplot2 object.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' plot <- passFilterTilePlot(bclQC,metricToPlot="Yield")
#' @export
passFilterTilePlot.basecallQC <- function(object,metricToPlot="Yield"){
  groupByS <- unique(c("Lane","Sample","Tile","Filter"))
  #groupByG <- unique(c(groupBy))
  toPlot <- object@baseCallMetrics$convStatsProcessed %>%
    group_by_(.dots=as.list(groupByS)) %>%
    filter(Sample != "all") %>%
    summarise_(.dots = setNames(list(interp( ~sum(as.numeric(var)),
                                             var=as.name(metricToPlot))
    )
    ,metricToPlot)) %>%
    spread_("Filter",metricToPlot) %>%
    mutate(Ff=Raw-Pf) %>%
    dplyr::select(-Raw) %>%
    tbl_df %>%
    gather_(key_col="PassFilter",value_col=as.name(metricToPlot),c("Ff","Pf")) %>%
    mutate(Surface=str_sub(Tile,1,1),Box=str_sub(Tile,2,2),Pos=str_sub(Tile,3))
    pPf <- filter(toPlot,PassFilter=="Pf") %>%
    ggplot(aes(x=Box,y=Pos))+geom_tile(aes_string(fill=metricToPlot))+facet_grid(~Lane)+scale_fill_gradient2(low = "white", high = "darkblue")+theme_bw()
    pFf <- filter(toPlot,PassFilter=="Ff") %>%
    ggplot(aes(x=Box,y=Pos))+geom_tile(aes_string(fill=metricToPlot))+facet_grid(~Lane)+scale_fill_gradient2(low = "white", high = "darkblue")+theme_bw()
  return(list(PassFilter=pPf,FailFilter=pFf))
}

setGeneric("passFilterTilePlot", function(object="basecallQC",metricToPlot="character") standardGeneric("passFilterTilePlot"))

#' @rdname passFilterTilePlot
#' @export
setMethod("passFilterTilePlot", signature(object="basecallQC"), passFilterTilePlot.basecallQC)


#' Boxplot of Illumina demultiplexing statistics.
#'
#' Produces a boxplot of for demultiplexing statistics of reads with perfect/mismatched barcode.
#'
#' @usage
#' \S4method{demuxBoxplot}{baseCallQC}(object,groupBy)
#' @docType methods
#' @name demuxBoxplot
#' @rdname demuxBoxplot
#' @aliases demuxBoxplot demuxBoxplot,baseCallQC-method
#' @author Thomas Carroll and Marian Dore
#' @param object baseCallQC A  basecall QC object
#' @param groupBy Character vector of lane and/or Sample
#' @return ggPlotObject A ggplot2 object.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' plot <- demuxBoxplot(bclQC)
#' @export
demuxBoxplot.basecallQC <- function(object,groupBy=c("Lane")){
  metricToPlot <- "Count"
  groupByS <- unique(c("Lane","Sample","Project","Barcode","BarcodeStat"))
  groupByG <- unique(c(groupBy))

  toPlot <- object@demultiplexMetrics$demuxStatsProcessed %>%
    group_by_(.dots=as.list(groupByS)) %>%
    filter(Sample != "all") %>%
    summarise_(.dots = setNames(list(interp( ~sum(as.numeric(var)),
                                             var=as.name(metricToPlot))
    )
    ,metricToPlot)) %>%
    spread_("BarcodeStat",metricToPlot) %>%
    mutate(mismatchedBarcodeCount=BarcodeCount-PerfectBarcodeCount) %>%
    dplyr::select(-BarcodeCount) %>%
    tbl_df %>%
    gather_(key_col="BarcodeCount",value_col=as.name(metricToPlot),c("mismatchedBarcodeCount","PerfectBarcodeCount"))
    p <- ggplot(data=toPlot,aes_string(x=groupByG,y=metricToPlot,fill="BarcodeCount"))+geom_violin(scale = "width")+ coord_flip()+facet_grid(BarcodeCount~.)
  return(p)
}

setGeneric("demuxBoxplot", function(object="basecallQC",groupBy="character") standardGeneric("demuxBoxplot"))

#' @rdname demuxBoxplot
#' @export
setMethod("demuxBoxplot", signature(object="basecallQC"), demuxBoxplot.basecallQC)

#' Barplot of Illumina demultiplexing statistics.
#'
#' Produces a barplot of for demultiplexing statistics of reads with perfect/mismatched barcode.
#'
#' @usage
#' \S4method{demuxBarplot}{baseCallQC}(object,groupBy)
#' @docType methods
#' @name demuxBarplot
#' @rdname demuxBarplot
#' @aliases demuxBarplot demuxBarplot,baseCallQC-method
#' @author Thomas Carroll and Marian Dore
#' @param object baseCallQC A  basecall QC object
#' @param groupBy Character vector of lane and/or Sample
#' @return ggPlotObject A ggplot2 object.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' plot <- demuxBarplot(bclQC)
#' @export
demuxBarplot.basecallQC <- function(object,groupBy=c("Lane")){
  metricToPlot <- "Count"
  groupByS <- unique(c("Lane","Sample","Project","Barcode","BarcodeStat"))
  groupByG <- unique(c(groupBy))

  toPlot <- object@demultiplexMetrics$demuxStatsProcessed %>%
    group_by_(.dots=as.list(groupByS)) %>%
    filter(Sample != "all") %>%
    summarise_(.dots = setNames(list(interp( ~sum(as.numeric(var)),
                                             var=as.name(metricToPlot))
    )
    ,metricToPlot)) %>%
    spread_("BarcodeStat",metricToPlot) %>%
    mutate(mismatchedBarcodeCount=BarcodeCount-PerfectBarcodeCount) %>%
    dplyr::select(-BarcodeCount) %>%
    tbl_df %>%
    gather_(key_col="BarcodeCount",value_col=as.name(metricToPlot),c("mismatchedBarcodeCount","PerfectBarcodeCount"))
  p <- ggplot(data=toPlot,aes_string(x=groupByG,y=metricToPlot,fill=groupByG))+geom_bar(stat="identity")+ coord_flip()+facet_grid(BarcodeCount~.)
  return(p)
}

setGeneric("demuxBarplot", function(object="basecallQC",groupBy="character") standardGeneric("demuxBarplot"))

#' @rdname demuxBarplot
#' @export
setMethod("demuxBarplot", signature(object="basecallQC"), demuxBarplot.basecallQC)

demuxBarplot.list <- function(object,groupBy=c("Lane")){
  metricToPlot <- "Count"
  groupByS <- unique(c("Lane","Sample","Project","Barcode","BarcodeStat"))
  groupByG <- unique(c(groupBy))
  
  toPlot <- object$demuxStatsProcessed %>%
    group_by_(.dots=as.list(groupByS)) %>%
    filter(Sample != "all") %>%
    summarise_(.dots = setNames(list(interp( ~sum(as.numeric(var)),
                                             var=as.name(metricToPlot))),
                                metricToPlot)) %>%
    spread_("BarcodeStat",metricToPlot) %>%
    mutate(mismatchedBarcodeCount=BarcodeCount-PerfectBarcodeCount) %>%
    dplyr::select(-BarcodeCount) %>%
    tbl_df %>%
    gather_(key_col="BarcodeCount",value_col=as.name(metricToPlot),c("mismatchedBarcodeCount","PerfectBarcodeCount"))
  p <- ggplot(data=toPlot,aes_string(x=groupByG,y=metricToPlot,fill=groupByG))+geom_bar(stat="identity")+ coord_flip()+facet_grid(BarcodeCount~.)
  return(p)
}

#' @rdname demuxBarplot
#' @export
setMethod("demuxBarplot", signature(object="list"), demuxBarplot.list)


  
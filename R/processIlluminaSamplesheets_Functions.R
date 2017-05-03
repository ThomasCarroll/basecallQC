if(getRversion() >= "2.15.1")  utils::globalVariables(c("BarcodeCount","BarcodeStat","Box","Count","ID",
                                                        "IndexRead1","IndexRead2","Name","PassFilter",
                                                        "PerfectBarcodeCount","Pf","Pos","Project",
                                                        "ApplicationName","ApplicationVersion","Barcode","ChemistryVersion",
                                                        "ClusterCount","ClusterDensity","ComputerName","Cycle","ExperimentName",
                                                        "Flowcell","NumberOfClusters","NumberOfClustersPF","NumberOfReads",
                                                        "NumberOfReadsPF","PhasingForRead1","PhasingForRead2","PrePhasingForRead1",
                                                        "PrePhasingForRead2","RTAVersion","Read","Reads","RunID","RunMode","Yield30",
                                                        "Yield30Sum","YieldSum","code","distinct","full_join","meanClusterDensity",
                                                        "meanPhasingForRead1","meanPhasingForRead2","meanPrePhasingForRead1",
                                                        "meanPrePhasingForRead2","percent_PF_Clusters","percent_PF_ClustersSD","q30",
                                                        "sd","sdClusterDensity","spread","starts_with","str_replace_all","summarise_all",
                                                        "ungroup","value",
                                                        "Raw","Read1","Read2","Sample","SampleID",
                                                        "SampleName","Sample_ID","Sample_Name",
                                                        "Sample_Project","Tile","Yield","bcl2fastqparams",
                                                        "index","index2","Index","Index2","Lane","basemask",
                                                        "index1Mask","index2Mask","indexLength","indexLength2",
                                                        "read1Mask","read2Mask","SampleRef","."))



#' Illumina sample sheet cleaning and updating for
#' bcl2Fastq versions >= 2.1.7
#'
#' Parses an Illumina bcl2Fastq sample sheet  to create a
#' standardised/updated sample sheet for bcl2Fastq >= Version 2.1.7
#'
#'
#' @docType methods
#' @name validateBCLSheet
#' @rdname validateBCLSheet
#'
#' @author Thomas Carroll and Marian Dore
#' @param sampleSheet File location of a sample sheet for Illumina basecalling using bcl2Fastq (See vignette for more details).
#' @param param A BCL2FastQparams object
#' @return cleanedSampleSheet A data.frame containing the cleaned sample sheet for
#'  Illumina basecalling using bcl2Fastq versions >= 2.1.7.
#' @import stringr XML RColorBrewer methods raster BiocStyle lazyeval
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#' cleanedSampleSheet <- validateBCLSheet(sampleSheet,param=bcl2fastqparams)
#'
#' @export
validateBCLSheet <- function(sampleSheet,param=bcl2fastqparams){
  #runParam <- runParameters(param)
  fread(sampleSheet,sep=",",header=TRUE,stringsAsFactors=FALSE,skip="Sample") %>%
    tbl_df %>%
    {if(exists('Project', where = .) & !exists('Sample_Project', where = .)) dplyr::rename(.,Sample_Project = Project) else .} %>%
    {if(exists('SampleID', where = .) & !exists('Sample_ID', where = .)) dplyr::rename(.,Sample_ID = SampleID) else .} %>%
    {if(exists('ID', where = .) & !exists('Sample_ID', where = .)) dplyr::rename(.,Sample_ID = ID) else .} %>%
    {if(exists('SampleName', where = .) & !exists('Sample_Name', where = .)) dplyr::rename(.,Sample_Name = SampleName) else .} %>%
    {if(exists('Name', where = .) & !exists('Sample_Name', where = .)) dplyr::rename(.,Sample_Name = Name) else .} %>%
    {if(exists('SampleRef', where = .) & !exists('Sample_Name', where = .)) dplyr::rename(.,Sample_Name = SampleRef) else .} %>%
    {if(exists('index', where = .) & !exists('Index', where = .)) dplyr::rename(.,Index = index) else .} %>%
    {if(exists('index2', where = .) & !exists('Index2', where = .)) dplyr::rename(.,Index2 = index2) else .} %>%
    {if(!exists('Index2', where = .)) tidyr::separate(.,Index, c("Index", "Index2"), "-",fill="right") else .} %>%
    mutate(Sample_Project = if (exists('Sample_Project', where = .)) Sample_Project else NA,
           Lane = if (exists('Lane', where = .)) Lane else NA,
           Sample_ID = if (exists('Sample_ID', where = .)) Sample_ID else "",
           Sample_Name = if (exists('Sample_Name', where = .)) Sample_Name else "",
           Index = if (exists('Index', where = .)) Index else NA,
           Index2 = if (exists('Index2', where = .)) Index2 else NA) %>%
    tbl_df %>%
    dplyr::select(Sample_Project,Lane,Sample_ID,Sample_Name,Index,Index2,everything()) %>%
    mutate(Sample_ID=replace(as.character(Sample_ID),grep("^\\d",Sample_ID),paste0("Sample_",Sample_ID[grep("^\\d",Sample_ID)]))) %>% 
    mutate(Sample_Name=replace(as.character(Sample_Name),which(Sample_Name == "" | is.na(Sample_Name)),Sample_ID[which(Sample_Name == "" | is.na(Sample_Name))])) %>% 
    mutate(Sample_Name=replace(as.character(Sample_Name),grep("^\\d",Sample_Name),paste0("Sample_",Sample_Name[grep("^\\d",Sample_Name)]))) %>% 
    mutate(Sample_ID=validNames(Sample_ID,prefix="Sample_")) %>%
    mutate(Sample_Name=validNames(Sample_Name,prefix="Sample_")) %>%
    mutate(Index2=replace(as.character(Index2),is.na(Index2),"")) %>% 
    mutate(Index=str_to_upper(Index)) %>%
    mutate(Index2=str_to_upper(Index2)) %>%
    mutate(Sample_ID=gsub("\\?|\\(|\\)|\\[|\\]|\\\\|/|\\=|\\+|<|>|\\:|\\;|\"|\'|\\*|\\^|\\||\\&|\\.","_",Sample_ID)) %>%
    mutate(Sample_Name=gsub("\\?|\\(|\\)|\\[|\\]|\\\\|/|\\=|\\+|<|>|\\:|\\;|\"|\'|\\*|\\^|\\||\\&|\\.","_",Sample_Name)) %>%
    mutate(Index=str_trim(Index,"both"),
           Index2=str_trim(Index2,"both"))    %>%
    mutate(Index=str_sub(Index,1,as.numeric(indexlengths(param)$IndexRead1)),    
           Index2=str_sub(Index2,1,as.numeric(indexlengths(param)$IndexRead2)))
}

#' Function to create basemasks for basecalling from Illumina samplesheet (for bcl2Fastq versions >= 2.1.7).
#'
#' Parses the Illumina sample sheet for versions >= 2.1.7 and creates basemasks.
#'
#'
#' @docType methods
#' @name createBasemasks
#' @rdname createBasemasks
#'
#' @author Thomas Carroll and Marian Dore
#' @param cleanedSampleSheet Data.frame of cleaned samplesheet for Illumina basecalling using bcl2Fastq versions >= 2.1.7 (see vignette for more details)
#' @param param A BCL2FastQparams object
#' @return A data.frame containing basecall masks per lane for reads and indexes as well as per lane complete basemasks.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#'

#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#'
#' cleanedSampleSheet <- validateBCLSheet(sampleSheet,param=bcl2fastqparams)
#' basemasks <- createBasemasks(cleanedSampleSheet,param=bcl2fastqparams)
#'
#' @export
createBasemasks <- function(cleanedSampleSheet,param){
  indexCombinations <- cleanedSampleSheet %>%
    mutate(Index2=ifelse(is.na(Index2), "", Index2),Index=ifelse(is.na(Index), "", Index)) %>%
    mutate(indexLength=str_length(Index),indexLength2=str_length(Index2)) %>%
    dplyr::group_by(Lane) %>% dplyr::count(indexLength,indexLength2)

  if(nrow(indexCombinations) == length(unique(indexCombinations$Lane))){
    baseMasks <- indexCombinations %>%
      mutate(index1Mask = str_c(str_dup("I",indexLength),
                                str_dup("N",indexlengths(param)$IndexRead1-indexLength)),
             index2Mask = str_c(str_dup("I",indexLength2),
                                str_dup("N",indexlengths(param)$IndexRead2-indexLength2))) %>%
      mutate(read1Mask = str_c(str_dup("Y",as.numeric(readlengths(param)$Read1))),
             read2Mask = str_c(str_dup("Y",as.numeric(readlengths(param)$Read2)))) %>%
      mutate(read1Mask = str_replace(read1Mask,"Y$","N"),
             read2Mask = str_replace(read2Mask,"Y$","N")) %>%
      mutate(basemask = str_c(read1Mask,index1Mask,index2Mask,read2Mask,sep=",")) %>%
      mutate(basemask = str_c(Lane,":",basemask)) %>%
      mutate(basemask = str_replace(basemask,",,",",")) %>%
      mutate(basemask = str_replace(basemask,",$","")) %>%
      tbl_df %>%
      dplyr::select(Lane,basemask,read1Mask,index1Mask,index2Mask,read2Mask)
      }
}

#' Function to create command for Illumina basecalling/demultiplexing using bcl2fastq versions >= 2.1.7.
#'
#' Creates the command to be used for basecalling/demultiplexing with bcl2fastq versions >= 2.1.7
#'
#'
#' @docType methods
#' @name createBCLcommand
#' @rdname createBCLcommand
#'
#' @author Thomas Carroll and Marian Dore
#' @param bcl2fastqparams A BCL2FastQparams object.
#' @param cleanedSampleSheet Data.frame of cleaned samplesheet for Illumina basecalling/demultiplexing using bcl2fastq versions >= 2.1.7 (see vignette for more details)
#' @param baseMasks A data.frame of basemasks as created by createBasemasks() function
#' @return A character vector containing the command for Illumina basecalling using bcl2fastq versions >= 2.1.7
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#'

#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' bcl2fastqparams <- BCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#'
#' cleanedSampleSheet <- validateBCLSheet(sampleSheet,param=bcl2fastqparams)
#' baseMasks <- createBasemasks(cleanedSampleSheet,param=bcl2fastqparams)
#' toSubmit <- createBCLcommand(bcl2fastqparams,cleanedSampleSheet,baseMasks)
#' @export
createBCLcommand <- function(bcl2fastqparams,cleanedSampleSheet,baseMasks){
  sampleSheetLocation <- paste0(file.path(bcl2fastqparams@RunDir,bcl2fastqparams@RunParameters$runParams$Barcode),".csv")
  bclPath <- bcl2fastqparams@RunParameters$configParams[bcl2fastqparams@RunParameters$configParams$name == "configureBclToFastq","value"]
  write.table("[DATA]",file=sampleSheetLocation,sep="",quote=FALSE,row.names=FALSE)
  write.table(cleanedSampleSheet,file=sampleSheetLocation,sep=",",quote=FALSE,row.names=FALSE,append = TRUE)
  baseMasksToUse <- str_c("--use-bases-mask ",dplyr::select(tbl_df(baseMasks),basemask)$basemask,collapse = " ")
  bclcommand <- str_c(as.vector(bclPath$value),"--output-dir ",bcl2fastqparams@OutDir,"--sample-sheet",sampleSheetLocation,baseMasksToUse,sep=" ")
  return(bclcommand)
}

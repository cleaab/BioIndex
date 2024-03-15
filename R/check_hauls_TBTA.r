
#' Check presence in TB of the hauls in TA
#' @description
#' The function check the presence of the TB (catch data table) hauls in the TA (haul data table)
#'
#' @param DataTA data frame of TA table
#' @param DataTB data frame of TB table
#' @param year reference year for the analysis
#' @param wd working directory
#' @param suffix name of the logfile
#' @param verbose boolean. If TRUE messages are promted in the console
#' @export
check_hauls_TBTA <- function(DataTA,DataTB,year,wd=NA,suffix, verbose=FALSE){

  if (FALSE){
    #library(MEDITS)
    wd <- tempdir()
    suffix=paste(as.character(Sys.Date()),format(Sys.time(), "_time_h%Hm%Ms%OS0"),sep="")
    # DataTA <- RoME::TA # read.csv("~/GitHub/RoME/data/TA_GSA18_1994-2018.csv", sep=";") # DataTA[DataTA$YEAR == 2018, ]
    # DataTB <- RoME::TB # read.csv("~/GitHub/RoME/data/TB_GSA18_1994-2018.csv", sep=";") # DataTB[DataTB$YEAR == 2018, ]
    year=2008
    # check_hauls_TBTA(DataTA,DataTB,year,wd,suffix)
  }

  if (is.na(wd)) {
    wd = tempdir()
    if (verbose){
      message(paste("Missing working directory. Results are saved in the temporary folder: ",wd,sep=""))
    }
  }

  if (!file.exists(file.path(wd, "Logfiles"))){
    dir.create(file.path(wd, "Logfiles"), recursive = TRUE, showWarnings = FALSE)
  }
  if (!exists("suffix")){
    suffix=paste(as.character(Sys.Date()),format(Sys.time(), "_time_h%Hm%Ms%OS0"),sep="")
  }
  numberError = 0
  Errors <- file.path(wd,"Logfiles",paste("Logfile_",suffix,".dat",sep=""))
  if (!file.exists(Errors)){
    file.create(Errors)
  }

  ### FILTERING DATA FOR THE SELECTED YEAR
  arg <- "year"
  if (!exists(arg)) {
    stop(paste0("'", arg, "' argument should be provided"))
  } else if (length(year) != 1) {
    stop(paste0("only one value should be provided for '", arg, "' argument"))
  } else if (is.na(year)) {
    stop(paste0(arg, " argument should be a numeric value"))
  }
  DataTA <- DataTA[DataTA$YEAR == year, ]
  DataTB <- DataTB[DataTB$YEAR == year, ]
  ########################################

  ResultTB = DataTB
  write(paste("\n----------- check presence in TA of TB hauls - ", ResultTB$YEAR[1]), file = Errors, append = TRUE)
  ResultTA = DataTA


  ResultTB <- unique(ResultTB$HAUL_NUMBER)
if (length(ResultTB) != 0) {
  j <- 1
  for (j in 1:length(ResultTB)) {
    ResultTA_temp <- ResultTA[which(ResultTA$HAUL_NUMBER == ResultTB[j]), ]
    if (nrow(ResultTA_temp) == 0) {
      write(paste("No haul", ResultTB[j], "in TA"), file = Errors, append = TRUE)
      numberError <- numberError + 1
    }
  }
}
  if (numberError ==0) {
    write(paste("No error occurred"), file = Errors, append = TRUE)
  }
   if (file.exists(file.path(tempdir(), "Logfiles"))){
  unlink(file.path(tempdir(),"Logfiles"),recursive=T)
  }
  if (file.exists(file.path(tempdir(), "Graphs"))){
  unlink(file.path(tempdir(),"Graphs"),recursive=T)
    }
	if (file.exists(file.path(tempdir(), "files R-Sufi"))){
  unlink(file.path(tempdir(),"files R-Sufi"),recursive=T)
    }
  if (numberError ==0) {
    return(TRUE)
  } else { return(FALSE) }
}

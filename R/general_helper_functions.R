#' Function to insert latex equations into .md dox
#'
#' Function copied directly from here:
#' https://github.com/STAT545-UBC/Discussion/issues/102
#` This is a workaround for the fact that markdown doesn't read latex so
#' equations don't render properly
#'
#' @param latex latex code to generate equation
#' @export
#'

latexImg = function(latex){

    link = paste0('http://latex.codecogs.com/gif.latex?',
           gsub('\\=','%3D',URLencode(latex)))

    link = gsub("(%..)","\\U\\1",link,perl=TRUE)
    return(paste0('![](',link,')'))
}

#' Function to generate log transformed sequence evenly distributed across broad orders of magnitude
#'
#' Functions just like `seq` but evenly distributes values across the full range
#' rather than for instance `seq(0.000001, 10000, length.out = 100)` returning values that are all >100
#'
#' @param min minimum value in the sequence
#' @param max maximum value in the sequence
#' @param seq.length length of the sequence
#'
#' @return numeric vector spanning min and max with n = seq.length entries
#' @export
#'
#'
exp_seq <- function(min, max, seq.length){
  exp(seq(log(min), log(max), length.out = seq.length))
}

#' Function that takes a google drive file id and loads the file into your local R environment
#'
#' You can find file ids by going to the sharing settings of a particular csv and
#' copy/pasting the alphanumeric string between "https://drive.google.com/file/d/" and "/view?usp=sharing"
#'
#' @param id googledrive id
#'
#' @return R object from csv in googledrive
#' @export

load_csv_from_googledrive <- function(id){
  require(readr)
  require(googledrive)

  temp <- tempfile(fileext = ".csv")
  dl <- drive_download(as_id(id), path = temp, overwrite = TRUE)

  out <- read_csv(temp)

  file.remove(temp)

  return(out)
}

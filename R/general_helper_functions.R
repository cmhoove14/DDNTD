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

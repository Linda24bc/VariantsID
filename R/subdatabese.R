#' Title
#'
#' @param HbDatabase
#' @param Mshift
#' @param error_Da_start
#' @param error_Da_end
#' @import dplyr
#' @return
#' @export

SubDatabase <- function(HbDatabase, Mshift= -29.97, error_Da_start=-0.06, error_Da_end=0.05){
  s1 <- Mshift + error_Da_start
  e1 <- Mshift + error_Da_end
  subHbBeta1 <- dplyr::filter(HbDatabase, Delta.mass <=e1 & Delta.mass >= s1)
  return(subHbBeta1)
}

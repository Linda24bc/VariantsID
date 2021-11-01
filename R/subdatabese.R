#' SubDatabase
#'
#' @param HbDatabase list of original database or updated database
#' @param Mshift observed mass shift at MS1 level
#' @param error_Da_start mass error tolerance in Da
#' @param error_Da_end mass error tolerance in Da
#' @importFrom  dplyr filter


SubDatabase <- function(HbDatabase, Mshift= -29.97, error_Da_start=-0.06, error_Da_end=0.05){
  s1 <- Mshift + error_Da_start
  e1 <- Mshift + error_Da_end
  subHbBeta1 <- dplyr::filter(HbDatabase, Delta.mass <=e1 & Delta.mass >= s1)
  return(subHbBeta1)
}

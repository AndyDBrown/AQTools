#' Return example of receptor dataframe
#'
#' receptor dataframe
#' @param none
#' @return receptro dataframe example
#' @examples
#' receptorExample();
#'
#' @export
#'
receptorExample <- function() {
  return(data.frame(ID = c("R1", "R2", "R4"), RoadNOx = c(12.1, 14.65, 13.4),
                    No2Backgd = c(19, 24.5, 17.84), LAs = c("Armagh Banbridge and Craigavon")))
}

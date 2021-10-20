#' Return explanation for receptor dataframe
#'
#' receptor dataframe explanation
#' @param none
#' @return receptor dataframe explanation
#' @examples
#' receptorExplained();
#'
#' @export
#'

receptorExplained <- function() {
  message("In the example shown, 'ID' column represents a unique ID for each row which represets known locations.
  The 'RoadNOx' column represents the Road NOx concentration at the receptor defined in the 'ID' column.
  The 'NO2Backgd' column represents the background NO2 concentration for the receptor deinfed in the 'ID' column.")
}

#' #' Calculate the density of the estimated locations.
#' #'
#' #' @param fit ACRE fit object
#' #' @param id 
#' #' @param infotypes 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' calculate_density <- function (fit, id, infotypes) {
#'   combined_density <- 1
#'   
#'   # Calculate respective densities
#'   # Note that if an info-type was not used in the fitted model, 
#'   # get_..._density(fit) will return a vector of 0's
#'   location_density <- get_location_density(fit)
#'   capture_density <- get_capture_density(fit)
#'   dist_density <- get_dist_density(fit)
#'   bearing_density <- get_bearing_density(fit)
#'   ss_density <- get_ss_density(fit)
#'   toa_density <- get_toa_density(fit)
#'   
#'   if ("location" %in% infotypes) combined_density *= location_density
#'   if ("capt" %in% infotypes) combined_density *= capture_density
#'   if ("dist" %in% infotypes) combined_density *= dist_density
#'   if ("bearing" %in% infotypes) combined_density *= bearing_density
#'   if ("ss" %in% infotypes) combined_density *= ss_density
#'   if ("toa" %in% infotypes) combined_density *= toa_density
#'   
#'   return(combined_density)
#' }
#' 
#' get_location_density <- function(fit) {
#'   # return fit$output.tmb$report$location_density
#' }
#' 
#' get_capture_density <- function(fit) {
#'   # return fit$output.tmb$report$capture_density
#' }
#' 
#' get_dist_density <- function(fit) {
#'   # return fit$output.tmb$report$dist_density
#' }
#' 
#' get_bearing_density <- function(fit) {
#'   # return fit$output.tmb$report$bearing_density
#' }
#' 
#' get_ss_density <- function(fit) {
#'   # return fit$output.tmb$report$ss_density
#' }
#' 
#' get_toa_density <- function(fit) {
#'   # return fit$output.tmb$report$toa_density
#' }
#' 
#' 
#' 
#' 
#' 
#' 
#' 

#' COVID-19 epidemic curve for the Republic of Ireland
#'
#' Recorded cases of COVID-19 in the Republic of Ireland, ordered by
#' epidemiological date.
#'
#' The epidemiological date is the earliest recorded date
#' associated with a confirmed case of COVID-19. This is either the date of
#' onset of symptoms, date of diagnosis, laboratory specimen collection date,
#' laboratory received date, laboratory reported date or the notification date.
#' Sorting cases by their epidemiological date strips out some random effects on
#' the epidemic curve introduced by reporting delays
#'
#' @format A data frame with 363 observations of 2 variables:
#' \describe{
#'   \item{date}{Epidemiological date.}
#'   \item{count}{Recorded COVID-19 cases}
#' }
#' @source{https://www.cso.ie/en/releasesandpublications/br/b-cdc/covid-19deathsandcasesseries19/}
"covid_incidence_roi_epidemiological_date"

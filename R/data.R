#' Sample Base Dataset
#'
#' A sample matrix containing subject-level demographic information and the number of observed events for each subject.
#'
#' @format A matrix with $n$ rows and the following variables:
#' \describe{
#' \item{patid}{Subject ID.}
#' \item{MaleInd}{Indicator for male sex.}
#' \item{CityAInd}{Indicator for residence in region A.}
#' \item{CityBInd}{Indicator for residence in region B.}
#' \item{nmvst}{Number of observed events.}
#' \item{p1stvt}{Row index of the subject's first observed event in the
#' \code{my.visit} matrix.}
#' \item{birthl}{Lower endpoint of the birthdate interval, represented as a Julian date.}
#' \item{birthr}{Upper endpoint of the birthdate interval, represented as a Julian date.}
#' \item{birthJulian}{Exact birthdate, represented as a Julian date.}
#' }
"my.base"

#' Sample Visit Dataset
#'
#' A sample matrix containing event-level information for each subject.
#'
#' @format A matrix with $n$ rows and the following variables:
#' \describe{
#' \item{patid}{Subject ID.}
#' \item{reptind}{Indicator for a subsequent event: 0 for the first observed
#' event and 1 for subsequent events for the subject.}
#' \item{date}{Julian date of the event.}
#' \item{age}{Integer age at the date of the event.}
#' }
"my.visit"


#' Sample Population Aggregate Dataset
#'
#' A sample data frame containing aggregated population census information.
#'
#' @format A data frame with $n$ rows and the following variables:
#' \describe{
#' \item{AGE}{Integer age.}
#' \item{MaleInd}{Indicator for male sex.}
#' \item{CityAInd}{Indicator for residence in region A.}
#' \item{CityBInd}{Indicator for residence in region B.}
#' \item{pop}{Number of subjects in the corresponding demographic group.}
#' }
"my.agge.info"

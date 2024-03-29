% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_functions.R
\name{incidence_year2day}
\alias{incidence_year2day}
\title{Converts annual incidence to daily}
\usage{
incidence_year2day(incidence)
}
\arguments{
\item{incidence}{annual incidence (cases per 1000 person year)}
}
\value{
Scalar representing daily incidence (cases per 1 person day).
h=incidence/1000/365
}
\description{
Converts annual incidence to daily
}

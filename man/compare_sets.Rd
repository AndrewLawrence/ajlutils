% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{compare_sets}
\alias{compare_sets}
\title{Compare two sets}
\usage{
compare_sets(x, y, list = FALSE)
}
\arguments{
\item{x}{a vector of elements to compare with y}

\item{y}{a vector of elements to compare with x}

\item{list}{boolean should matching and non-matching elements be returned?
Or just a count?}
}
\value{
a list of set sizes,
    or a list containing both the set sizes and lists of elements
}
\description{
Will convert input to sets of unique elements.
}

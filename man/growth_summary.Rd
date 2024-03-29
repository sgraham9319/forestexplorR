% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/growth_summary.R
\name{growth_summary}
\alias{growth_summary}
\title{Calculate annual growth}
\usage{
growth_summary(data)
}
\arguments{
\item{data}{A dataframe containing repeated measurements of growth. Must 
have a column named "tree_id" that contains unique tree ID values, a column
named "year" containing the year of the measurement, and a column named
"dbh" containing the dbh measurements.See built-in dataset \code{tree} for
an example.}
}
\value{
A dataframe containing growth rate measurements for each tree:
\describe{
  \item{tree_id}{unique identification code of the tree}
  \item{stand_id}{name of stand in which tree is located - only appears if 
  input dataframe contains a \code{stand_id} column}
  \item{species}{species identity of tree as a four letter code}
  \item{first_record}{year of first size measurement of the tree}
  \item{last_record}{year of last size measurement of the tree}
  \item{begin_size}{first size measurement as diameter at breast height,
  in cm}
  \item{final_size}{last size measurement as diameter at breast height,
  in cm}
  \item{mean_size}{average size across measurements as diameter at breast
  height, in cm}
  \item{midpoint_size}{average of maximum and minimum size across
  measurements as diameter at breast height, in cm}
  \item{annual_growth}{annual growth rate as average yearly increase in 
  diameter at breast height, in cm/year}
  \item{annual_bai}{annual basal area increment as average yearly increase
  in area at breast height, in cm^2/year}
  \item{size_corr_growth}{square root of \code{annual_growth} divided by
  \code{begin_size} to give somewhat normally distributed growth rates}
  \item{size_corr_growth_basal}{square root of \code{annual_bai} divided by
  \code{begin_size} to give somewhat normally distributed growth rates}
}
}
\description{
Requires a data frame where separate measurements of the same tree appear
in different rows and returns a data frame with a single row for each tree
containing annual growth measurements.
}
\examples{
growth_summary(tree)
}

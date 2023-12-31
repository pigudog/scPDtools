% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCP-plot.R
\name{segementsDf}
\alias{segementsDf}
\title{Shorten and offset the segment}
\usage{
segementsDf(data, shorten_start, shorten_end, offset)
}
\arguments{
\item{offset}{}
}
\description{
Shorten and offset the segment
}
\examples{
library(ggplot2)
tempNodes <- data.frame("x" = c(10, 40), "y" = c(10, 30))
data <- data.frame("x" = c(10, 40), "y" = c(10, 30), "xend" = c(40, 10), "yend" = c(30, 10))
ggplot(tempNodes, aes(x = x, y = y)) +
  geom_point(size = 12) +
  xlim(0, 50) +
  ylim(0, 50) +
  geom_segment(data = data, aes(x = x, xend = xend, y = y, yend = yend))

ggplot(tempNodes, aes(x = x, y = y)) +
  geom_point(size = 12) +
  xlim(0, 50) +
  ylim(0, 50) +
  geom_segment(
    data = segementsDf(data, shorten_start = 2, shorten_end = 3, offset = 1),
    aes(x = x, xend = xend, y = y, yend = yend)
  )
}

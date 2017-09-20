\name{rain}
\docType{data}
\alias{rain}

\title{Precipitation Observations and Forecasts for Innsbruck}

\description{
Accumulated 18-30 hour precipitation ensemble forecasts and corresponding observations at Innsbruck. The dataset includes GEFS
reforecasts (Hamill et al. 2013) and observations from SYNOP station Innsbruck
Airport (11120) from 2000-01-02 to 2016-01-01.
}

\usage{data("temp")}

\format{
  A data frame with 2749 rows. The first column (\code{rain}) are 12-hour
  accumulated precipitation observations. Columns 2-12 (\code{rainfc})
  are 18-30 hour accumulated precipitation forecasts from the individual
  ensemble members.
}

\source{
Observations: http://www.ogimet.com/synops.phtml.en

Reforecasts: http://www.esrl.noaa.gov/psd/forecasts/reforecast2/
}

\references{
  Hamill TM, Bates GT, Whitaker JS, Murray DR, Fiorino M, Galarneau Jr TJ,
    Zhu Y, Lapenta W (2013). NOAA's Second-Generation Global Medium-Range
    Ensemble Reforecast Data Set. \emph{Bulletin of the American Meteorological
    Society}, 94(10), 1553-1565.

  Vannitsem S, Wilks DS, Messner JW (2017). Statistical Postprocessing of
    Ensemble Forecasts, \emph{Elsevier}, to appear.
}

\examples{
## Diagnostic plots similar to Figure 8 in Vannitsem et al. ##

## load and prepare data
data("rain")

rain <- sqrt(rain)

rain$ensmean <- apply(rain[,2:12], 1, mean)
rain$enssd <- apply(rain[,2:12], 1, sd)

## Scatterplot of precipitation by ensemble mean
plot(rain~ensmean, rain, col = gray(0.2, alpha = 0.4), 
  main = "Scatterplot")
abline(0, 1, lty = 2)

## Verification rank histogram
rank <- apply(rain[,1:12], 1, rank)[1,]
hist(rank, breaks = 0:12 + 0.5, main = "Verification Rank Histogram")

## Spread skill relationship
sdcat <- cut(rain$enssd, quantile(rain$enssd, seq(0, 1, 0.2)))
boxplot(abs(rain-ensmean)~sdcat, rain, ylab = "absolute error",
xlab = "ensemble standard deviation", main = "Spread-Skill")

## Histogram
hist(rain$rain, xlab = "square root of precipitation", main = "Histogram")
}

\keyword{datasets}

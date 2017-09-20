\name{temp}
\docType{data}
\alias{temp}

\title{Minimum Temperature Observations and Forecasts for Innsbruck}

\description{
18-30 hour minimum temperature ensemble forecasts and corresponding observations at Innsbruck. The dataset includes GEFS reforecasts (Hamill et al. 2013) and observations from the SYNOP station Innsbruck Airport (11120) from 2000-01-02 
to 2016-01-01.
}

\usage{data("temp")}

\format{
  A data frame with 2749 rows. The first column (\code{temp}) are 12-hour
  minimum temperature observations. Columns 2-12 (\code{tempfc})
  are 18-30 hour minimum temperature forecasts from the individual
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
## Diagnostic plots similar to Figure 1 and 3 in Vannitsem et al. ##

## load and prepare data
data("temp")

temp$ensmean <- apply(temp[,2:12], 1, mean)
temp$enssd <- apply(temp[,2:12], 1, sd)

## Scatterplot of minimum temperature observation by ensemble mean
plot(temp~ensmean, temp, main = "Scatterplot")
abline(0, 1, lty = 2)

## Verification rank histogram
rank <- apply(temp[,1:12], 1, rank)[1,]
hist(rank, breaks = 0:12 + 0.5, main = "Verification Rank Histogram")

## Spread skill relationship
sdcat <- cut(temp$enssd, breaks = quantile(temp$enssd, seq(0, 1, 0.2)))
boxplot(abs(temp-ensmean)~sdcat, temp, ylab = "absolute error",
xlab = "ensemble standard deviation", main = "Spread-Skill")

## Histogram
hist(temp$temp, xlab = "minimum temperature", main = "Histogram")
}

\keyword{datasets}

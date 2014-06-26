\name{aemet}
\Rdversion{1.1}
\alias{aemet}
\docType{data}
\title{aemet data}
\description{
Series of daily summaries of 73 spanish weather stations selected for the period 1980-2009. 
The dataset contains geographic information of each station and the average for the period 1980-2009 of daily temperature, daily precipitation and daily wind speed.\cr
Meteorological State Agency of Spain (AEMET), \url{http://www.aemet.es/}.\cr
Government of Spain.\cr
}
\usage{data(aemet)}
\format{
Elements of aemet:\cr
 \code{..$df:} Dataframe with information of each wheather station:\cr
 \code{ind:}{ Indicated weather station.\cr}
 \code{name:}{ Station Name.  36 marked UTF-8 strings \cr}
 \code{province:}{  Province (region) of Spain. 36 marked UTF-8 strings \cr}
 \code{altitude:}{  Altitude of the station (in meters).\cr}
 \code{year.ini:}{ Start year.\cr}
\code{year.end:}{ End year.\cr} 
\code{longitude:}{ x geographic coordinate of the station  (in decimal degrees) .\cr}
\code{latitude:}{ y geographic coordinate of the station (in decimal degrees) .\cr}
\cr       
 The functional variables:  \cr
  \code{...$temp}: mean curve of the average daily temperature for the period 1980-2009 (in degrees Celsius, marked with UTF-8 string).\cr
  \code{...$wind.speed}: mean curve of the average daily wind speed  for the period 1980-2009 (in m/s).\cr
 \code{...$logprec}: mean curve of the log precipitation for the period 1980-2009 (in log mm).\cr
}
\details{
It marks 36 UTF-8 string of names of stations and 3 UTF-8 string names of provinces through the function \code{\link{iconv}}.\cr

In leap years temperatures for February 28 and 29 were averaged.\cr
Negligible precipitation (less than 1 tenth of mm) is replaced  by \code{0.05} and 
no precipitation (0.0 mm) is replaced  by \code{0.01}.  Then the logarithm is applied.
}
\source{ The data were obtained from the FTP of AEMET in 2009.  
}

\author{Manuel Febrero Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}}
\examples{
data(aemet)
names(aemet)
names(aemet$df)
par(mfrow=c(3,1))
plot(aemet$temp)
plot(aemet$wind.speed)
plot(aemet$logprec)}
\keyword{datasets}
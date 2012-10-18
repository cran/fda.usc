\name{phoneme}
\Rdversion{1.1}
\alias{phoneme}
\docType{data}
\title{
phoneme data
}
\description{
Phoneme curves.
}
\usage{data(phoneme)}
\format{
Elements of phoneme:\cr
  \code{..$learn}: learning sample of curves. \code{fdata} class object with: i.- \code{"data"}:  Matrix of class \code{fdata} with 250 curves (rows) discretized in 150 points or argvals (columns).\cr, ii.- \code{"argvals"}, iii.- \code{"rangeval"}: range(\code{"argvals"}), iv.- \code{"names"} list with:  \code{main} an overall title "Phoneme learn",	\code{xlab} title for \code{x} axis "frequencies" and \code{ylab} title for \code{y} axis  "log-periodograms".\cr
  \cr
  \code{..$test}: testing sample of curves. \code{fdata} class object with: i.- \code{"data"}: Matrix of class \code{fdata} with 250 curves (rows) discretized in 150 points or argvals (columns).\cr, ii.- \code{"argvals"}, iii.- \code{"rangeval"}: range(\code{"argvals"}), iv.- \code{"names"} list with:  \code{main} an overall title "Phoneme learn",	\code{xlab} title for \code{x} axis "frequencies" and \code{ylab} title for \code{y} axis  "log-periodograms".\cr
  \cr
  \code{..$classlearn}:learning class numbers (as factor). Factor levels: "sh" 1, "iy" 2, "dcl" 3, "aa" 4 and "ao" 5.\cr
  \cr
  \code{..$classtest}: testing class numbers (as factor). Factor levels: "sh" 1, "iy" 2, "dcl" 3, "aa" 4 and "ao" 5.\cr
}
\details{
The following instructions have been used file: \cr
\url{http://www.math.univ-toulouse.fr/staph/npfda/npfda-phondiscRS.txt}\cr
of   \code{Phoneme dataset} file.
}
\source{
\url{http://www.math.univ-toulouse.fr/staph/npfda/npfda-datasets.html}
}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{NPFDA in practice}. Free access on line at
\url{http://www.lsp.ups-tlse.fr/staph/npfda/}
}
\author{
Manuel Febrero-Bande and Manuel Oviedo de la Fuente <manuel.oviedo@usc.es>
}
\examples{
data(phoneme)
names(phoneme)
names(phoneme$learn)
class(phoneme$learn)
dim(phoneme$learn)
table(phoneme$classlearn)
}
\keyword{datasets}
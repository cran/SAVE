\name{Spotweld}
\alias{spotweldfield}
\alias{spotweldmodel}
\docType{data}
\title{
Spotweld data
}
\description{
These data relate to an engineering application where two metal sheets are compressed by water-cooled copper electrodes under an applied load L and the resistance spot welding is measure. Two data sets are provided, one containing observations from the real experiment (\code{spotweldfield}) and another one containing observations from the simulator (i.e. the computer model, \code{spotweldmodel}).

Details about this data set can be found in Bayarri et al (2007). }
\usage{
	data(spotweldfield)
	
	data(spotweldmodel)
}
\format{
  Two data frame with observations on the following variables.
  \describe{
    \item{\code{L}}{Controllable input (Load)}
    \item{\code{C}}{Controllable input (The direct current of magnitude supplied to the sheets)}
    \item{\code{G}}{Controllable input}
    \item{\code{N}}{Response variable (Resistance)}
    \item{\code{t}}{Calibration input (hence only in \code{spotweldmodel})}
  }
}
\references{
Bayarri MJ, Berger JO, Paulo R, Sacks J, Cafeo JA, Cavendish J, Lin CH, Tu J (2007). A Framework for Validation of Computer Models. Technometrics, 49, 138-154.
}

\keyword{datasets}

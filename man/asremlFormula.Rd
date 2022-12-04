\name{asremlFormula}
\alias{asremlFormula}

\title{
Obtain the best formula for asreml designed for the Fully Automated MET Pipeline
}
\description{
Based on data, an asreml formula and constrainst we find the best asreml formula to fit and avoid singularities.
}
\usage{
asremlFormula(fixed, random, rcov, dat, minRandomLevels=NULL, minResidualLevels=NULL, exchangeRandomEffects=NULL, 
              exchangeResidualEffects=NULL,customRandomLevels=NULL, customResidualLevels=NULL,
              xCoordinate= "ROWf",yCoordinate ="RANGEf", doubleConstraintRandom=c("ROWf","RANGEf"))
}

\arguments{
  \item{fixed}{
the fixed effects formula of the form y ~ x
}
\item{random}{
 the random effects formula of the form ~ diag(f):x or ~ x
}
\item{rcov}{
the residual effects formula of the form ~ diag(f):units or ~ units
}
\item{dat}{
data frame containing the response and effects specified in the formulas
}
\item{minResidualLevels}{
a list specifying the minimum number of levels that need to be present for each residual effect.
                         By default at least 2 levels need to be present to fit such residuals effects. I.e. minResidualLevels=list(Rowf=2) 
}
\item{minRandomLevels}{
a list specifying the minimum number of levels that need to be present for each random effect.
                         By default at least 2 levels need to be present to fit such residuals effects. I.e. minRandomLevels=list(Rowf=2) 
}
\item{exchangeRandomEffects}{
a list specifying if the constraints to be applied to two random effects should be applied in an inverse way.  
}
\item{exchangeResidualEffects}{
a list specifying if the constraints to be applied to two random effects should be applied in an inverse way.
}
\item{customRandomLevels}{
a list coming in the format of the "used" slot from the same function to indicated the levels of each random effect of the formula and force the results. 
}
\item{customResidualLevels}{
a list coming in the format of the "used" slot from the same function to indicated the levels of each residual effect of the formula and force the results.  
}
\item{xCoordinate}{
character string specifying the column name in the dataset corresponding to the x coordinate to know the spatial form.  
}
\item{yCoordinate}{
character string specifying the column name in the dataset corresponding to the y coordinate to know the spatial form.  
}
\item{doubleConstraintRandom}{
a list specifying the pair of random effects that should be constrained simultaneously. 
}

}
\details{

At high marker density, the relationship matrix is estimated as \eqn{A=W W'/c}, where \eqn{W_{ik} = X_{ik} + 1 - 2 p_k} and \eqn{p_k} is the frequency of the 1 allele at marker k.  By using a normalization constant of \eqn{c = 2 \sum_k {p_k (1-p_k)}}, the mean of the diagonal elements is \eqn{1 + f} (Endelman and Jannink 2012).

}
\value{
If return.imputed = FALSE, the \eqn{n \times n} additive relationship matrix is returned.

If return.imputed = TRUE, the function returns a list containing
\describe{
\item{$out}{A character vector of same length as \code{Source} in which the \code{Search}-values have been replaced by the corresponding \code{Replace}-values.}
}
}

\references{
None
}

\examples{

# library(plyr)
# head(soyMetData)
# subdata <- soyMetData
# 
# iTrait <- "YIELD"
# mde <- asremlFormula(fixed=as.formula(paste(iTrait,"~ 1")),
#                      random=~ at(FIELDINST):TRIAL + at(FIELDINST):ROWf + at(FIELDINST):RANGEf,
#                      rcov=~at(FIELDINST):id(ROWf):id(RANGEf),
#                      dat=droplevels(subdata[!is.na(subdata[,iTrait]),]),
#                      
#                      minRandomLevels=list(TRIAL=2,ROWf= 6, RANGEf=6),
#                      minResidualLevels=list(ROWf=6, RANGEf=6),
#                      
#                      exchangeRandomEffects=list(ROWf="RANGEf", RANGEf="ROWf"),
#                      
#                      exchangeResidualEffects=list(ROWf="RANGEf", RANGEf="ROWf"),
#                      
#                      customRandomLevels=NULL, customResidualLevels=NULL,
#                      
#                      xCoordinate= "ROWf",yCoordinate ="RANGEf",
#                      doubleConstraintRandom=c("ROWf","RANGEf"))
# as.formula(mde$fixed)
# as.formula(mde$random)
# as.formula(mde$rcov)

}
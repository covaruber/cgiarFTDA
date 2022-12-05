\name{geneticFormulation}
\alias{geneticFormulation}

\title{
Expands the genetic formula for an ASReml model
}
\description{
Given the parameters it expands the genetic formula for a factor analytic model (rr), diagonal model (dg) or a compound simmetry model (cs) in a additive + non.additive, a pure additive model, or single genetic model.
}
\usage{
geneticFormulation(mymodel,noFa=c(1,1),overlayColumnName="half", materialNameColumnName="material_id",
                               femaleMaterialNameColumnName="female_material_id",
                               maleMaterialNameColumnName="male_material_id",
                               pedargName="pedarg")
}

\arguments{
  \item{mymodel}{
a character vector with 4 terms:
               the model to fit: one of the 3 available "rr" (factor analytic), "dg" (diagonal), or "cs" (compound simmetry)
               level for fitting: one of the 2 available "sib" (per se level), "prt" (parental)
               additive term: one of the 2 available "vm" (include additive) or "na" (single genetic)
               non-additive term: one of the 2 available "ide" (include non-additive) or "na" (don't include non-additive)
}
\item{noFa}{
a numeric vector with 2 values indicating the number of factors to add to the additive and non-aaditive terms in the factor analytic models.
}
\item{overlayColumnName}{
character term indicating the name of the column that is used for the overlayed model.
}
\item{materialNameColumnName}{
character term indicating the name of the column that contains the ID for the material evaluated per se.
}
\item{femaleMaterialNameColumnName}{
character term indicating the name of the column that contains the female material ID of the material evaluated per se. 
}
\item{maleMaterialNameColumnName}{
character term indicating the name of the column that contains the male material ID of the material evaluated per se. 
}
\item{pedargName}{
name of pedigree argument used in your modeling to include everything at once. If the value is NA it assumes no pedarg to be included in the formula (ASReml3) 
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

# fa model sample
as.formula(geneticFormulation(mymodel=c("rr","sib","vm","ide"),noFa=c(2,1)))
as.formula(geneticFormulation(mymodel=c("rr","prt","vm","ide"),noFa=c(2,1)))
as.formula(geneticFormulation(mymodel=c("rr","sib","na","na"),noFa=c(2,1)))
as.formula(geneticFormulation(mymodel=c("rr","prt","na","na"),noFa=c(2,1)))
# dg model sample
as.formula(geneticFormulation(mymodel=c("dg","sib","vm","ide"),noFa=c(2,1)))
as.formula(geneticFormulation(mymodel=c("dg","prt","vm","ide"),noFa=c(2,1)))
as.formula(geneticFormulation(mymodel=c("dg","sib","na","na"),noFa=c(2,1)))
as.formula(geneticFormulation(mymodel=c("dg","prt","na","na"),noFa=c(2,1)))
# cs model sample
as.formula(geneticFormulation(mymodel=c("cs","sib","vm","ide"),noFa=c(2,1)))
as.formula(geneticFormulation(mymodel=c("cs","prt","vm","ide"),noFa=c(2,1)))
as.formula(geneticFormulation(mymodel=c("cs","sib","na","na"),noFa=c(2,1)))
as.formula(geneticFormulation(mymodel=c("cs","prt","na","na"),noFa=c(2,1)))
# ASReml4 vs ASReml3 formula
as.formula(geneticFormulation(mymodel=c("rr","sib","vm","ide"),noFa=c(2,1), pedargName="pedarg"))
as.formula(geneticFormulation(mymodel=c("rr","sib","ped","ide"),noFa=c(2,1), pedargName=NA))


}
\name{pheno}
\alias{pheno}
\docType{data}
\title{
A dummy phenotype dataset that provides an example of the input phenotype data used by the package
}
\description{
A dummy dataset of 5 phenotypes measured in 150 individuals. The data has been generated to yield significant results for SNP1 and SNP2 of the snps dataset.  The first two columns have been generated as alpha + beta1*snp + beta2*snp2 + error (with different alphas, betas and errors for each phenotype), the third has been generated as alpha + beta1*testPheno2 + beta2*snp3 + error, the fourth column is the results of sample of a binomial distributioni correlated with testPheno3, and the final column is the 1st PC of the principal component analysis of the snps matrix.
}
\format{
  A matrix with 150 phenotype observations.
 \describe{
    \item{\code{testPheno1}}{a numeric vector}
    \item{\code{testPheno2}}{a numeric vector}
    \item{\code{testPheno3}}{a numeric vector}
    \item{\code{testPheno4}}{a numeric vector}
    \item{\code{testPheno5}}{a numeric vector}
  }
}
\details{
Please note the following IMPORTANT issue: the 'pheno' matrix has both column names and row names! the column names MUST be the names of the phenotypes and the row names MUST be the codes representing each individual in the pheno matrix, one individual for each row. Both row names and column names are extracted by the main function and are therefore mandatory}
\examples{
data(pheno)
head(pheno)
dimnames(pheno)[[1]] # the row names
dimnames(pheno)[[2]] # the column names
}
\keyword{datasets}

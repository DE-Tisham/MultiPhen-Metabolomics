\name{snps}
\alias{snps}
\docType{data}
\title{
A dummy snp dataset that provides an example of the input snp data used by the package
}
\description{
A dummy dataset of three SNPs, as a matrix of 3 column and 150 rows. The genotypes are in 0/1/2 format (0 for ``AA'', 1 for ``Aa'' and 2 for ``aa'', where A and a correspond, arbitrarily, to the two alleles). The data has been randomly generated, for instructional purposes only, and do not yield a significant association with any of the example phenotypes.
}
\usage{data(snps)}
\format{
  A data frame with 150 genotype observations.
  \describe{
    \item{\code{rsID1}}{a numeric vector}
    \item{\code{rsID2}}{a numeric vector}
    \item{\code{rsID3}}{a numeric vector}
  }
}
\details{
The 150 genotypes in 'snp' correspond to the phenotype data on 150 individuals in `pheno', i.e. one individual for each line. Please note the following important points: genotype data must be in matrix format, with one row for each individual and as many columns for each SNP. In the case of one single genotype the data must still conform to this format, as a matrix of as many rows as individuals and one single column for the one genotype present. A second important point is that the column names must be the rsID of the SNP for genotypes in the 0/1/2 format. Further options of genotype format (incorporating raw genotype data, and CNV genotypes) will be available and documented in future releases.
}
\examples{
data(snps)
dim(snps)
colnames(snps)
}
\keyword{datasets}

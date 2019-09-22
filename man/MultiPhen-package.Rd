\name{MultiPhen-package}
\alias{MultiPhen-package}
\alias{MultiPhen}
\docType{package}
\title{
MultiPhen, a package for the genetic association testing of multiple phenotypes
}
\description{
MultiPhen is the package containing the function mPhen, which performs association testing between genetic variants (SNPs; CNVs to be added soon) and multiple phenotypes. The primary purpose is for modelling and testing multiple phenotypes jointly by performing an ordinal regression where SNPs are treated as the outcome and multiple phenotypes are predictors; this can have large increases in statistical power to detect genotype-phenotype associations over the univariate approach. However, mPhen can also be used to perform standard univariate linear regression (SNP as predictor) and univariate ordinal regression (SNP as outcome) on the phenotypes under study. mPhen can be applied to genotyped or imputed data. From version 0.4 the option ``multiPhenTest'' is now called ``JointModel'', and its default is now ``TRUE''
}
\details{
\tabular{ll}{
Package: \tab MultiPhen\cr
Type: \tab Package\cr
Version: \tab 2.0.0\cr
Date: \tab 2013-12-09\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Lachlan Coin, Federico Calboli, Clive Hoggart, Paul O'Reilly, Yotsawat Pomyen.

Maintainer: Federico Calboli <f.calboli@imperial.ac.uk>}
\references{
O'Reilly et al. 2012. MultiPhen: Joint model of multiple phenotypes can increase discovery in GWAS. http://dx.plos.org/10.1371/journal.pone.0034861
}
\keyword{ package }
\keyword{ genetics }
\keyword{ GWAS }
\keyword{ regression }
\keyword{ multiple phenotypes }
\keyword{ association }


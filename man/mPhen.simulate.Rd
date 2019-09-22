\name{mPhen.simulate}
\alias{mPhen.simulate}
\title{Simulates phenotypes according to a correlation structure.}
\description{This function simulates phenotypes based on a pre-defined correlation structure (which can also be obtained from 
mPhen.sampleCovar), and a genetic effect x.  The function works by sampling a phenotype from a correlation matrix in a linearly 
transformed space such that the genetic effect direction is only in the direction of the x-axis, then transforming back into the
 original space.  If inverse is TRUE,  then the phenotypes are sampled first with no genetic effect, then the genotype is sampled 
according to the effect direction.
}
\usage{mPhen.simulate(x,sample_names, covar,effDir,
varexp,inverse=FALSE, geno.link="gaussian",
effDirInReverseEigenspace=FALSE, freq = 0.1)}
\arguments{
\item{x}{A vector of genotype effect. If a single SNP has an effect, this will just be genotypes at this SNP}
\item{sample_names}{Vector of sample names, should have same length as x}
\item{covar}{Covariance of phenotypes, should be an n by n matrix, where n is the number of phenotypes to simulate.}
\item{effDir}{Direction in phenotype space in which to simulate the effect}
\item{varexp}{The proportion of variance of phenotype variation in the target direction explained by the genotypic effect overall.}
\item{inverse}{If TRUE, then simulates correlated phenotypes, and then simulates genotypes from phenotypes in specified direction.  Otherwise, simulates correlated phenotypes with direction of effect based on input genotypes.}
\item{geno.link}{Only applicable if inverse = TRUE, in which case it specifies a link function for genotypes.  Can be binomial, gaussian or ordinal.}
\item{effDirInReverseEigenspace}{If TRUE, then effDir is interpreted as eigenvector weights, ordered from eigenvector with smallest eigenvalue to eigenvector with biggest eigenvector (i.e. in reverse direction).  This is useful if you want to simulate directs which are in the least variable axis of variation.}
\item{freq}{If inverse = TRUE, then this is target allele frequency.}
}
\value{Returns a list of three elements: pheno - a matrix of phenotype values, with phenotypes in columns and samples in rows; and geno- a vector genotypes, either sampled if inverse=TRUE, or the original genetic effect x; and limit - this is a default list of phenotypes to include in subsequent association analysis, as well as covariates (empty) and residuals (empty). This object can then be used by mPhen.preparePheno}

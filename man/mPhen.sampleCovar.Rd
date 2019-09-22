
\name{mPhen.sampleCovar}
\alias{mPhen.sampleCovar}
\title{Generates a covariance matrix.}
\description{This function can be use to sample covariance matrices.  This is useful when simulating data to test Multiphenotype based association strategies. This function lets the user decide on the orthoganality within 'blocks' and between 'blocks' of correlated variables/}
\usage{mPhen.sampleCovar(noPhenos,blockSize, orthogAll = c(0.9,0.5),
dirichletScale = 50,  resample = FALSE, 
sd = rgamma(noPhenos,shape=10,rate = 10))}
\arguments{
\item{noPhenos}{The number of phenotypes to simulate}
\item{blockSize}{The number of phenotypes per covariance block}
\item{orthogAll}{The orthogonality relationships between and within blocks expressed as a number on the interval (0,1).  A number closer to one indicates closer to orthogonality, whereas 0 indicates non-orthogonality.  First number is orthogonality between blocks, second is orthogonality within blocks.}
\item{dirichletScale}{When sampling off diagonal elements of the cholesky decomposition, how much deviation from uniform to allow. Should be a number in interval (0,+Inf).  Smaller value leads to greater variation}
\item{resample}{Whether to randomly shuffle phenotype columns after sampling. }
\item{sd}{Standard deviation for each phenotype}
}
\value{Simulated covariance matrix}


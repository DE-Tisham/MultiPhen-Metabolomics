\name{mPhen.sampleGeno}
\alias{mPhen.sampleGeno}
\title{Sample genotypes}
\description{...}
\usage{mPhen.sampleGeno(n = 100, sampSize = 100, chr="0",pos = 1:n, 
snpids = paste(chr,pos,sep="_"),meanAlleleFreq=0.2, mu = 10,
samples =paste("id",1:sampSize,sep="_"),imputed = FALSE, dirichlet = 1)}
\arguments{
\item{n}{Number of genotypes to sample}
\item{sampSize}{Number of individuals to sample}
\item{chr}{Name of chromosome}
\item{pos}{Positions of genotypes on chromosome}
\item{snpids}{Ids of genotypes}
\item{meanAlleleFreq}{The mean allele frequency to simulate}
\item{mu}{A weight parameter which controls how close to the meanAlleleFreq the allele frequencies are sampled, via a beta distribution.  A higher number implies allele frequencies stay closer to mean}
\item{samples}{The sample ids}
\item{imputed}{Whether to simulate imputed data}
\item{dirichlet}{The weight of a dirichlet distribution used to simulated imputed data}
}
\value{Returns matrix of genotypes, with individuals by rows, and snps by column, or a 3 dimensional array if imputed is TRUE}


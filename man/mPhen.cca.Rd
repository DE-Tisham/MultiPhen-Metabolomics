\name{mPhen.cca}
\alias{mPhen.cca}
\title{Sparse canonical correlation analysis}
\description{Carries out a sparse canonical correlation analysis}
\usage{mPhen.cca(genoData, phenoObject, opts =mPhen.options("regression"), 
		         subinds = 1:(dim(genoData)[1]),
		vs.G = opts$mPhen.variable.selection,
                 vs.P = opts$mPhen.variable.selection
)}
\arguments{
\item{genoData}{A 2 dimensional array. The first dimension (rows) corresponds to individuals, and row.names are inidividual IDs. The second dimension corresponds to SNPS, with col.names equal to the snp identifiers. The entries are either genotypes, or expected genotypes.}
\item{phenoObject}{A phenotype object prepared by mPhen.preparePheno.}
\item{opts}{A list of options, which is obtained from mPhen.options("regression").  To get more information about these options, type 
mPhen.options("regression",descr=TRUE).  }
\item{vs.G}{If true performs variable selection on genotypes.  Is equal to  opts$mPhen.variable_selection by default}
\item{vs.P}{If true performs variable selection on phenotypes.  Is equal to  opts$mPhen.variable_selection by default}
\item{subinds}{This indicates the indices of individuals to include in the analysis. It is possible to have repeat indices (i.e. for bootstrap)}
}
\value{A list with following elements
\item{betasp}{Phenotype weights}
\item{betasg}{Genotype weights}
\item{resultsGeno}{Single (combined) phenotype analysis results against all genotypes}
\item{resultsPheno}{Multiple phenotype analysis results against single (combined) genotype}
 }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetics }
\keyword{ GWAS }
\keyword{ regression }
\keyword{ multiple phenotypes } 
\keyword{ association }
\keyword{ file }
\keyword{ connection }

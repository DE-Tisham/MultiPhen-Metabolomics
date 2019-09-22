\name{mPhen.assoc}
\alias{mPhen.assoc}
\title{
A function for the genetic association testing of multiple phenotypes
}
\alias{mPhen.assoc}
\description{This function is called by mPhen.  If you are doing association on multiple batches 
of genotype data, it is more efficient to use this function, and to pre-prepare a 'phenoObject' object once and then use this function
}
\usage{mPhen.assoc(genoData, phenoObject,opts = mPhen.options("regression"),
		 subinds =1:(dim(genoData)[1]))}
\arguments{
\item{phenoObject}{
A phenotype object prepared by mPhen.preparePheno
}
\item{genoData}{
This can be obtained from mPhen.readGenoConnection(...)$genoData.
It is either a matrix (for directly measured genotypes) or 3 dimensional array (for imputed genotypes). The first dimension (rows) corresponds to individuals, and row.names are inidividual IDs. The second dimension corresponds to SNPS, with col.names equal to the snp identifiers. For directly measured genotypes, the value in each cell is a numeric genotype( i.e. AA = 0, AB=1, BB = 2).  For imputed genotypes, the 3rd dimension corresponds to genotypes (with dimnames(genoData)[[3]]) equal to a numeric vector corresponding to genotype values. The values in these cells are the probability of each genotype multiplied by 1000.  For copy number genotypes the numeric values correspond to numbers of copies.  An example provided by 'snps'.}
\item{opts}{A list of options, which is obtained from mPhen.options("regression").  To get more information about these options, type 
mPhen.options("regression",descr=TRUE)}
\item{subinds}{This indicates the indices of individuals to include in the analysis. It is possible to have repeat indices (i.e. for bootstrap)}
}
\value{
Returns a list, with two items.  The first item (Results) is a  Results is a 4 dimensional matrix, with
dimensions [strata, snps, phenotypes, result_type], where result_type includes beta, pvalue and Nobs.  The second item is a vector of minor allele frequencies.  
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

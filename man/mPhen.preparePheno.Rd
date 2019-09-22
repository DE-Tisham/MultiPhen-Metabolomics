\name{mPhen.preparePheno}
\alias{mPhen.preparePheno}
\title{Prepare phenotype data for analysis}
\description{This harmonises the phenotype data with genotype data, and
also extract the relevant columns from a larger phenotype matrix, and also pre-calculates stratification indices and residuals.  This is called by mPhen function, but quicker to do this just once for batched genotype data.
}
\usage{mPhen.preparePheno(phenoData,
pcs = NULL, indiv = if (is.null(pcs)) rownames(phenoData) else rownames(pcs), 
    opts = mPhen.options("regression"))}
\arguments{
\item{phenoData}{
This is typically the output of mPhen.readPhenoFiles(...) or mPhen.simulate(...). It is a list containing the elements phenoData$pheno and phenoData$limit.  The element phenoData$pheno is a matrix containing phenotype data, where each row corresponds to an individual and row.names are individual IDs. Each column contains data on a certain phenotype across the sample of individuals (can be quantitative, case/control or ordinal. Must be numeric); the column header provides the phenotype name. An example is provided by 'pheno'.   The element phenoData$limit is a list which specifies which of the phenotypes to use as covariates, variables to associate, etc, in the following way:

 limit$phenotypes - vector of phenotypes to be tested. If set to 'all' then all phenotypes are used.
 limit$covariates - vector of phenotypes to be considered as covariates to be controlled for in the regression ,
 limit$resids - vector of phenotypes to be considered as residuals, which is an alternative way to adjust for covariates, which pre-calculates offset terms to use in the per SNP regression 
 limit$strats - statification vector (i.e. cases/controls, exposed/not exposed, male/female etc).
 limit$excls - Exclusion vector, ie.  names of phenotypes which should be used as  exclusion criteria respectively.  Rows will be excluded if the value in any of exclusion columns is NA or 1

Alternatively, phenoData can simply be a matrix containing phenotype data, in which case, the default value of limit used is limit = list(phenotypes="all")
}
\item{pcs}{This is the genotype pcs which should be used in the analysis.  If specified it should be in the same sample order as 'inidv'.  The user still needs to specify in the phenoData$limit$covariate or phenoData$limit$resids the names of the PCs they wish to include (i.e. covariate = c("PC1","PC2")).  These would typically be obtained from mPhen.readGenotypes.}
\item{indiv}{individuals to be used.  If unspecified the fucntion defaults to using all individuals}  
\item{opts}{A list of options, which is obtained from mPhen.options("regression").  To get more information about these options, type 
mPhen.options("regression",descr=TRUE)}
}
\value{A list object, which can then be used in mPhen.assoc.}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetics }
\keyword{ GWAS }
\keyword{ regression }
\keyword{ multiple phenotypes } 
\keyword{ association }
\keyword{ file }
\keyword{ connection }

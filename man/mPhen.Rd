\name{mPhen}
\alias{mPhen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for the genetic association testing of multiple phenotypes
}
\description{
mPhen performs association testing between genetic variants (SNPs; CNVs) and multiple phenotypes. The primary purpose is for modelling and testing multiple phenotypes jointly by performing an ordinal regression where SNPs are treated as the outcome and multiple phenotypes are predictors; this can have large increases in statistical power to detect genotype-phenotype associations over the univariate approach (method described in O'Reilly et al. 2012, see below). However, mPhen can also be used to perform standard univariate linear regression (SNP as predictor) and univariate ordinal regression (SNP as outcome) on the phenotypes under study. mPhen can be applied to directly genotyped or imputed data.
}
\usage{
mPhen(genoData, phenoData, phenotypes = "all",
      covariates = NULL, resids = NULL, strats = NULL, 
      opts = mPhen.options(c("regression","pheno.input")))
}
\arguments{
  \item{genoData}{
Either a matrix (for directly measured genotypes) or 3 dimensional array (for imputed genotypes). The first dimension (rows) corresponds to individuals, and row.names are inidividual IDs. The second dimension corresponds to SNPS, with col.names equal to the snp identifiers. For directly measured genotypes, the value in each cell is a numeric genotype( i.e. AA = 0, AB=1, BB = 2).  For imputed genotypes, the 3rd dimension corresponds to genotypes (with dimnames(genoData)[[3]]) equal to a numeric vector corresponding to genotype values. The values in these cells are the probability of each genotype multiplied by 1000.  For copy number genotypes the numeric values correspond to numbers of copies.  An example provided by 'snps' and 'snps.imputed'
}
\item{phenoData}{
Matrix containing phenotype data, where each row corresponds to an individual and row.names are individual IDs. Each column contains data on a certain phenotype across the sample of individuals (can be quantitative, case/control or ordinal. Must be numeric); the column header provides the phenotype name. An example is provided by 'pheno'. 
}
  \item{phenotypes}{
Vector of phenotype names, to be tested.  If value is 'all' then all phenotypes are included after removing covariates and residuals.
}
  \item{covariates}{
Vector of phenotypes, from phenoData, to be considered as covariates to be controlled for in the regression (Default is no covariates).
}
  \item{resids}{
Vector of residuals, from phenoData, alternative way to adjust for covariates, which pre-calculates offset terms to use in the per SNP regression (Default is no residuals).
}
  \item{strats}{
Statification vector (i.e. cases/controls, exposed/not exposed, male/female etc), from phenoData (Default is no stratification).
}
%  \item{excls}{
% Exclusion vector, ie.  names of columns of phenoData which should be used as  exclusion criteria respectively.  Rows will be excluded if the value in any of exclusion columns is NA or 1
%}  
\item{opts}{A list of options, which is obtained from mPhen.options(c("regression","pheno.input")).  To get more information about these options, type 
mPhen.options(c("regression","pheno.input"),descr=TRUE)}
}
\value{
Returns a list, with two items.  The first item (Results) is a  Results is a 4 dimensional matrix, with
dimensions [strata, snps, phenotypes, result_type], where result_type includes beta, pvalue and Nobs.  The second item is a vector of minor allele frequencies.
}
\references{
O'Reilly et al. 2012. MultiPhen: Joint model of multiple phenotypes can increase discovery in GWAS. http://dx.plos.org/10.1371/journal.pone.0034861}
\author{
Lachlan Coin, Federico Calboli, Clive Hoggart, Paul O'Reilly, Yotsawat Pomyen.

Maintainer, Federico Calboli \email{f.calboli@imperial.ac.uk} 
}
\note{
The user should remember that the genotype data file is always a matrix of at least a column, hence if taking a subset of 1 SNP in the non-imputed genotype data matrix, the option drop = FALSE should be used (see the example below)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(snps); data(snps.imputed); data(pheno)
opts = mPhen.options(c("regression","pheno.input"))
res = mPhen(snps, pheno, phenotypes = all,
      covariates = c('testPheno3', 'testPheno4'),opts = opts) 
# performs a MultiPhen analysis, with snp as outcome,
# and phenotypes testPheno1, testPheno2 as predictors, 
#with testPheno3 and testPheno4 as covariates using ordinal regression

res = mPhen(snps, pheno, phenotypes = c('testPheno1', 'testPheno2'), 
      covariates = c('testPheno3', 'testPheno4'), resids = 'testPheno5', opts = opts) 
# the same as above, with the fifth phenotype as residual 

res = mPhen(snps[,2, drop = FALSE], pheno, phenotypes = c('testPheno1', 'testPheno2'), 
      covariates = 'testPheno3',  opts = opts) 
# please note the use use of drop = FALSE if analysing only one SNP


res = mPhen(snps.imputed, pheno, phenotypes = c('testPheno1', 'testPheno2'), 
      covariates = 'testPheno3',  opts = opts) 
# for imputed data
}









% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetics }
\keyword{ GWAS }
\keyword{ regression }
\keyword{ multiple phenotypes } 
\keyword{ association }
% __ONLY ONE__ keyword per line




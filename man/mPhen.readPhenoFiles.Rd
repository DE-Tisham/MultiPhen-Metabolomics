\name{mPhen.readPhenoFiles}
\alias{mPhen.readPhenoFiles}
\title{Read and merge phenotype files}
\description{This helper function merges multiple phenotype files into a single phenotype matrix, and applies missing value and exclusion criteria}
\usage{mPhen.readPhenoFiles(phenoFiles,
	limitFile = getOption("mPhen.limitFile","./limit.txt"),
	 excludeFile =getOption("mPhen.excludeFile","./exclude.txt"),	
		opts = mPhen.options("pheno.input"))}
\arguments{
\item{phenoFiles}{A list of paths to phenotype files (can be more than 1)}
\item{excludeFile}{ A path to a file which lists ids to exclude from further analysis, or
alternatively is a two column file, with the first column of ids and a second column
of numberical values which are used in conjunction with opts$quantileThresh}
\item{limitFile}{
 As an alternative to specifying covariates, resids,strats and excl in mPhen.preparePheno(..), you can also specify this information via a limitfile, which is tab delimited file in which the first column specifies the type of variable to set (pheno,covar,resid,strat,excl), the second 
column specifies the phenotype name, and the third column optionally specifies a transformation
 Different lines can then be used for different values.  
The transformation syntax includes 'quantile' and 'factor', and also 'thresh_x_y' in which values less than x are coded 0 and greater than y are coded 1; and also 'toptail_x_y' where values less than x percentile are coded 0 and greater than y percentile are coded 1.
}
\item{opts}{A list of options, which is obtained from mPhen.options("pheno.input").  To get more information about these options, type 
mPhen.options("pheno.input",descr=TRUE)}
}
\value{An object consisting of a  single merged phenotype matrix, and also a 'limit' object which specifies phenotypes to include in analyses.  The limit is a list with the following entries.
	 phenotypes - vector of phenotypes to be tested. If set to 'all' then all phenotypes are used.
         covariates - vector of phenotypes to be considered as covariates to be controlled for in the regression ,
	 resids - vector of phenotypes to be considered as residuals, which is an alternative way to adjust for covariates, which pre-calculates offset terms to use in the per SNP regression 
       strats - statification vector (i.e. cases/controls, exposed/not exposed, male/female etc).

       excls - Exclusion vector, ie.  names of phenotypes which should be used as  exclusion criteria respectively.  Rows will be excluded if the value in any of exclusion columns is NA or 1

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

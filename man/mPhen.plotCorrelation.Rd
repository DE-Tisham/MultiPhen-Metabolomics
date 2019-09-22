\name{mPhen.plotCorrelation}
\alias{mPhen.plotCorrelation}
\title{Plots correlation between phenotype values conditional on genotype}
\description{Plots the corrleation between phenotype values, with different genotypes coloured differently.  Note that this will plot (dim(pheno_to_plot)[2] -1 )* dim(geno)[2] plots}
\usage{mPhen.plotCorrelation(pheno_to_plot,geno,title="",cex=0.25, cols =c(1, 2, 3))}
\arguments{
\item{pheno_to_plot}{the phenotype to plot}
\item{geno}{Genotypes to use for stratifying samples}
\item{title}{Title of plots}
\item{cex}{Scaling of points}
\item{cols}{thre colours to be used in the plot}
}
\value{None}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetics }
\keyword{ GWAS }
\keyword{ regression }
\keyword{ multiple phenotypes } 
\keyword{ association }
\keyword{ file }
\keyword{ connection }


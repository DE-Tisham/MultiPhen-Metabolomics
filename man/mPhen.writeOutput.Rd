\name{mPhen.writeOutput}
\alias{mPhen.writeOutput}
\title{Prepares output files and plots from MultiPhen results}
\description{Writes output to files defined in mPhen.openOutputConnection, and extracts pvalues and betas for
further plots.
}
\usage{mPhen.writeOutput(results,
output = getOption("mPhen.resultsName","resultsDir/"),geno = NULL,
towrite = list(long.txt = getOption("mPhen.writeLong",TRUE),
qc.txt =  getOption("mPhen.writeQC",FALSE), 
wide.txt = getOption("mPhen.writeWide",TRUE)),
toplot = list(.manh = TRUE, .qq = TRUE,.heatm = TRUE,
              %.heatmBeta = FALSE,
              .fprint = !is.null(geno)),
opts = mPhen.options("plot"))}

\arguments{
\item{results}{Output of mPhen.assoc}
\item{output}{Directory to write results, or object returned by mPhen.writeOutput(..)}
\item{towrite}{List specifying which formats to write output - long.txt and wide.txt for standard results; qc.txt for per-sample qc output.}
\item{toplot}{List specifying which formats to plot  output - .qq for qq plot, .manh for manhattan; .heatm for pvalue heatmap; .fprint for fingerprint plot}
\item{geno}{Genotype matrix.  Note that attr(geno,"closeConnection") controls whether plots are produced, as this indicates whether all batches of genotype data have been analysed}
\item{opts}{A list of options, which is obtained from mPhen.options("plot").  To get more information about these options, type
mPhen.options("plot",descr=TRUE)}
}
\value{Returns an outputConnection, which can be used to write further results.}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetics }
\keyword{ GWAS }
\keyword{ regression }
\keyword{ multiple phenotypes } 
\keyword{ association }
\keyword{ file }
\keyword{ connection }

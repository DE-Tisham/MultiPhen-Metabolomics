\name{mPhen.options}
\alias{mPhen.options}
\title{Retrieves default mPhen options, and descriptions.}
\description{ This command is used to get options which can be modified to control the behaviour of MultiPhen commands.
It provides a list of options which are relevant to a particular command.  For example, mPhen.assoc() has its behavioiur controlled by options
in mPhen.options("regression").  In order to get a list of all options, you can type mPhen.options(descr=TRUE).
}
\usage{mPhen.options(type=c("regression", "plot", "geno.input",
"pheno.input","meta.analysis","misc"), descr = FALSE)}
\arguments{
\item{type}{A value which can take any of the following values: "regression", "plot","geno.input","pheno.input","meta.analysis","misc".}
\item{descr}{If set to TRUE, then returns descriptions of all the options.  If FALSE, then returns the values of all the options.}
}
\value{A list of default option values.  Note, the default value for opts = mPhen.options("regression"), has opts$inverseRegress =TRUE, opts$JointModel = TRUE and opts$geno.link = "ordinal", which is the standard multiPhen model. }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetics }
\keyword{ GWAS }
\keyword{ regression }
\keyword{ multiple phenotypes } 
\keyword{ association }
\keyword{ file }
\keyword{ connection }

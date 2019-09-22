\name{mPhen.defineOptions}
\alias{mPhen.defineOptions}
\title{Defines and modifies  options. }
\description{ This reads options from the command-line, if provided.  Also it can replace references to system variables in option values, such as '*' or \$WORK, with fully qualified values.  Also translates coordinates which use Gb, Mb, Kb to integer values.  All defined options of the form mPhen.xxx are examined and modified.} 
\usage{mPhen.defineOptions(file = NULL, getOptionsFromCommandLine = TRUE)}
\arguments{
\item{file}{specified if the options are in a script file and not set manually}
\item{getOptionsFromCommandLine}{If running from a script using command Rscript, then this will read in command line options, such as \'--mPhen.logp=FALSE\'}
%\item{outDir}{If not equal to NULL, then all specified options are saved to outDir/options.txt.}
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

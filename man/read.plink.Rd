\name{read.plink}
\alias{read.plink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to read (small) binary PLINK binary files in a R session
}
\description{
\code{read.plink} is a convenience function designed to read PLINK binary files (i.e. files that end with the suffix ``.bed'') in a R session.  Please be aware that binary PLINK files are binary for a reasosn, i.e. to store genotype data in a compact way.  Once they are imported in R they exist in R in a un ``unpacked'' form, and can therefore be very big.  If the .bed file is big, or very big, the result will be that R will run out of memory and crash, or make the whole system slow or unresponsive.  It is MANDATORY that in the directory containing the binary file also reside two accesory files, with the same name as the binary file but with extensions .fam and .bim, both produced by PLINK.
}
\usage{
read.plink(root,indiv = NULL, opts = mPhen.options("geno.input"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{root}{
filename of the dataset in PLINK binary format, WITHOUT the .bed extension.  }
\item{indiv}{List of individuals, results will be in this order}
\item{opts}{List of options, use mPhen.options("geno.input",descr=TRUE) for more details of each option.}
}
\details{
Please note that, if the binary file is listed a ``mydata.bed'', the filename is ``mydata'', and the extension is ``.bed''.  In this case ``mydata'' would be used as \code{root} value.
}
\value{
A matrix of dimesions n by m, with n rows corresponding to the n individuals in the dataset, and m columns corresponding to the m markers.  The colnames are retrived from the .fam file, and (should) correspond to the markers' names.
}
\references{
The plink homepage is at: \url{http://pngu.mgh.harvard.edu/~purcell/plink/}
}
\author{
Federico Calboli \email{f.calboli@imperial.ac.uk}  
}
\note{
Please do note that the concept of a ``big'' binary file, or a binary file that is ``too big'' is purely dependent on the computer on which the code is running.  A computer with 512MB of RAM will stop being able to read in a whole binary file well before a 16GB RAM machine.
}

%\examples{
%## given a directory containing the following files:
%## ~\$> ls
%## chrom.bed  chrom.bim  chrom.dat  chrom.fam
%## the code to read the bed file into R is

%\dontrun{setwd("/path/to/the/binary/file/")}
%\dontrun{my.chrom = read.plink("chrom")}

%## alternatively, without using setwd()

%\dontrun{my.chrom = read.plink("/path/to/the/binary/file/chrom")}

%## both methods would read the data in chrom.bed into the object called my.chrom



%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ connection }
\keyword{ genetics }
\keyword{ data imput }


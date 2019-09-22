
\name{mPhen.readGenotypes}
\alias{mPhen.readGenotypes}
\title{Open, and read from a read connection to a genotype file}
\description{Opens a read connection to a list of files which have a VCF-like format. Can read a .gz file. Also supports
plink bed format. Also supports cnvPipe format.  Also supports a .zip file format used by cnvHap.
}
\usage{mPhen.readGenotypes(genoConnection, indiv = NULL,opts =mPhen.options("geno.input" ))}
\arguments{
\item{genoConnection}{
A list of paths to genotype files.  Each file might represent a different cohort with overlapping snp set. 
These files can be .vcf .vcf.gz and .zip files.  They can also be plink bed files.  In the case that bed files are used, the root name of the file should be given (i.e. hapmap2 instead of hapmap2.bed).

Alternatively, this can also be the object returned by  mPhen.readGenotypes(...).  This means that the connection only has to be established once, and can be read from multiple times in batches.
}
\item{indiv}{A list of individual ids.  If provided the results will be in this order.  If not provided the results will be in the order given by the genotype file}
\item{opts}{A list of options, which is obtained from mPhen.options("geno.input").  To get more information about these options, type 
mPhen.options("geno.input",descr=TRUE)}
}
\value{Returns a list object which can be used by subsequent calls to readGenotypes(..).  This list also includes a value genoData, which contains the genotypes which have been read.  This  includes a connection to the underlying file (conn), a list of sample ids, and a flag indicating whether the file is zipped.  It also includes 'pcs', which is a matrix of genotype pcs.  These are only calculated is opts$mPhen.numGenoPCs>0, and are only calculated once all batches have been read (but will include all genotypes which have been read previously using the same connection).}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetics }
\keyword{ GWAS }
\keyword{ regression }
\keyword{ multiple phenotypes } 
\keyword{ association }
\keyword{ file }
\keyword{ connection }

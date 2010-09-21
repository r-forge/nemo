% 
% TODO: Manual for class StaticNEM 
% 
\name{StaticNEMData}
\alias{StaticNEMData}
\title{
	StaticNEMData class
}
\description{
	This is an class defining the data format and operations for class StaticNEM.  
}
\usage{
	StaticNEMData(design=NULL,data=NULL,process_type=NULL)
}
\arguments{
	\item{design}{
		the design table including sample name, group name and replicate id
	}
	\item{data}{
		data table with column orders corresponding to the design table
	}
	\item{process_type}{
		specifies the pre-processing to compute E-gene effects
	}
}
\value{
	a StaticNEMData object
}
\details{
	This class includes pre-processing steps for both continuous (not finished) and discrete E-gene effect computing.   
}
\author{
	Xin Wang
}
\examples{
	# demo("static_nem", "nemo")
}


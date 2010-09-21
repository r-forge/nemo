% 
% TODO: Manual for class NEM 
% 
\name{NEM}
\alias{NEM}
\title{
	NEM class
}
\description{
	This is an abstract class defining a Nested Effects Model. A NEM is a graphical model for signalling pathway reconstruction.  
}
\usage{
	NEM(graph=NULL,dat=NULL,para=NULL)
}
\arguments{
	\item{graph}{
		a graph object specifying the network structure of a NEM
	}
	\item{dat}{
		a data object specifying the input data to NEM
	}
	\item{para}{
		a parameter object specifying parameters for NEM inference
	}
}
\value{
	a NEM object
}
\details{
	This is an abstract class which will be inherited by class StaticNEM.
}
\author{
	Xin Wang
}
\examples{
	# demo("static_nem", "nemo")
}


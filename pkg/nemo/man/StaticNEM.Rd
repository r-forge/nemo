% 
% TODO: Manual for class StaticNEM 
% 
\name{StaticNEM}
\alias{StaticNEM}
\title{
	StaticNEM class
}
\description{
	This is an class defining a Static Nested Effects Model for inference of static signalling pathways. It inherits from the abstract class NEM.  
}
\usage{
	StaticNEM(static_nem_para=NULL,static_nem_data=NULL)
}
\arguments{
	\item{static_nem_data}{
		a StaticNEMData object specifying the input data to StaticNEM
	}
	\item{static_nem_para}{
		a StaticNEMPara object specifying parameters for StaticNEM
	}
}
\value{
	a StaticNEM object
}
\details{
	This is the central class for inference of static signalling pathway. According to input static systematic perturbation data set and parameters, this class will 
	infer the static signalling pathway by invoking the StaticNEMBayesianModel.
}
\author{
	Xin Wang
}
\examples{
	# demo("static_nem", "nemo")
}


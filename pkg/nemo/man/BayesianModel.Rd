% 
% TODO: Manual for class BayesianModel 
% 
\name{BayesianModel}
\alias{BayesianModel}
\title{
	BayesianModel class
}
\description{
	This is an abstract class defining a Bayesian Model.   
}
\usage{
	BayesianModel(observation=NULL,model=NULL,data_type=NULL,bayesian_infer_engine=NULL,prior=NULL)
}
\arguments{
	\item{observation}{
		downstream E-gene effects after systematic perturbations
	}
	\item{model}{
		input model/models candidates
	}
	\item{data_type}{
		the data type of input observations--discrete or continuous
	}
	\item{bayesian_infer_engine}{
		maximum likelihood estimation or maximum a posteriori
	}
	\item{prior}{
		prior distribution 
	}
}
\value{
	a BayesianModel object
}
\details{
	This is an abstract class which will be inherited by class StaticNEMBayesianModel.
}
\author{
	Xin Wang
}
\examples{
	# demo("static_nem", "nemo")
}


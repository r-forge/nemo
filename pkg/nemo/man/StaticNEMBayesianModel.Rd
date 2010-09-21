% 
% TODO: Manual for class BayesianModel 
% 
\name{StaticNEMBayesianModel}
\alias{StaticNEMBayesianModel}
\title{
	StaticNEMBayesianModel class
}
\description{
	This is a class inheriting from the abstract class BayesianModel. It include detailed inference methods for static NEMs.   
}
\usage{
	StaticNEMBayesianModel(observation=NULL,model=NULL,data_type=NULL,bayesian_infer_engine=NULL,egene_pos_infer_engine=NULL,
									ss_prior=NULL,ss_prior_type=NULL,se_prior=NULL,alpha=NULL,beta=NULL,lambda=NULL)
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
	\item{ss_prior}{
		prior distribution of edges
	}
	\item{ss_prior_type}{
		the type of edge prior
	}
	\item{egene_pos_infer_engine}{
		marginalize over all E-gene positions or do MAP estimation
	}
	\item{se_prior}{
		prior information about S-gene-E-gene regulations
	}
	\item{alpha}{
		false positive error rate alpha
	}
	\item{beta}{
		false negative error rate beta
	}
	\item{lambda}{
		regularization weight
	}
}
\value{
	a StaticNEMBayesianModel object
}
\details{
	The class StaticNEMBayesianModel acts as the core inference engine for static NEMs.  
}
\author{
	Xin Wang
}
\examples{
	# demo("static_nem", "nemo")
}


% 
% TODO: Manual for class StaticNEM 
% 
\name{StaticNEMPara}
\alias{StaticNEMPara}
\title{
	StaticNEMPara class
}
\description{
	This is an class defining the data format and operations for class StaticNEM.  
}
\usage{
	StaticNEMPara(nem_infer_engine=NULL, bayesian_infer_engine=NULL, data_type=NULL, egene_pos_infer_engine=NULL, graphic_engine=NULL, alpha=NULL, beta=NULL, 
	lambda=NULL, ss_prior=NULL, ss_prior_type=NULL,	se_prior=NULL, greedy_initial=NULL, triplet_thresh=NULL, transitive_close=NULL, delta=NULL)
}
\arguments{
	\item{nem_infer_engine}{
		the inference engine for NEM
	}
	\item{bayesian_infer_engine}{
		maximum likelihood estimation or maximum a posteriori
	}
	\item{data_type}{
		the type of input data--discrete or continuous
	}
	\item{egene_pos_infer_engine}{
		marginalize over all E-gene positions or do MAP estimation
	}
	\item{graphic_engine}{
		the name of selected graphic engine to visualize
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
	\item{ss_prior}{
		prior distribution of edges
	}
	\item{ss_prior_type}{
		the type of edge prior
	}
	\item{se_prior}{
		prior information about S-gene-E-gene regulations
	}
	\item{greedy_initial}{
		the starting model of greedy inference
	}
	\item{triplet_thresh}{
		the threshold for triplet inference
	}
	\item{transitive_close}{
		transitively closed graph or not
	}
	\item{delta}{
		a parameter for adjustment of E-gene position prior when input data type is continuous
	}
}
\value{
	a StaticNEMPara object
}
\details{
	Parameters nem_infer_engine, bayesian_infer_engine, data_type, egene_pos_infer_engine and graphic_engine must be provided, because they are essential to run NEM.
	Parameters alpha and beta must be provided or generated only when data_type is "DISC".
	Parameters lambda, ss_prior, ss_prior_type and	se_prior are needed only when bayesian_infer_engine is "MAP".
	Parameter greedy_initial is needed only when nem_infer_engine is "GREEDY".
	Parameter triplet_thresh should be provided only when nem_infer_engine is "TRIPLETS".
	Parameter delta should be specified only when data_type is in c("CONT_PVAR","CONT_RATIO") and egene_pos_infer_engine is "MAP".     
}
\author{
	Xin Wang
}
\examples{
	# demo("static_nem", "nemo")
}


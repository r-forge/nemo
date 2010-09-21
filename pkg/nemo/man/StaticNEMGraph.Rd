% 
% TODO: Manual for class StaticNEM 
% 
\name{StaticNEMGraph}
\alias{StaticNEMGraph}
\title{
	StaticNEMGraph class
}
\description{
	This is an class defining the network structure and operations for NEM.  
}
\usage{
	StaticNEMGraph(s_genes=NULL,e_genes=NULL,n_sgenes=NULL,n_egenes=NULL,ss_edges=NULL,se_edges=NULL,adj=NULL, graphic_engine=NULL)
}
\arguments{
	\item{s_genes}{
		S-gene names
	}
	\item{e_genes}{
		E-gene names
	}
	\item{n_sgenes}{
		the number of S-genes
	}
	\item{n_egenes}{
		the number of E-genes
	}
	\item{ss_edges}{
		S-gene interaction table
	}
	\item{se_edges}{
		S-gene-E-gene interaction table
	}
	\item{adj}{
		adjacency matrix for S-gene Graph
	}
	\item{graphic_engine}{
		the name of selected graphic engine to visualize
	}
}
\value{
	a StaticNEMGraph object
}
\details{
	This class defines the network structure using multiple ways such as adjacency matrix and interaction tables, and supports their conversions to each other. In some
	cases when users don't know the exact names of S-genes, they can still input the number of S-genes and do some graphical operations such as graph generation.   
}
\author{
	Xin Wang
}
\examples{
	# demo("static_nem", "nemo")
}


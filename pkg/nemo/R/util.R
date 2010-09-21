###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 13:45:07, on 9 Aug 2010
###############################################################################
Subsets <- function(n, r, v = 1:n, set = TRUE) {
	if(r < 0 || r > n) stop("invalid r for this n")
	if(set) {
		v <- unique(sort(v))
		if (length(v) < n) stop("too few different elements")
	}
	v0 <- vector(mode(v), 0)
	sub <- function(n, r, v) { ## Inner workhorse
		if(r == 0) v0 else
		if(r == n) matrix(v, 1, n) else
			rbind(cbind(    v[1],
							Recall(n-1, r-1, v[-1])),
					Recall(n-1, r, v[-1]))
	}
	sub(n, r, v[1:n])
}
PhiDistr <- function(Phi, Pm, a=1, b=0.5){		
	d <- abs(Phi - Pm)
	pPhi = a/(2*b)*(1 + d/b)^(-a-1)	
	sum(log(pPhi))
}
#paste one vector to a string
CollapseArrayToString <- function(arr=NULL,sep=',') {
	if(is.null(arr))
		stop("[NEMO util]: No array has been inputed!")
	lapply(list(arr),paste,collapse=sep)[[1]]
}

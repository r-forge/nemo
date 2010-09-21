###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 14:57:07, on 19 Jul 2010
###############################################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Constructor~~~~~~~~~~~~~~~~~~~~~~~~~~~##
setConstructorS3(
		name		=	"BayesianModel", 
		definition	=	function(observation=NULL,model=NULL,data_type=NULL,bayesian_infer_engine=NULL,prior=NULL) {
			
							extend(
									Object(),
									"BayesianModel",
									.obs					=	observation,
									.model					=	model,
									.data_type				=	data_type,
									.bayesian_infer_engine	=	bayesian_infer_engine,
									.likelihood				=	NULL,
									.prior					=	prior
							)
						},
		protected	=	TRUE,
		abstract	=	TRUE
)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Methods~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
###Marginal Likelihood calculation
setMethodS3(
		name		=	"Likelihood",
		class		=	"BayesianModel",
		definition	=	function(this,...) {
							cat("[BayesianModel]: Likelihood calculation! \n")
							
						}
)
###Posterior calculation
setMethodS3(
		name		=	"Posterior",
		class		=	"BayesianModel",
		definition	=	function(this,...) {
							cat("[BayesianModel]: Posterior calculation! \n")
						}
)
###Bayesian inference
setMethodS3(
		name		=	"Infer",
		class		=	"BayesianModel",
		definition	=	function(this,...) {
							cat("[BayesianModel]: Inferece! \n")
						}
)
###Maximum likelihood inference
setMethodS3(
		name		=	"ML",
		class		=	"BayesianModel",
		definition	=	function(this,...) {
							cat("[BayesianModel]: Inferece by Maximum Likelihood Estimate! \n")
						}
)
###Maximum a posterior inference
setMethodS3(
		name		=	"MAP",
		class		=	"BayesianModel",
		definition	=	function(this,...) {
							cat("[BayesianModel]: Inference by Maximum a posteriori Estimate! \n")
						}
)


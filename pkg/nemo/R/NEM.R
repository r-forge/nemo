###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 13:16:03, on 24 Aug 2010
###############################################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Constructor~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#NEM is a abstract class which is inherited by all other NEM classes;
#it defines the framework of a Nested Effect Model: 
#member variables: a graphical structure, observations, parameters and results achieved after NEM inference
#member functions: Infer(), Show(), Save() 
setConstructorS3(
		name="NEM", 
		definition	=	function(graph=NULL,dat=NULL,para=NULL) {
							#[R.oo routine]
							extend(
									Object(),
									"NEM",
									.graph		=	graph,						#specifies the graphical structure
									.dat		=	dat,						#input data
									.para		=	para,						#parameters
									.statistics	=	NULL						#statistics got after inference
							)
						},
		protected	=	TRUE,
		abstract	=	TRUE
)
###Need to be filled up in detailed class that inherits NEM
setMethodS3(
		name		=	"Infer",
		class		=	"NEM",
		definition	=	function(this, ...) {
			cat("[NEM]: Inference! \n ")
		}
)
###Need to be filled up in detailed class that inherits NEM
setMethodS3(
		name		=	"Show",
		class		=	"NEM",
		definition	=	function(this, ...) {
			cat("[NEM]: Plot the signalling pathway inferred by NEM! \n ")
		}
)
###Need to be filled up in detailed class that inherits NEM
setMethodS3(
		name		=	"Save",
		class		=	"NEM",
		definition	=	function(this, ...) {
			cat("[NEM]: Save NEM object! \n ")
		}
)
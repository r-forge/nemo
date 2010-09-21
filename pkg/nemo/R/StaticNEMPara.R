###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 10:44:52, on 19 Jul 2010
###############################################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Constructor~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#infer_engine: 			SEARCH, GREEDY, PAIRWISE, TRIPLETS, MODULE_NETWORK
#bayesian_infer_engine:	ML, MAP, FULLML
#data_type:				DISC, CONT_PVAR, CONT_RATIO, CONT_DENS
#egene_pos_infer_engines:	MAP, NONE
setConstructorS3(
		name		=	"StaticNEMPara", 
		definition	=	function(
									nem_infer_engine=NULL, bayesian_infer_engine=NULL, data_type=NULL, egene_pos_infer_engine=NULL, graphic_engine=NULL,	#must be provided
									alpha=NULL, beta=NULL,											#needed only when data_type					=	"DISC" 
									lambda=NULL, ss_prior=NULL, ss_prior_type=NULL,					#needed only when bayesian_infer_engine		=	"MAP"
									se_prior=NULL,													#needed only when egene_pos_infer_engine	=	"MAP"
									greedy_initial=NULL,											#needed only when nem_infer_engine			=	"GREEDY"
									triplet_thresh=NULL,											#needed only when nem_infer_engine			=	"TRIPLETS"
									transitive_close=NULL,											#
									delta=NULL														#needed only when data_type %in% c("CONT_PVAR","CONT_RATIO") & egene_pos_infer_engine=="MAP"
								) {
							#R.oo routine
							extend(
									Object(),
									"StaticNEMPara",
									.nem_infer_engine		=	nem_infer_engine,
									.bayesian_infer_engine	=	bayesian_infer_engine,
									.data_type				=	data_type,
									.egene_pos_infer_engine	=	egene_pos_infer_engine,
									.alpha					=	alpha,
									.beta					=	beta,
									.lambda					=	lambda,
									.delta					=	delta,
									.ss_prior				=	ss_prior,
									.ss_prior_type			=	ss_prior_type,
									.se_prior				=	se_prior,
									.greedy_initial			=	greedy_initial,
									.triplet_thresh			=	triplet_thresh,
									.transitive_close		=	transitive_close,
									.graphic_engine			=	graphic_engine
							)
						}
)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Methods~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
###check parameters
setMethodS3(
		name		=	"CheckPara",
		class		=	"StaticNEMPara",
		definition	=	function(this, ...) {
							cat("[StaticNEMPara]: checking parameters!\n")
							#check input
							all_nem_infer_engine		<-	c("SEARCH", "GREEDY", "PAIRWISE", "TRIPLETS", "MODULE_NETWORK")
							all_bayesian_infer_engines	<-	c("ML","MAP","FULLML")
							all_data_types				<-	c("DISC", "CONT_EFFECT", "CONT_RATIO", "CONT_DENS")
							all_egene_pos_infer_engine	<-	c("MAP","BAYES","NONE")
							all_ss_prior_type			<-	c("Regularization", "Bayesian_Model_AVG")
							all_graphic_engine			<-	c("Rgraphviz","dyNet")
							
							if(is.null(this$.nem_infer_engine)) {
								stop("[StaticNEMPara]: Must specify an inference engine! \n ")
							} else if(!(this$.nem_infer_engine %in% all_nem_infer_engine)) {
								stop(paste("[StaticNEMPara]: Must input one of the following inference engines-- ", CollapseArrayToString(all_nem_infer_engine,','),sep=""))
							}	
							if(is.null(this$.bayesian_infer_engine)) {
								stop("[StaticNEMPara]: Must specify a bayesian model! \n ")
							} else if(!(this$.bayesian_infer_engine %in% all_bayesian_infer_engines)) {
								stop(paste("[StaticNEMPara]: Must input one of the following bayesian models-- ", CollapseArrayToString(all_bayesian_models,','),sep=""))
							}	
							if(is.null(this$.data_type)) {
								stop("[StaticNEMPara]: Must specify the data type! \n")
							} else if(!(this$.data_type %in% all_data_types)) {
								stop(paste("[StaticNEMPara]: Must input one of the following data types-- ", CollapseArrayToString(all_data_types,','),sep=""))
							}
							if(is.null(this$.egene_pos_infer_engine)) {
								stop("[StaticNEMPara]: Must specify an inference engine for E-gene position!\n ")
							} else if(!(this$.egene_pos_infer_engine %in% all_egene_pos_infer_engine)) {
								stop(paste("[StaticNEMPara]: Must input one of the following infer engine for E-gene positions-- ",	CollapseArrayToString(all_egene_pos_infer_engine,','),sep=""))
							}			
							if(is.null(this$.graphic_engine)) {
								stop("[StaticNEMPara]: must specify a graphic engine")
							} else if(!(this$.graphic_engine %in% all_graphic_engine)) {
								stop("[StaticNEMPara]: graphic_engine not supported! ")
							} 
								
							#check logics among inputs
							##when data_type=='DISC', alpha, beta must be provided
							if(this$.data_type=='DISC') {
								if(is.null(this$.alpha) | is.null(this$.beta))
									stop("[StaticNEMPara]: error rates alpha and beta must be provided when data_type=='DISC'! \n")
								if(this$.alpha>1 | this$.alpha <0 | this$.beta>1 | this$.beta<0)
									stop("[StaticNEMPara]: Error rates must be in (0, 1)")
							}
							##when bayesian_infer_engine=='MAP', lambda, ss_prior, ss_prior_type must be provided
							if(this$.bayesian_infer_engine=='MAP') {
								if(is.null(this$.lambda) | is.null(this$.ss_prior) | is.null(this$.ss_prior_type)) {
									stop("[StaticNEMPara]: lambda, ss_prior and ss_prior_type must be provided when bayesian_infer_engine=='MAP'! \n")
								}
								if(!(this$.ss_prior_type %in% all_ss_prior_type))
									stop(paste("[StaticNEMPara]: Must input one of the following types to incorporate edge prior--", CollapseArrayToString(all_ss_prior_type,','),sep=""))	
								if(this$.lambda<0)
									stop("[StaticNEMPara]: lambda must be > 0")
							}
							##when nem_infer_engine=='GREEDY'
							if(this$.nem_infer_engine=='GREEDY') {
								if(is.null(this$.greedy_initial)) {
									cat("[StaticNEMPara]: greedy_initial is not provided, a diag matrix is used as the initial model for greedy inference! \n")
								}
							}
							##when nem_infer_engine=='TRIPLETS'
							if(this$.nem_infer_engine=='TRIPLETS') {
								if(is.null(this$.triplet_thresh)) {
									cat("[StaticNEMPara]: triplet_thresh is not provided, so it is set as 1 by default ! \n")
									this$.triplet_thresh<-0.5
								} else {
									if(this$.triplet_thresh<0)
										stop("[StaticNEMPara]: triplet_thresh should be > 0")	
								}
							}
							##when transitive_close==NULL
							if(is.null(this$.transitive_close)) {
								cat("[StaticNEMPara]: set transitive_close as TRUE by default! \n")
								this$.transitive_close <- TRUE
							} else {
								if(class(this$.transitive_close) != "logical")
									stop("[StaticNEMPara]: transitive_close must be logical!")
							}
							
							##
							if(this$.data_type %in% c("CONT_EFFECT","CONT_RATIO") & this$.egene_pos_infer_engine=="MAP") {
								if(!is.null(this$.delta)) {
									if(this$.delta<0)
										stop("[StaticNEMPara]: delta must be > 0")
								} else {
									cat("[StaticNEMPara]: set delta as 1 by default! \n")
									this$.delta <- 1	
								}
							}
						}
)
###
setMethodS3(
		name		=	"SetPara",
		class		=	"StaticNEMPara",
		definition	=	function(this, ...) {
							cat("[StaticNEMPara]: ")
						}
)


















###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 14:42:05, on 19 Jul 2010
###############################################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Constructor~~~~~~~~~~~~~~~~~~~~~~~~~~~##
setConstructorS3(
		name		=	"StaticNEMBayesianModel", 
		definition	=	function(	observation=NULL,model=NULL,data_type=NULL,bayesian_infer_engine=NULL,egene_pos_infer_engine=NULL,
									ss_prior=NULL,ss_prior_type=NULL,se_prior=NULL,alpha=NULL,beta=NULL,lambda=NULL) {

							extend(
									BayesianModel(
											observation				=	observation,
											model					=	model,
											data_type				=	data_type,
											bayesian_infer_engine	=	bayesian_infer_engine,
											prior					=	ss_prior
									),
									"StaticNEMBayesianModel",
									.egene_pos_infer_engine		=	egene_pos_infer_engine,
									.ss_prior_type				=	ss_prior_type,
									.se_prior 					=	se_prior,
									.alpha						=	alpha,
									.beta						=	beta,
									.lambda						=	lambda,
									.score						=	NULL,
									.result						=	NULL
							)
						}
)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Methods~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
###Bayesian inference for static NEMs
setMethodS3(
		name		=	"Infer",
		class		=	"StaticNEMBayesianModel",
		definition	=	function(this,...) {
							#cat("[StaticNEMBayesianModel]: Bayesian Inference! \n")
							if(this$.bayesian_infer_engine=="ML") 
								result<-this$ML()				#here result is a list including likelihood,...
							else if(this$.bayesian_infer_engine=="MAP")
								result<-this$MAP()				#here result is a list including posterior,...
						}
)
###Maximum likelihood inference
setMethodS3(
		name		=	"ML",
		class		=	"StaticNEMBayesianModel",
		definition	=	function(this,...) {
							#cat("[StaticNEMBayesianModel]: Maximum Likelihood inferece! \n")
							this$.score<-sapply(this$.model,this$Likelihood)
							this$MaximumScore()
						}
)
###Maximum a posterior inference
setMethodS3(
		name		=	"MAP",
		class		=	"StaticNEMBayesianModel",
		definition	=	function(this,...) {
							#cat("[StaticNEMBayesianModel]: Maximum a posterior inference! \n")
							this$.score<-sapply(this$.model,this$Posterior)
							this$MaximumScore()
						}
)
###
setMethodS3(
		name		=	"MaximumScore",
		class		=	"StaticNEMBayesianModel",
		definition	=	function(this,...) {
							s_genes		<- colnames(this$.obs$D1)
							n_s_genes	<- length(s_genes)
							s       	<- unlist(this$.score["mLL",])
							ep      	<- this$.score["pos",]
							map     	<- this$.score["mappos",] 	
							LLperGene = this$.score["LLperGene",]    
							
							ppost = exp(s - max(s))
							ppost <- ppost/sum(ppost)  	# posterior model probability
							
							if(length(ppost) > 1){
								sort_ppost <- sort(ppost,decreasing=TRUE)		
								#cat("[StaticNEMBayesianModel]: (probabilities of best and second best model:", sort_ppost[1],",",sort_ppost[2],")\n")
							} 
				
							# winning model       
							winner_id <- which.max(s)
							winner <- this$.model[[winner_id]]  
							diag(winner) <- 0  
							PHI <- matrix(0,ncol=n_s_genes,nrow=n_s_genes)
							dimnames(PHI) <- list(s_genes,s_genes)  
							for(i in 1:length(this$.model))
								PHI <- PHI + this$.model[[i]]*ppost[i]   
							selected = this$.score["mappos",winner_id][[1]]  
							selected = unique(names(unlist(selected[s_genes])))  
							# output  
							this$.result <- list(graph=winner, mLL=s, ppost=ppost, avg=PHI, pos=ep[[winner_id]], mappos=map[[winner_id]], selected=selected, LLperGene=LLperGene[[winner_id]])
						}
)
###Marginal Likelihood calculation
setMethodS3(
		name		=	"Likelihood",
		class		=	"StaticNEMBayesianModel",
		definition	=	function(this,Phi,...) {
							#cat("[StaticNEMBayesianModel]: Marginal Likelihood calculation! \n")
							score<-list()
							local_prior_incorporation<-function(L,se_prior) {	
								#if(!is.null(se_prior))
								#	LP <- L + se_prior
								#else
								#	LP <- L	
								#rss <- rowSums(exp(LP-min(LP)))
								#LLperGene = log(rss) + min(LP)
								#ep <- LP - LLperGene
								#Theta = apply(ep,1,function(e) e ==max(e))
								#s  <- sum(LLperGene)
								#map = apply(Theta,1,which)	
								#list(mLL=s,pos=ep,mappos=map,LLperGene=LLperGene)
								
								if(!is.null(se_prior))
									LP <- L*se_prior
								else
									LP <- L	
								LLperGene = log(rowSums(LP))
								ep <- LP/(rowSums(LP))  	  			
								Theta = apply(ep,1,function(e) e ==max(e))
								s  <- sum(LLperGene)
								map = apply(Theta,1,which)	
								list(mLL=s,pos=ep,mappos=map,LLperGene=LLperGene)
							}
							if(this$.data_type == "DISC") {				#Discrete
								#L  <- (this$.obs$D1 %*% (1-Phi)) 	* 	log(this$.alpha) 		+ 
								#	  (this$.obs$D0 %*% (1-Phi)) 	* 	log((1-this$.alpha)) 	+ 
								#	  (this$.obs$D1 %*% Phi) 		* 	log((1-this$.beta)) 	+ 
								#	  (this$.obs$D0 %*% Phi) 		* 	log(this$.beta)	
							  	L	<-	this$.alpha^(this$.obs$D1 %*% (1-Phi)) * 
										(1-this$.alpha)^(this$.obs$D0 %*% (1-Phi)) * 
										(1-this$.beta)^(this$.obs$D1 %*% Phi) * 
										this$.beta^(this$.obs$D0 %*% Phi)
							  	score<-local_prior_incorporation(L,this$.se_prior)
							} else if(this$.data_type %in% c("CONT_EFFECT") & this$.egene_pos_infer_engine=="NONE") {#Continuous
								L <- exp(log(this$.obs$D1)%*%Phi + log((1-this$.obs$D1))%*%(1-Phi)) 
								score<-local_prior_incorporation(L,this$.se_prior)
							} else if(this$.data_type %in% c("CONT_RATIO","CONT_DENS") & this$.egene_pos_infer_engine=="BAYES") {#Continuous--integrate out E-gene positions
								L <- exp(this$.obs$D1%*%Phi)
								score<-local_prior_incorporation(L,this$.se_prior)
							} else if(this$.data_type %in% c("CONT_DENS","CONT_RATIO") & this$.egene_pos_infer_engine=="MAP") {	#Continuous--estimate E-gene positions with MAP inference
								Phi2 = cbind(Phi, double(ncol(Phi)))
								colnames(Phi2)[ncol(Phi2)] = "null"
								ep = this$.obs$D1%*%Phi2 + log(this$.se_prior)	
								Theta = apply(ep,1,function(e) e ==max(e))
								L = t((Phi2%*%(Theta*1)>0)*1)*this$.obs$D1
								LLperGene=rowSums(L)		
								s = sum(LLperGene)			
								map = apply(Theta,1,which)	
								score <- list(mLL=s,pos=ep,mappos=map,LLperGene=LLperGene)
							}
							if(!is.null(rownames(this$.obs$D1)))
								score$map = sapply(score$map, names)    
							return(score)
						}
)
###Posterior calculation
setMethodS3(
		name		=	"Posterior",
		class		=	"StaticNEMBayesianModel",
		definition	=	function(this,Phi,...) {
							#cat("[StaticNEMBayesianModel]: Posterior calculation! \n")
							score<-this$Likelihood(Phi)
							if(!is.null(this$.prior)){
								if(this$.ss_prior_type == "Regularization"){
									score$mLL<-this$Regularize(Phi,score$mLL)  
								} else if(this$.ss_prior_type == "Bayesian_Model_AVG") {
									score$mLL<-this$Bayesian_Model_AVG(Phi,score$mLL)
								}
							}
							return(score)
						}
)
###Regularization
setMethodS3(
		name		=	"Regularize",
		class		=	"StaticNEMBayesianModel",
		definition	=	function(this,Phi,likelihood,...) {
							#cat("[StaticNEMBayesianModel]: Regularization \n")		
							posterior <- likelihood- this$.lambda*sum(abs(Phi - this$.prior)) + ncol(Phi)^2*log(this$.lambda*0.5)  
							return(posterior)
						}
)
###Bayesian model averaging
setMethodS3(
		name		=	"Bayesian_Model_AVG",
		class		=	"StaticNEMBayesianModel",
		definition	=	function(this,Phi,likelihood,...) {
							#cat("[StaticNEMBayesianModel]: Bayesian model averaging \n")
							posterior <- likelihood + PhiDistr(Phi, this$.prior, a=1, b=0.5)
							return(posterior)
						}
)




















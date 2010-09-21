###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 21:56:13, on 18 Jul 2010
###############################################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Constructor~~~~~~~~~~~~~~~~~~~~~~~~~~~##
setConstructorS3(
		name		=	"StaticNEM", 
		definition	=	function(static_nem_para=NULL,static_nem_data=NULL) {
			
							if(!is.null(static_nem_para) & !is.null(static_nem_data)) {
								#find all S-genes from input data
								s_genes	<-	setdiff(unique(static_nem_data$.static_nem_design$.design[,"group"]),c("ctr","ctr_pos","ctr_neg"))
								#find all E-genes from input data
								e_genes	<-	rownames(static_nem_data$.pro_dat)
								graphic_engine	<-	static_nem_para$.graphic_engine
							} else {
								static_nem_data	<-	NULL
								static_nem_para	<-	NULL
								s_genes	<-	NULL
								e_genes	<-	NULL
								graphic_engine	<-	NULL
							}
							
							#[R.oo routine]
							extend(
									NEM(
											graph	=	StaticNEMGraph(
																s_genes		=	s_genes,
																e_genes		=	e_genes,
																n_sgenes	=	length(s_genes),
																n_egenes	=	length(e_genes),
																ss_edges	=	NULL,
																se_edges	=	NULL,
																adj			=	NULL,
																graphic_engine	=	graphic_engine),							#Initialize StaticNEMGraph
											dat		=	static_nem_data,
											para	=	static_nem_para
									),
									"StaticNEM"
							)
						},
		protected	=	TRUE
)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Methods~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
###The central member function for NEM Inference. 
###It invokes specific inference member function according to user-specified nem_infer_engine
setMethodS3(
		name		=	"Infer",
		class		=	"StaticNEM",
		definition	=	function(this,...) {
							#checking parameters
							this$.para$CheckPara()
							#cat("[StaticNEM]: Inference \n")
							if(this$.para$.nem_infer_engine=="SEARCH")
								this$InferExhaustiveSearch()
							else if(this$.para$.nem_infer_engine=="GREEDY")
								this$InferGreedy()
							else if(this$.para$.nem_infer_engine=="PAIRWISE")
								this$InferPairwise()
							else if(this$.para$.nem_infer_engine=="TRIPLETS")
								this$InferTriplets()
							else if(this$.para$.nem_infer_engine=="MODULE_NETWORK")
								this$InferModuleNetwork()
						}
)
###Summarize E-gene effects based on input discretized data
setMethodS3(
		name		=	"EgeneEffectSummarize",
		class		=	"StaticNEM",
		definition	=	function(this, s_genes_sel, ...) {
							if(!all(s_genes_sel %in% this$.graph$.s_genes))
								stop("[StaticNEM]: Input correct s_genes! \n")
							obs<-list()
							if(this$.para$.data_type == "DISC" & this$.para$.bayesian_infer_engine %in% c("ML","MAP","FULLML")){
								obs$D1 = sapply(s_genes_sel, function(s) rowSums(this$.dat$.pro_dat[,this$.dat$.static_nem_design$GetSampleByGroup(s),drop=FALSE]))  
								obs$D0 = sapply(s_genes_sel, function(s) length(this$.dat$.static_nem_design$GetSampleByGroup(s))) - obs$D1 
							}
							else{
								obs$D1 = this$.dat$.pro_dat
								obs$D0 = NULL
							}
							return(obs)
						}
)
###Set up local prior of E-gene positions
setMethodS3(
		name		=	"SetEgenePositionPrior",
		class		=	"StaticNEM",
		definition	=	function(this, e_genes_sel, s_genes_sel, ...) {
							#check input S-genes
							if(!all(s_genes_sel %in% this$.graph$.s_genes))
								stop("[StaticNEM]: Input correct S_genes! \n")
							if(!all(e_genes_sel %in% this$.graph$.e_genes))
								stop("[StaticNEM]: Input correct E_genes! \n")
							n_s_genes_sel<-length(s_genes_sel)
							n_e_genes_sel<-length(e_genes_sel)
							#set prior
							if(is.null(this$.para$.se_prior)) {
								se_prior <- matrix(1/n_s_genes_sel,nrow=n_e_genes_sel,ncol=n_s_genes_sel)
								colnames(se_prior) <- s_genes_sel
								rownames(se_prior) <- e_genes_sel
							} else {
								se_prior <- this$.para$.se_prior[e_genes_sel, s_genes_sel]	
							}          
							#adjust prior for Continuous data types
							if(this$.para$.data_type %in% c("CONT_DENS","CONT_RATIO") & this$.para$.egene_pos_infer_engine=="MAP") {
								se_prior <- cbind(se_prior, double(this$.graph$.n_egenes))
								se_prior[,ncol(se_prior)] <- this$.para$.delta/n_s_genes_sel
								se_prior <- se_prior/rowSums(se_prior)
							}
							return(se_prior) 
						}
)
###Set up edge prior
setMethodS3(
		name		=	"SetSSPrior",
		class		=	"StaticNEM",
		definition	=	function(this, s_genes_sel, ...) {
							if(this$.para$.bayesian_infer_engine=="MAP") {
								if(is.null(this$.para$.ss_prior)){
									if(this$.para$.lambda != 0) {
										cat("[StaticNEM]: Generating sparsity prior automatically! \n")
										ss_prior = diag(length(s_genes_sel))	
										dimnames(ss_prior) <- list(s_genes_sel, s_genes_sel)
									}
								} else {
									ss_prior <- this$.para$.ss_prior[s_genes_sel, s_genes_sel]
								}	
							} else {
								ss_prior<-NULL
							}
							return(ss_prior)
						}
)
#########################StaticNEM Inference Methods##########################
###Exhaustive searching
setMethodS3(
		name		=	"InferExhaustiveSearch",
		class		=	"StaticNEM",
		definition	=	function(this, ...) {
							cat("[StaticNEM]: Inference by exhaustive searching! \n")
							#preparations
							models		<-	this$.graph$EnumerateSSGraph(trans.close=TRUE)
				
							obs			<-	this$EgeneEffectSummarize(this$.graph$.s_genes)
							
							se_prior	<-	this$SetEgenePositionPrior(this$.graph$.e_genes, this$.graph$.s_genes)
							ss_prior	<-	this$SetSSPrior(this$.graph$.s_genes)
							
							#construct a Static NEM Baysian Model
							sbm<-StaticNEMBayesianModel(
									observation					=	obs,
									model						=	models,
									data_type					=	this$.para$.data_type,
									bayesian_infer_engine		=	this$.para$.bayesian_infer_engine,
									egene_pos_infer_engine		=	this$.para$.egene_pos_infer_engine,
									ss_prior					=	ss_prior,
									ss_prior_type				=	this$.para$.ss_prior_type,
									se_prior					=	se_prior,
									alpha						=	this$.para$.alpha,
									beta						=	this$.para$.beta,
									lambda						=	this$.para$.lambda)
							#Do inference
							sbm$Infer()
							#Save results
							this$.graph$.adj			<-	sbm$.result$graph
							this$.graph$.se_edges 		<-	sbm$.result$mappos
							this$.graph$.e_genes_sel 	<-	sbm$.result$selected
							this$.statistics			<-	list(score=sbm$.result$mLL, ppost=sbm$.result$ppost, 
																avg=sbm$.result$avg, opt_ep=sbm$.result$pos)
						}
)
###Greedy hill-climbing
setMethodS3(
		name		=	"InferGreedy",
		class		=	"StaticNEM",
		definition	=	function(this, ...) {
							cat("[StaticNEM]: Inference by greedy hill-climbing! \n")
							#If no specified initial graph, then generate a graph without any edge
							if(is.null(this$.para$.greedy_initial))
								Phi <- matrix(0,nrow=this$.graph$.n_sgenes,ncol=this$.graph$.n_sgenes)		
							else
								Phi = this$.para$.greedy_initial		
							diag(Phi) 		<- 	1	
							dimnames(Phi) 	<- 	list(this$.graph$.s_genes,this$.graph$.s_genes)	
							obs				<-	this$EgeneEffectSummarize(this$.graph$.s_genes)
							se_prior		<-	this$SetEgenePositionPrior(this$.graph$.e_genes, this$.graph$.s_genes)
							ss_prior		<-	this$SetSSPrior(this$.graph$.s_genes)
							#Initiation
							sbm				<-	StaticNEMBayesianModel(
													observation					=	obs,
													model						=	list(Phi),
													data_type					=	this$.para$.data_type,
													bayesian_infer_engine		=	this$.para$.bayesian_infer_engine,
													egene_pos_infer_engine		=	this$.para$.egene_pos_infer_engine,
													ss_prior					=	ss_prior,
													ss_prior_type				=	this$.para$.ss_prior_type,
													se_prior					=	se_prior,
													alpha						=	this$.para$.alpha,
													beta						=	this$.para$.beta,
													lambda						=	this$.para$.lambda)
							
							sbm$Infer()
							sco 			<- 	sbm$.result$mLL
							finished 		<- 	FALSE
							#Greedy hill-climbing
							while(!finished){
								idx = which(Phi == 0)
								if(length(idx) > 0){
									models <- list()
									for(i in 1:length(idx)){ 			# test all possible new edges
										Phinew = Phi
										Phinew[idx[i]] = 1
										Phinew = StaticNEMGraph()$TransitiveClosure(g=Phinew,loop=TRUE)	
										models[[i]] <- Phinew
									}
									models 		<- 	unique(models)
									sbm$.model	<-	models
									sbm$Infer()
									if(max(sbm$.result$mLL) > sco){
										sco 	<- 	max(sbm$.result$mLL)
										Phi 	<- 	sbm$.result$graph			
									} else 								# otherwise no improving edge could be inserted
										finished <- TRUE
								} else
									finished <- TRUE	
							}
							sbm$.model	<-	list(Phi)
							sbm$Infer()
							#Save results
							this$.graph$.adj			<-	sbm$.result$graph
							this$.graph$.se_edges 		<-	sbm$.result$mappos
							this$.graph$.e_genes_sel 	<-	sbm$.result$selected
							this$.statistics			<-	list(	score=sbm$.result$mLL, ppost=sbm$.result$ppost, 
																 	avg=sbm$.result$avg, opt_ep=sbm$.result$pos)
						}
)
###Pairwise inference
setMethodS3(
		name		=	"InferPairwise",
		class		=	"StaticNEM",
		definition	=	function(this,...) {
							#preparations
							cat("[StaticNEM]: Pairwise inference! \n")
							obs				<-	this$EgeneEffectSummarize(this$.graph$.s_genes)			
							se_prior		<-	this$SetEgenePositionPrior(this$.graph$.e_genes, this$.graph$.s_genes)
							n_test			<-	this$.graph$.n_sgenes*(this$.graph$.n_sgenes)/2
							cat("[StaticNEM]: Decomposing S-gene-S-gene graph to ", n_test ," edges! \n")
							ss_prior		<-	this$SetSSPrior(this$.graph$.s_genes)
							ss_local_prior	<-	rep(0.25,4)
							#initiation
							graph				<-	diag(this$.graph$.n_sgenes)
							dimnames(graph)		<-	list(this$.graph$.s_genes,this$.graph$.s_genes)
							scores				<-	matrix(nrow=n_test,ncol=5)
							dimnames(scores) 	<-	list(as.character(1:n_test),c("..","->","<-","<->","support")) 
							ith					<-	1
							#traversal
							pair_graph			<-	StaticNEMGraph(n_sgenes=2)					#for later use of pair model generation
							pair_sbm			<-	StaticNEMBayesianModel(
														observation					=	NULL,
														model						=	NULL,
														data_type					=	this$.para$.data_type,
														bayesian_infer_engine		=	this$.para$.bayesian_infer_engine,
														egene_pos_infer_engine		=	this$.para$.egene_pos_infer_engine,
														ss_prior					=	NULL,
														ss_prior_type				=	this$.para$.ss_prior_type,
														se_prior					=	NULL,
														alpha						=	this$.para$.alpha,
														beta						=	this$.para$.beta,
														lambda						=	this$.para$.lambda)
							#infer each pair of S-genes
							for(sr in 1:(this$.graph$.n_sgenes-1)) {
								for(dt in (sr+1):this$.graph$.n_sgenes) {
									#filter data
									x				<-	this$.graph$.s_genes[sr]
									y				<-	this$.graph$.s_genes[dt]
									sgene_sel		<-	which(colnames(obs$D1)==x | colnames(obs$D1)==y)
									egene_sel		<-	which(rowSums(obs$D1[,sgene_sel])>0)
									obs_local		<-	list(D1=obs$D1[egene_sel, sgene_sel], D0=obs$D0[egene_sel, sgene_sel])
									#local prior
									pair_sel		<-	c(sr, dt)
				
									pair_se_prior	<-	this$SetEgenePositionPrior(names(egene_sel), c(x,y))
									pair_se_prior[rowSums(pair_se_prior)==0]	<-	1e-10
									pair_ss_prior	<-	ss_prior[pair_sel, pair_sel, drop=FALSE]
									support			<-	length(egene_sel)
									#scoring 4 models per pair of S-genes
									pair_graph$.s_genes	<-	c(x,y)
									pair_models		<-	pair_graph$EnumerateSSGraph(trans.close=TRUE)
									if(support > 0) {
										#change parameters according to current pair of S-genes
										pair_sbm$.obs		<-	obs_local
										pair_sbm$.model		<-	pair_models
										pair_sbm$.prior		<-	pair_ss_prior
										pair_sbm$.se_prior	<-	pair_se_prior
										#do inference for current pair of S-genes
										pair_sbm$Infer()
										#posterior
										pair_post					<-	exp(pair_sbm$.result$mLL) * ss_local_prior
										pair_post					<-	pair_post/sum(pair_post)
										pair_post[is.na(pair_post)]	<-	0
									} else {
										pair_post	<-	matrix(0, nrow=1, ncol=4)
									}
									#identify pair model with highest score
									pair_winner		<-	pair_models[[which.max(pair_post)]]
									graph[sr, dt]	<-	pair_winner[1, 2]
									graph[dt, sr]	<-	pair_winner[2, 1]
									#save score
									scores[ith, ]			<-	c(pair_post, support)
									rownames(scores)[ith]	<-	paste(c(x, y), collapse="~")
									ith						<-	ith+1
								}
							}
							#estimate effect positions
							cat("[StaticNEM]: Estimating effect positions \n")
							if(this$.para$.transitive_close)
								graph	<-	StaticNEMGraph()$TransitiveClosure(g=graph)
							sbm		<-	StaticNEMBayesianModel(
											observation					=	obs,
											model						=	list(graph),
											data_type					=	this$.para$.data_type,
											bayesian_infer_engine		=	this$.para$.bayesian_infer_engine,
											egene_pos_infer_engine		=	this$.para$.egene_pos_infer_engine,
											ss_prior					=	ss_prior,
											ss_prior_type				=	this$.para$.ss_prior_type,
											se_prior					=	se_prior,
											alpha						=	this$.para$.alpha,
											beta						=	this$.para$.beta,
											lambda						=	this$.para$.lambda)
							sbm$Infer()
							#Save results
							this$.graph$.adj			<-	sbm$.result$graph
							this$.graph$.se_edges 		<-	sbm$.result$mappos
							this$.graph$.e_genes_sel 	<-	sbm$.result$selected
							this$.statistics			<-	list(score=sbm$.result$mLL, ppost=sbm$.result$ppost, 
									avg=sbm$.result$avg, opt_ep=sbm$.result$pos)
						}
)
###Triplets inference
setMethodS3(
		name		=	"InferTriplets",
		class		=	"StaticNEM",
		definition	=	function(this,...) {
							cat("[StaticNEM]: Triplets inference! \n")
							n_test			<-	choose(this$.graph$.n_sgenes, 3)
							cat("[StaticNEM]: Decomposing S-gene-S-gene graph to ", n_test ," triplets! \n")
							obs				<-	this$EgeneEffectSummarize(this$.graph$.s_genes)
							se_prior		<-	this$SetEgenePositionPrior(this$.graph$.e_genes, this$.graph$.s_genes)
							ss_prior		<-	this$SetSSPrior(this$.graph$.s_genes)
							ss_local_prior	<-	rep(1/29, 29)
							#enumerate all (ordered) triples for the given S-genes
							triples			<-	Subsets(this$.graph$.n_sgenes, 3)
							#29 models for each triple of S-genes
							local_models		<-	StaticNEMGraph(s_genes=c('a', 'b', 'c'), n_sgenes=3)$EnumerateSSGraph(trans.close=TRUE)
							local_opt_models	<- 	list()					#to store optimal model for each triple of S-genes
							
							tri_sbm			<-	StaticNEMBayesianModel(
									observation					=	NULL,
									model						=	local_models,
									data_type					=	this$.para$.data_type,
									bayesian_infer_engine		=	this$.para$.bayesian_infer_engine,
									egene_pos_infer_engine		=	this$.para$.egene_pos_infer_engine,
									ss_prior					=	NULL,
									ss_prior_type				=	this$.para$.ss_prior_type,
									se_prior					=	NULL,
									alpha						=	this$.para$.alpha,
									beta						=	this$.para$.beta,
									lambda						=	this$.para$.lambda)
							for(i_tri in 1:nrow(triples)){
								tri_sel				<-	which(colnames(obs$D1) %in% this$.graph$.s_genes[triples[i_tri,]])
								tri_obs				<-	list(D1=obs$D1[, tri_sel], D0=obs$D0[, tri_sel])
								colnames(tri_obs$D1)[colnames(tri_obs$D1) == this$.graph$.s_genes[triples[i_tri,1]]] = "a"
								colnames(tri_obs$D1)[colnames(tri_obs$D1) == this$.graph$.s_genes[triples[i_tri,2]]] = "b"
								colnames(tri_obs$D1)[colnames(tri_obs$D1) == this$.graph$.s_genes[triples[i_tri,3]]] = "c"	
								if(!is.null(tri_obs$D0))
									colnames(tri_obs$D0)<-colnames(tri_obs$D1)
								tri_se_prior		<-	this$SetEgenePositionPrior(this$.graph$.e_genes, this$.graph$.s_genes[triples[i_tri,]])
								tri_ss_prior		<-	ss_prior[tri_sel, tri_sel, drop=FALSE]
								
								tri_sbm$.obs		<-	tri_obs
								tri_sbm$.prior		<-	tri_ss_prior
								tri_sbm$.se_prior	<-	tri_se_prior
								
								tri_sbm$Infer()
								#prior
								tri_post			<-	tri_sbm$.result$mLL	+ log(ss_local_prior)
								tri_winner			<-	local_models[[which.max(tri_post)]]
								local_opt_models[[i_tri]]					<- 	list(graph = tri_winner, posterior = tri_post)
								dimnames(local_opt_models[[i_tri]]$graph) 	<- 	list(this$.graph$.s_genes[triples[i_tri,]], this$.graph$.s_genes[triples[i_tri,]])
							}
							
							#
							A			<-	matrix(NA, ncol = this$.graph$.n_sgenes, nrow = this$.graph$.n_sgenes)
							dimnames(A)	<-	list(this$.graph$.s_genes, this$.graph$.s_genes)
							diag(A)		<-	1
							for(sr in 1:(this$.graph$.n_sgenes-1)) {
								for(dt in (sr+1):this$.graph$.n_sgenes) {
									##=== find contributing triples
									inds = apply(triples,1,function(x) sum( x %in% c(sr,dt)) == 2)
									inds = which(inds)
									##=== count edges and set A_ij accordingly
									tmp = list()
									for(k in 1:length(inds)){ ##=== for each contributing triple
										tmp[[k]] = local_opt_models[[inds[k]]]$graph #=== get adjacency matrix
									}
									##=== mean number of edges from i to j (including doubles...)
									A[sr,dt] = mean(unlist(lapply(tmp,function(x) x[this$.graph$.s_genes[sr],this$.graph$.s_genes[dt]])))
									##=== mean number of edges from j to i (including doubles...) 
									A[dt,sr] = mean(unlist(lapply(tmp,function(x) x[this$.graph$.s_genes[dt],this$.graph$.s_genes[sr]])))
								}
							}
							B	=	(A>=this$.para$.triplet_thresh)*1 - diag(nrow(A))
							## 3. estimate effect positions
							if(this$.para$.transitive_close)
								B = StaticNEMGraph()$TransitiveClosure(B)
							diag(B) = 1								#?
							sbm		<-	StaticNEMBayesianModel(
									observation					=	obs,
									model						=	list(B),
									data_type					=	this$.para$.data_type,
									bayesian_infer_engine		=	this$.para$.bayesian_infer_engine,
									egene_pos_infer_engine		=	this$.para$.egene_pos_infer_engine,
									ss_prior					=	ss_prior,
									ss_prior_type				=	this$.para$.ss_prior_type,
									se_prior					=	se_prior,
									alpha						=	this$.para$.alpha,
									beta						=	this$.para$.beta,
									lambda						=	this$.para$.lambda)
							sbm$Infer()
							#Save results
							this$.graph$.adj			<-	sbm$.result$graph
							this$.graph$.se_edges 		<-	sbm$.result$mappos
							this$.graph$.e_genes_sel 	<-	sbm$.result$selected
							this$.statistics			<-	list(score=sbm$.result$mLL, ppost=sbm$.result$ppost, 
									avg=sbm$.result$avg, opt_ep=sbm$.result$pos)
						}
)
###!Module Network
setMethodS3(
		name		=	"InferModuleNetwork",
		class		=	"StaticNEM",
		definition	=	function(this,...) {
							cat("[StaticNEM]: Inference by Module Network! \n")			
						}
)
#######################Visualization and File operations######################
###The central member function for Static NEM Visualization
setMethodS3(
		name		=	"Show",
		class		=	"StaticNEM",
		definition	=	function(this, what="SSGraph", save_filename=NULL, save_pdf=FALSE, ...) {
							#check input
							if(class(save_pdf)!="logical")
								stop("[StaticNEM]: save_pdf must be logical!")
							if(save_pdf)
								if(is.null(save_filename))
									stop("[StaticNEM]: save_filename needed when save_pdf=TRUE!")
								else {
									pdf(file=save_filename)
								}		
							#check if statistics have been generated
							if(is.null(this$.statistics))
								stop("[StaticNEM]: inference has not been done or was not successful!")
							#invoke specific functions
							if(what=="SSGraph")
								this$.graph$ShowSSGraph()
							else if(what=="FULLGraph")
								this$.graph$ShowFULLGraph()
							else if(what=="ScoreDistribution")
								this$ShowScoreDistribution()
							else if(what=="EgenePosPosterior")
								this$ShowEgenePosPost()
							#close dev
							if(save_pdf) dev.off()
						}
)
###Show the scoring order of top models
setMethodS3(
		name		=	"ShowScoreDistribution",
		class		=	"StaticNEM",
		definition	=	function(this, ...) {
							par(cex=1.3)
							ss <- sort(unique(this$.statistics$score),decreasing=TRUE)[1:min(30,length(this$.statistics$score))]
							if(this$.para$.bayesian_infer_engine=="MAP")
								ylab <- "Log-posterior"
							else 
								ylab <- "Marginal log-likelihood"
							plot(x=1:length(ss), y=ss, pch=19, main="Score distribution",
									xlab=paste(length(ss),"top ranked models"),
									ylab=ylab, 
									ylim=c(ss[length(ss)]-10,ss[1]+10)
							)
							points(1,max(unique(this$.statistics$score)),pch=21,cex=1.7,lwd=2)
						}
)
###Show a heatmap of E-gene posterior
setMethodS3(
		name		=	"ShowEgenePosPost",
		class		=	"StaticNEM",
		definition	=	function(this, ...) {
							egenes<-rownames(this$.statistics$opt_ep)
							this$.statistics$opt_ep[is.na(this$.statistics$opt_ep)] = 0
							par(las=2,mgp=c(5.5,1,0),mar=c(6.7,7,4,1),cex.lab=1.3,cex.main=1.7)        	
							image(x=1:ncol(this$.statistics$opt_ep),
									y=1:nrow(this$.statistics$opt_ep),
									z = t(this$.statistics$opt_ep),
									main = "Posterior of E-gene positions",
									xlab="Perturbed signalling gene",
									xaxt="n",
									ylab="Effect reporter gene",
									yaxt="n",
									col=gray(seq(.95,0,length=10))
							)
							abline(v=(1:(ncol(this$.statistics$opt_ep)-1))+.5)
							axis(1,1:ncol(this$.statistics$opt_ep),colnames(this$.statistics$opt_ep))        
							axis(2,1:length(egenes),egenes)
						}
)

###
setMethodS3(
		name		=	"Save",
		class		=	"StaticNEM",
		definition	=	function(this, save_path,save_filename, ...) {
							cat("[StaticNEM]: Saving inference results! \n")
							R.utils::saveObject(object=this, file=save_filename,path=save_path)
						}
)

#######################Get members in this class######################

###
setMethodS3(
		name		=	"GetGraph",
		class		=	"StaticNEM",
		definition	=	function(this, ...) {
							return(this$.graph)
						}
)
###
setMethodS3(
		name		=	"GetData",
		class		=	"StaticNEM",
		definition	=	function(this, ...) {
							return(this$.dat)
						}
)
###
setMethodS3(
		name		=	"GetPara",
		class		=	"StaticNEM",
		definition	=	function(this, ...) {
							return(this$.para)
						}
)
###
setMethodS3(
		name		=	"GetStatistics",
		class		=	"StaticNEM",
		definition	=	function(this, ...) {
							return(this$.statistics)
						}
)



































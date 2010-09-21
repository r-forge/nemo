###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 20:26:24, on 17 Jul 2010
###############################################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Constructor~~~~~~~~~~~~~~~~~~~~~~~~~~~##
setConstructorS3(
		name		=	"StaticNEMGraph", 
		definition	=	function(s_genes=NULL,e_genes=NULL,n_sgenes=NULL,n_egenes=NULL,ss_edges=NULL,se_edges=NULL,adj=NULL, graphic_engine=NULL) {
							#check input
							if(!is.null(se_edges) & !is.null(ss_edges) & is.null(adj)) {							#If input ss_edges and se_edges
								if(class(ss_edges)!="data.frame" | class(se_edges)!="data.frame")
									stop("[StaticNEMGraph]: ss_edges and se_edges must be data frames!\n")
								t.sgenes<-unique(as.character(se_edges[,"s_source"]))
								t.egenes<-unique(as.character(se_edges[,"e_dest"]))
								if(!is.null(s_genes) & !(length(setdiff(s_genes,t.sgenes))==0 & length(setdiff(t.sgenes,s_genes))>=0)) 
									stop("[StaticNEMGraph]: s_genes and se_edges are not corresponding with each other! \n")
								if(!is.null(e_genes) & !(length(setdiff(e_genes,t.egenes))==0 & length(setdiff(t.egenes,e_genes))>=0)) 
									stop("[StaticNEMGraph]: e_genes and se_edges are not corresponding with each other! \n")	
							} else if(is.null(ss_edges) & is.null(se_edges) & !is.null(adj)) {						#If input adjacent matrix
								if(class(adj)!="matrix")
									stop("[StaticNEMGraph]: input adj must be a matrix! \n")
								adj.cols<-colnames(adj)
								adj.rows<-rownames(adj)
								if(is.null(adj.cols) | is.null(adj.rows) | !setequal(adj.cols,adj.rows))
									stop("[StaticNEMGraph]: check colnames and rownames of input adjacent matrix! \n")
							} else if(is.null(ss_edges) & is.null(se_edges) & is.null(adj) & !is.null(n_sgenes)) {	#If input the number of S-genes
								if(is.null(s_genes)) {
									s_genes<-paste("S",1:n_sgenes,sep="")
								} else {
									if(length(s_genes)!=n_sgenes) 
										stop("StaticNEMGraph: check input n_sgenes and s_genes! \n")
								}
							} 
							#[R.oo routine]
							extend(
									Object(),
									"StaticNEMGraph",
									.s_genes		=	s_genes,
									.e_genes		=	e_genes,
									.e_genes_sel	=	NULL,
									.n_sgenes		=	n_sgenes,
									.n_egenes		=	n_egenes,
									.ss_edges		=	ss_edges,
									.se_edges		=	se_edges,
									.adj			=	adj,
									.graphic_engine =	graphic_engine
							)
						}
)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Methods~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
###Transfer between graph formats
setMethodS3(
		name		=	"Tab2Adj",
		class		=	"StaticNEMGraph",
		definition	=	function(this,...) {
			cat("[StaticNEMGraph]: Converting from a table to an adjacent matrix! \n")
		}
)
setMethodS3(
		name		=	"Adj2Tab",
		class		=	"StaticNEMGraph",
		definition	=	function(this,...) {
							cat("[StaticNEMGraph]: Converting from an adjacency matrix to a table! \n")
							if(is.null(this$.adj))
								stop("[StaticNEMGraph]: Adjacency matrix not available!")
							coln<-colnames(this$.adj)
							rown<-rownames(this$.adj)
							if(is.null(coln) | is.null(rown))
								stop("[StaticNEMGraph]: No column names or row names for adjacency matrix!")
							lr<-NULL
							lc<-NULL
							
							for(r in rown) {
								for(c in coln) {
									if(r!=c & this$.adj[r,c]!=0) {
										lr<-c(lr, r)
										lc<-c(lc, c)
									}
								}
							}
							this$.ss_edges<-data.frame(src=lr, dest=lc)
						}
)
###Enumerate graphs for a given number of S-genes
setMethodS3(
		name		=	"EnumerateSSGraph",
		class		=	"StaticNEMGraph",
		definition	=	function(this, trans.close=TRUE, ...) {
							#cat("StaticNEMGraph: Enumerate S-S graph! \n")
							# Sanity checks    	
							if (this$.n_sgenes==1) stop("StaticNEMGraph: choose n>1!")
							if (this$.n_sgenes>5)  stop("StaticNEMGraph: exhaustive enumeration not feasible with more than 5 perturbed genes")            
							if (this$.n_sgenes==5) cat ("StaticNEMGraph: this will take a while ... \n") 
							bc <- bincombinations(this$.n_sgenes*(this$.n_sgenes-1))
							fkt1 <- function(x,n,name) {
								M <- diag(n)
								M[which(M==0)]<-x
								dimnames(M) <- list(name,name)	
								if(trans.close)    
									M <- this$TransitiveClosure(M,mat=TRUE,loops=TRUE)    
								return(list(M))
							}
							
							models <- apply(bc,1,fkt1,this$.n_sgenes,this$.s_genes) 
							models <- unique(matrix(unlist(models),ncol=this$.n_sgenes*this$.n_sgenes,byrow=TRUE))
							
							fkt2 <- function(x,n,name){
								M <- matrix(x,n)
								dimnames(M) <- list(name,name)
								return(list(M))
							}
							models <- unlist(apply(models,1,fkt2,this$.n_sgenes,this$.s_genes),recursive=FALSE)
							cat("Generated",length(models),"unique models ( out of", 2^(this$.n_sgenes*(this$.n_sgenes-1)), ")\n")
							return(models)
						}
)
###Find the transitive closure of a given graph
setMethodS3(
		name		=	"TransitiveClosure",
		class		=	"StaticNEMGraph",
		definition	=	function(this,g,loops=TRUE,...) {
							#cat("StaticNEMGraph: Find the transitive closure of a given graph! \n")
							if (!(class(g)=="matrix")) 
								stop("[StaticNEMGraph]: Input must be either graphNEL object or adjacency matrix! \n")
							n <- ncol(g)
							matExpIterativ <- function(x,pow,y=x,z=x,i=1) {
								while(i < pow) {
									z <- z %*% x
									y <- y+z
									i <- i+1
								}
								return(y)
							}
							h <- matExpIterativ(g,n)
							h <- (h>0)*1   
							dimnames(h) <- dimnames(g)
							if (!loops) diag(h) <- rep(0,n) else diag(h) <- rep(1,n)
							return(h)
						}
)
###Plot this graph
setMethodS3(
		name		=	"ShowSSGraph",
		class		=	"StaticNEMGraph",
		definition	=	function(this,dynet.w=600, dynet.h=400, ...) {
			
							if(is.null(this$.adj))
								stop("[StaticNEMGraph]: no graph created! \n")
							cat("[StaticNEMGraph]: Show this graph! \n")
							if(this$.graphic_engine=="dyNet") {
								if(is.null(this$.adj))
									stop("[StaticNEMGraph]: no graph created!")
								if(is.null(this$.ss_edges))
									this$Adj2Tab()
								#build big network first, and set attributes specifically for different parts
								##node
								node.s_gene.label	<-	this$.s_genes
								
								node	<-	data.frame(
										label	=	node.s_gene.label,
										id		=	1:length(node.s_gene.label),
										alias	=	node.s_gene.label
								)
								##link
								link.ss.label	<-	paste("l",1:nrow(this$.ss_edges),sep="")
								link.ss	<-	data.frame(
										label	=	link.ss.label,
										source	=	match(this$.ss_edges[,1],this$.s_genes),
										target	=	match(this$.ss_edges[,2],this$.s_genes)
								)
								link	<-	link.ss
								##create RDyNetGraph 
								rdynet.graph	<-	RDyNetGraph(node, link)
								dynet.nem.SS_region		<-	list(w=dynet.w/3,h=dynet.h/4,x=dynet.w/2,y=dynet.h/8)
								##set up specific attributes
								###node
								rdynet.graph$.node$SetNodeAttr(
										label=node.s_gene.label,
										value=c(node.size="20.0",g.h="20.0",g.w="20.0",g.fill="#ffffff",g.width=3,g.outline="#ff0000",container="false",gradient="false")
								)
								###link
								rdynet.graph$.link$SetLinkAttr(
										label=link.ss.label,
										value=c(direction="1",g.width="3.0",g.fill="#eb9911")
								)
								####layout S-S graph
								rdynet.graph$Layout(label=node.s_gene.label, constraint=dynet.nem.SS_region, layout="spring")				
								rdynet<-RDyNet(rdynet_graph=rdynet.graph, save_filename="temp.S-gene-graph.xgmml")
								rdynet$GenXML()
								rdynet$InvokeDyNet()
							} else if(this$.graphic_engine=="Rgraphviz") {
								library(Rgraphviz)
								gR <- new("graphAM",adjMat=this$.adj,edgemode="directed")  
								gR <- as(gR,"graphNEL")
								Rgraphviz::plot(gR)	
							}	
						}
)
###Future extension
setMethodS3(
		name		=	"ShowFULLGraph",
		class		=	"StaticNEMGraph",
		definition	=	function(this, dynet.w=600, dynet.h=400, ...) {
							if(is.null(this$.adj))
								stop("[StaticNEMGraph]: no graph created!")
							cat("[StaticNEMGraph]: Show this graph! \n")
							if(this$.graphic_engine!="dyNet") 
								stop("[StaticNEMGraph]: only implemented in dyNet!")	
							if(is.null(this$.se_edges))
								stop("[StaticNEMGraph]: no S-E graph created!")
							if(is.null(this$.ss_edges))
								this$Adj2Tab()
							#build big network first, and set attributes specifically for different parts
							##node
							node.s_gene.label	<-	this$.s_genes
							node.e_gene.label	<-	unique(unlist(sapply(1:length(this$.se_edges),function(l) {names(this$.se_edges[[l]])})))
							node.e_gene_container.label	<-	paste("E genes regulated by ",this$.s_genes,sep="")
							if(length(intersect(node.s_gene.label,node.e_gene.label))>0)
								stop("[StaticNEMGraph]: E-genes and S-genes have overlapped names! \n ")
							node	<-	data.frame(
													label	=	c(node.s_gene.label, node.e_gene_container.label, node.e_gene.label),
													id		=	1:(length(node.s_gene.label) + length(node.e_gene_container.label) + length(node.e_gene.label)),
													alias	=	c(node.s_gene.label, node.e_gene_container.label, node.e_gene.label)
										)
							##link
							link.ss.label	<-	paste("l",1:nrow(this$.ss_edges),sep="")
							link.se.label	<-	paste("se",1:length(unlist(this$.se_edges)),sep="")
							link.ct.label	<-	paste("ct",1:length(unlist(this$.se_edges)),sep="")
							link.ss	<-	data.frame(
											label	=	link.ss.label,
											source	=	match(this$.ss_edges[,1],this$.s_genes),
											target	=	match(this$.ss_edges[,2],this$.s_genes)
										)
							SE_source	<-	NULL
							CT_source	<-	NULL
							SE_target	<-	NULL
							for(s in 1:length(this$.s_genes)) {
								if(!is.null(this$.se_edges[[s]])) {
									SE_target	<-	c(SE_target,names(this$.se_edges[[s]]))
									CT_source	<-	c(CT_source,rep(paste("E genes regulated by ",names(this$.se_edges)[s],sep=""),length(this$.se_edges[[s]])))
									SE_source	<-	c(SE_source,rep(names(this$.se_edges)[s],length(this$.se_edges[[s]])))
								}
								
							}			
							link.se	<-	data.frame(
											label	=	link.se.label,
											source	=	match(SE_source,this$.s_genes),
											target	=	match(SE_target,node.e_gene.label)+2*this$.n_sgenes
										)
							link.ct	<-	data.frame(
											label	=	link.ct.label,
											source	=	match(CT_source, node.e_gene_container.label)+this$.n_sgenes,
											target	=	match(SE_target,node.e_gene.label)+2*this$.n_sgenes
										)
							link	<-	rbind(link.ss, link.se, link.ct)
							
							##create RDyNetGraph 
							rdynet.graph	<-	RDyNetGraph(node, link)
							##partition of NEM regions
							###intervals
							dynet.nem.par.SC_h_interval	<-	200		#vertical interval between S-gene and container
							dynet.nem.par.CC_w_interval	<-	40		#horizontal interval between adjacent containers
							dynet.nem.par.CC_h_interval	<-	40
							###generate x, y, h, w for SS-graph and containers
							
							dynet.nem.CT_size		<-	200
							dynet.nem.CT_regions	<-	list()
							CT.x<-NULL
							CT.y<-NULL
							for(c in 1:length(node.e_gene_container.label)) {
								rowi<-floor((c-1)/4)
								colj<-(c-1)%%4+1
								dynet.nem.CT_regions[[node.e_gene_container.label[c]]]	<-	
										list(w=dynet.nem.CT_size, h=dynet.nem.CT_size, x=colj*dynet.nem.par.CC_w_interval+(colj-1)*dynet.nem.CT_size+(1/2)*dynet.nem.CT_size,
											 y=dynet.nem.par.SC_h_interval+rowi*dynet.nem.par.CC_h_interval+(rowi)*dynet.nem.CT_size+(1/2)*dynet.nem.CT_size)
							 	CT.x<-c(CT.x,dynet.nem.CT_regions[[node.e_gene_container.label[c]]]$x)
								CT.y<-c(CT.y,dynet.nem.CT_regions[[node.e_gene_container.label[c]]]$y)
							}
							
							dynet.nem.SS_region		<-	list(w=dynet.w/3,h=dynet.h/4,x=dynet.w/2,y=dynet.h/8)
							##set up specific attributes
							###node
							rdynet.graph$.node$SetNodeAttr(
									label=node.s_gene.label,
									value=c(node.size="20.0",g.h="20.0",g.w="20.0",g.fill="#ffffff",g.width=3,g.outline="#ff0000",container="false",gradient="false")
							)
							
							node.e_gene_container.attr<-data.frame(
									node.size=rep(as.character(dynet.nem.CT_size),length(node.e_gene_container.label)),
									g.h=rep(dynet.nem.CT_size,length(node.e_gene_container.label)),
									g.w=rep(dynet.nem.CT_size,length(node.e_gene_container.label)),
									g.x=CT.x,
									g.y=CT.y,
									g.fill=rep("#ccccff",length(node.e_gene_container.label)),
									g.width=rep(1,length(node.e_gene_container.label)),
									g.outline=rep("#ccccff",length(node.e_gene_container.label)),
									container=rep("true",length(node.e_gene_container.label)),
									gradient=rep("true",length(node.e_gene_container.label)),
									stringsAsFactors=FALSE)
							rownames(node.e_gene_container.attr)<-node.e_gene_container.label
							rdynet.graph$.node$SetNodeAttr(
									value=node.e_gene_container.attr
							)
							rdynet.graph$.node$SetNodeAttr(
									label=node.e_gene.label,
									value=c(node.size="10",g.h="10",g.w="10",
											g.fill="#4143f8",g.width=1,g.outline="#4143f8",container="false",gradient="false")		
							)
							###link
							rdynet.graph$.link$SetLinkAttr(
									label=link.ss.label,
									value=c(direction="1",g.width="3.0",g.fill="#eb9911")
									)
							rdynet.graph$.link$SetLinkAttr(
									label=link.se.label,
									value=c(direction="1",g.width="1.0",g.fill="#88b1f7")
							)
							rdynet.graph$.link$SetLinkAttr(
									label=link.ct.label,
									value=c(direction="0",g.width="1.0",g.fill="#898883")
							)
							###layout each part
							####layout S-S graph
							rdynet.graph$Layout(label=node.s_gene.label, constraint=dynet.nem.SS_region, layout="spring")
							
							####layout E-gene cluster in each container
							egene_cluster.label<-NULL
							for(c in 1:this$.n_sgenes) {
								temp.egene_cluster.label	<-	setdiff(names(this$.se_edges[[this$.s_genes[c]]]),egene_cluster.label)
								egene_cluster.label			<-	c(egene_cluster.label, temp.egene_cluster.label)
								dynet.nem.CT_regions[[node.e_gene_container.label[c]]]$w<-dynet.nem.CT_regions[[node.e_gene_container.label[c]]]$w*0.6
								dynet.nem.CT_regions[[node.e_gene_container.label[c]]]$h<-dynet.nem.CT_regions[[node.e_gene_container.label[c]]]$h*0.6
								rdynet.graph$Layout(label=temp.egene_cluster.label, constraint=dynet.nem.CT_regions[[node.e_gene_container.label[c]]], layout="random")
							}
							rdynet<-RDyNet(rdynet_graph=rdynet.graph, save_filename="temp.fullgraph.xgmml")
							rdynet$GenXML()
							rdynet$InvokeDyNet()
						}
)

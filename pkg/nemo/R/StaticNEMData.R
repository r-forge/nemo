###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 22:20:09, on 18 Jul 2010
###############################################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Constructor~~~~~~~~~~~~~~~~~~~~~~~~~~~##
setConstructorS3(
		name		=	"StaticNEMData", 
		definition	=	function(design=NULL,data=NULL,process_type=NULL) {
			
							if(!is.null(design) & !is.null(data)) {
								if (! all(colnames(data)==design[,"sample"])) stop("StaticNEMData: sample names and column names of data must be the same! \n")	
							} else {
								design 	<- 	NULL
								data	<-	NULL
							}
							
							extend(
									Object(),
									"StaticNEMData",
									.static_nem_design	=	StaticNEMDesign(design=design),
									.pre_dat			=	data,
									.pro_dat			=	NULL,
									.alpha				=	NULL,
									.beta				=	NULL
							)
						}
)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Methods~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
###Profile including p-values describing significance of differential expression
setMethodS3(
		name		=	"DiffProfile",
		class		=	"StaticNEMData",
		definition	=	function(this, ...) {
			
						}
)
###Discretization
setMethodS3(
		name		=	"Discretize",
		class		=	"StaticNEMData",
		definition	=	function(this, n_fold=2, cutoff=0:10/10, p_counts=20, empirical_pvar=0.05, ...) {
							cat("[StaticNEMData]: Discretization!\n")
							if("ctr" %in% this$.static_nem_design$.design[,"group"]) {
								cat("[StaticNEMData]: This is one-control dataset. Applying one-control discretization! \n")
								#Get control and sample
								ctr		<-	this$.pre_dat[,as.character(this$.static_nem_design$.design[which(as.character(this$.static_nem_design$.design[,"group"])=="ctr"),"sample"])]
								samp	<-	this$.pre_dat[,as.character(this$.static_nem_design$.design[which(as.character(this$.static_nem_design$.design[,"group"])!="ctr"),"sample"])]
								# empirical distr. function
								ctr_ecdf 	<- 	apply(ctr,1,ecdf) 
								samp_pvar 	<- 	matrix(0,ncol=ncol(samp),nrow=nrow(samp))   
								for (i in 1:nrow(ctr)){
									Pi 				<- 	ctr_ecdf[[i]](samp[i,])    
									samp_pvar[i,] 	<- 	ifelse(Pi<=.5,Pi,1-Pi)   
								} 
								this$.pro_dat 			<- 	(samp_pvar <= empirical_pvar)*1   
								dimnames(this$.pro_dat) <- 	dimnames(samp)
							} else if(all(c("ctr_neg","ctr_pos") %in% this$.static_nem_design$.design[,"group"])) {
								cat("[StaticNEMData]: This is two-control dataset. Applying two-control discretization! \n")
								# select diff - maybe as extra input function??
								ctr_pos		<-	this$.pre_dat[,as.character(this$.static_nem_design$.design[which(as.character(this$.static_nem_design$.design[,"group"])=="ctr_pos"),"sample"])]
								ctr_neg		<-	this$.pre_dat[,as.character(this$.static_nem_design$.design[which(as.character(this$.static_nem_design$.design[,"group"])=="ctr_neg"),"sample"])]
								samp		<-	this$.pre_dat[,as.character(this$.static_nem_design$.design[
														which(as.character(this$.static_nem_design$.design[,"group"])!="ctr_pos" & 
																		as.character(this$.static_nem_design$.design[,"group"])!="ctr_neg"), "sample"])]
								sel 		<- 	which(rowMeans(ctr_pos) - rowMeans(ctr_neg) > log(n_fold))
								samp.sel 	<- 	samp[sel,]
								ctr_pos.sel <- 	ctr_pos[sel,]
								ctr_neg.sel <- 	ctr_neg[sel,]
								# count false decisions for different cutoff levels
								count.false.decisions <- function(x){
									thrsh    		<- 	x*rowMeans(ctr_pos.sel) + (1-x)*rowMeans(ctr_neg.sel)
									ctr_pos.disc 	<- 	(ctr_pos.sel <= thrsh)*1
									ctr_neg.disc 	<- 	(ctr_neg.sel <= thrsh)*1        
									alpha    		<- 	round(sum(ctr_pos.disc)/length(ctr_pos.disc),2)
									beta     		<- 	1-round(sum(ctr_neg.disc)/length(ctr_neg.disc),2)                   
									return(c(alpha,beta))
								}
								false 			<- 	sapply(cutoff,count.false.decisions)
								dimnames(false) <- 	list(c("alpha","beta"),as.character(cutoff))
								mycutoff 		<- 	cutoff[which.min(false[2,])]
								# apply chosen cutoff
								thrsh  			<- 	mycutoff*rowMeans(ctr_pos.sel) + (1-mycutoff)*rowMeans(ctr_neg.sel)
								this$.pro_dat 	<- 	(samp.sel <= thrsh)*1
								ctr_pos.disc 	<- 	(ctr_pos.sel <= thrsh)*1
								ctr_neg.disc 	<- 	(ctr_neg.sel <= thrsh)*1        
								this$.alpha 	<-  round((sum(ctr_pos.disc)+p_counts)/(length(ctr_pos.disc)+p_counts),2)
								this$.beta  	<-  1-round( sum(ctr_neg.disc)        /(length(ctr_neg.disc)+p_counts),2)  
				
				#				# also includes downregulation ...
				#				ctr_pos	<-this$.pre_dat[,as.character(this$.static_nem_design$.design[which(as.character(this$.static_nem_design$.design[,"group"])=="ctr_pos"),"sample"])]
				#				ctr_neg	<-this$.pre_dat[,as.character(this$.static_nem_design$.design[which(as.character(this$.static_nem_design$.design[,"group"])=="ctr_neg"),"sample"])]
				#				samp	<-this$.pre_dat[,as.character(this$.static_nem_design$.design[
				#										which(as.character(this$.static_nem_design$.design[,"group"])!="ctr_pos" & 
				#														as.character(this$.static_nem_design$.design[,"group"])!="ctr_neg"), "sample"])]
				#				up_reg <- which(rowMeans(ctr_pos) - rowMeans(ctr_neg) > log(n_fold))
				#				down_reg <- which(rowMeans(ctr_pos) - rowMeans(ctr_neg) < -log(n_fold))
				#				up_samp.sel <- samp[up_reg,]
				#				up_ctr_pos.sel <- ctr_pos[up_reg,]
				#				up_ctr_neg.sel <- ctr_neg[up_reg,]
				#				down_samp.sel <- samp[down_reg,]
				#				down_ctr_pos.sel <- ctr_pos[down_reg,]
				#				down_ctr_neg.sel <- ctr_neg[down_reg,]
				#				# count false decisions for different cutoff levels
				#				count.false.decisions <- function(x){
				#					up.thrsh    <- x*rowMeans(up_ctr_pos.sel) + (1-x)*rowMeans(up_ctr_neg.sel)
				#					up_ctr_pos.disc <- (up_ctr_pos.sel <= up.thrsh)*1
				#					up_ctr_neg.disc <- (up_ctr_neg.sel <= up.thrsh)*1        
				#					down.thrsh    <- x*rowMeans(down_ctr_neg.sel) + (1-x)*rowMeans(down_ctr_pos.sel)
				#					down_ctr_pos.disc <- (down_ctr_pos.sel >= down.thrsh)*1
				#					down_ctr_neg.disc <- (down_ctr_neg.sel >= down.thrsh)*1  
				#					alpha    <- round(sum(sum(up_ctr_pos.disc),sum(down_ctr_pos.disc))/(length(up_ctr_pos.disc)+length(down_ctr_pos.disc)),2)
				#					beta     <- 1-round(sum(sum(up_ctr_neg.disc),sum(down_ctr_pos.disc))/(length(up_ctr_neg.disc)+length(down_ctr_neg.disc)),2)                   
				#					return(c(alpha,beta))
				#				}
				#				false<-sapply(cutoff,count.false.decisions)
				#				dimnames(false) <- list(c("alpha","beta"),as.character(cutoff))
				#				#mycutoff <- cutoff[which.min(false[1,]+false[2,])]
				#				mycutoff <- cutoff[which.min(false[2,])]				
				#				# apply chosen cutoff
				#				up.thrsh    <- mycutoff*rowMeans(up_ctr_pos.sel) + (1-mycutoff)*rowMeans(up_ctr_neg.sel)
				#				up_ctr_pos.disc <- (up_ctr_pos.sel <= up.thrsh)*1
				#				up_ctr_neg.disc <- (up_ctr_neg.sel <= up.thrsh)*1        
				#				down.thrsh    <- mycutoff*rowMeans(down_ctr_neg.sel) + (1-mycutoff)*rowMeans(down_ctr_pos.sel)
				#				down_ctr_pos.disc <- (down_ctr_pos.sel >= down.thrsh)*1
				#				down_ctr_neg.disc <- (down_ctr_neg.sel >= down.thrsh)*1  
				#				this$.alpha    <- round(sum(sum(up_ctr_pos.disc),sum(down_ctr_pos.disc))/(length(up_ctr_pos.disc)+length(down_ctr_pos.disc)+p_counts),2)
				#				this$.beta     <- 1-round(sum(sum(up_ctr_neg.disc),sum(down_ctr_pos.disc))/(length(up_ctr_neg.disc)+length(down_ctr_neg.disc)+p_counts),2)   
				#				this$.pro_dat <- rbind((up_samp.sel <= up.thrsh)*1,(down_samp.sel >= down.thrsh)*1)
							}
						}
)
###Get error rates alpha and beta
setMethodS3(
		name		=	"GetAlpha",
		class		=	"StaticNEMData",
		definition	=	function(this, ...) {
							if(is.null(this$.alpha))
								stop("[StaticNEMData]: Input data has not been processed! \n")
							return(this$.alpha)
						}
)
setMethodS3(
		name		=	"GetBeta",
		class		=	"StaticNEMData",
		definition	=	function(this, ...) {
							if(is.null(this$.beta))
								stop("[StaticNEMData]: Input data has not been processed! \n")
							return(this$.beta)
						}
)
#Get input data
setMethodS3(
		name		=	"GetPreData",
		class		=	"StaticNEMData",
		definition	=	function(this, ...) {
							if(is.null(this$.pre_data)) 
								stop("[StaticNEMData]: No data input! \n")
							return(this$.pro_data)
						}
)

#Get processed data
setMethodS3(
		name		=	"GetProData",
		class		=	"StaticNEMData",
		definition	=	function(this, ...) {
							if(is.null(this$.pro_data)) 
								stop("[StaticNEMData]: Input data has not been processed! \n")
							return(this$.pro_data)
						}
)


















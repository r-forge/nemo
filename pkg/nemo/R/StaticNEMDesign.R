###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 21:07:59, on 17 Jul 2010
###############################################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Constructor~~~~~~~~~~~~~~~~~~~~~~~~~~~##
setConstructorS3(
		name		=	"StaticNEMDesign", 
		definition	=	function(design=NULL) {
			
							if(!is.null(design)) {
								design.colnames<-colnames(design)
								#check if input design table has enough information
								if(length(design.colnames)<3) 
									stop("[StaticNEMDesign]: Missing column/columns in design table!\n")
								#check if input design table has the correct format
								if(design.colnames[1]!="sample" | design.colnames[2]!="group" | design.colnames[3]!="replicate") 
									stop("[StaticNEMDesign]: Reorder columns in design table as: sample, group, replicate! \n") 
								#check if has control samples
								if(!("ctr" %in% design[,"group"]) & !("ctr_pos" %in% design[,"group"] & "ctr_neg" %in% design[,"group"]))
									stop("[StaticNEMDesign]: Missing control samples!\n")
							}
							extend(
									Object(),
									"StaticNEMDesign",
									.design=design
							)
						}
)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define Methods~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
###Summarize design table
setMethodS3(
		name		=	"DesignSummary",
		class		=	"StaticNEMDesign",
		definition	=	function(this,...) {
							n_sample<-nrow(this$.design)
							pert_rep<-data.frame(table(this$.design[,"group"]),stringsAsFactors=FALSE)
							cat("[StaticNEMDesign]: A summary of input design table \n")
							cat("[StaticNEMDesign]: Experimental groups and their replicate numbers: \n")
							cat("Group","Replicate No.",sep="\t")
							for(g in 1:nrow(pert_rep)) {
								cat("\n")
								cat(as.character(pert_rep[g,1]),as.integer(as.character(pert_rep[g,2])),sep="\t\t")				
							}
							cat("\n")
						}
)
###Get design table
setMethodS3(
		name		=	"GetDesignTable",
		class		=	"StaticNEMDesign",
		definition	=	function(this,...) {
							cat("[StaticNEMDesign]: Get design table \n")
							return(this$.design)
						}
)
###Get sample name by group
setMethodS3(
		name		=	"GetSampleByGroup",
		class		=	"StaticNEMDesign",
		definition	=	function(this,group_name,...) {
							#cat("StaticNEMDesign: Get sample names by a group name \n")
							samples<-NULL
							if(class(group_name)!="character")
								stop("[StaticNEMDesign]: Input a group name! \n")
							for(name in group_name)
								samples<-c(samples, as.character(this$.design[which(this$.design[,"group"]==name),"sample"]))
							if(length(samples)==0)
								stop("[StaticNEMDesign]: No such group! \n")
							else
								return(samples)
						}
)
###Get group name by sample
setMethodS3(
		name		=	"GetGroupBySample",
		class		=	"StaticNEMDesign",
		definition	=	function(this,sample_name,...) {
							#cat("StaticNEMDesign: Get sample names by a group name \n")
							groups<-NULL
							if(class(sample_name)!="character")
								stop("[StaticNEMDesign]: Input a sample name! \n")
							for(name in sample_name)
								groups<-c(groups,as.character(this$.design[which(this$.design[,"sample"]==name),"group"]))
							if(length(groups)==0)
								stop("[StaticNEMDesign]: No such sample! \n")
							else
								return(groups)
						}
)
###Save design table
setMethodS3(
		name		=	"Save",
		class		=	"StaticNEMDesign",
		definition	=	function(this, file_name, ...) {
							cat("[StaticNEMDesign]: Save design table \n")
							design<-this$.design
							save(design, file=file_name)
						}
)





















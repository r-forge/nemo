###############################################################################
# TODO: 
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 21:53:11, on 18 Jul 2010
###############################################################################
library(nemo)
data(BoutrosRNAiExpression)
design <- read.table(file.path(.path.package("nemo"),"data","BoutrosRNAiExpression.design"), header=TRUE,sep="\t")
cat("*Create StaticNEMData & do Discretization of input data \n") 
nem.data<-StaticNEMData(design=design,data=BoutrosRNAiExpression,process_type="discretize")
nem.data$Discretize(cutoff=0.7)				
cat("*Create StaticNEMPara \n")
nem.para<-StaticNEMPara(nem_infer_engine="SEARCH",bayesian_infer_engine="ML",data_type="DISC",egene_pos_infer_engine="MAP",graphic_engine="dyNet",
alpha=nem.data$GetAlpha(), beta=nem.data$GetBeta(),lambda=NULL, ss_prior=NULL, ss_prior_type=NULL,se_prior=NULL,greedy_initial=NULL,
triplet_thresh=NULL,transitive_close=TRUE,delta=NULL)					
cat("*Create StaticNEM with StaticNEMPara and StaticNEMData & Do inference \n")
static.nem<-StaticNEM(static_nem_para=nem.para,static_nem_data=nem.data)
static.nem$Infer()
cat("*Visualization!\n ")
static.nem$Show(what="FULLGraph")
cat("*Over!\n ")
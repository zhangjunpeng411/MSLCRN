############################################################
## Example scripts of MSLCRN for inferring module-specific #
## lncRNA-mRNA causal regulatory network in human cancer   #
############################################################

source("MSLCRN.R")
lncR <- 1:9704
mR <- 9705:27986
datacsv_GBM <- "GBM.csv"
datacsv_LSCC <- "LSCC.csv"
datacsv_OvCa <- "OvCa.csv"
datacsv_PrCa <- "PrCa.csv"
ExpData_GBM <- expDataTranspose(datacsv_GBM)
ExpData_LSCC <- expDataTranspose(datacsv_LSCC)
ExpData_OvCa <- expDataTranspose(datacsv_OvCa)
ExpData_PrCa <- expDataTranspose(datacsv_PrCa)

module_GBM <- moduleWGCNA(ExpData_GBM)
module_LSCC <- moduleWGCNA(ExpData_LSCC)
module_OvCa <- moduleWGCNA(ExpData_OvCa)
module_PrCa <- moduleWGCNA(ExpData_PrCa)

Causal_Score_GBM <- causalScore(ExpData_GBM, module_GBM, lncR, mR, num.cores = 6)
Causal_Score_LSCC <- causalScore(ExpData_LSCC, module_LSCC, lncR, mR, num.cores = 6)
Causal_Score_OvCa <- causalScore(ExpData_OvCa, module_OvCa, lncR, mR, num.cores = 6)
Causal_Score_PrCa <- causalScore(ExpData_PrCa, module_PrCa, lncR, mR, num.cores = 6)

# The cutoff of absolute value of causal effects can be set from 0.1 to 0.6 with step 0.05. For compromise, the cutoff is set to 0.45. 
Causal_lncRmR_GBM <- lapply( seq_along(Causal_Score_GBM), function(i) 
    cbind(colnames(Causal_Score_GBM[[i]])[which(abs(Causal_Score_GBM[[i]])>=0.45, arr.ind = TRUE)[, 2]], 
    rownames(Causal_Score_GBM[[i]])[which(abs(Causal_Score_GBM[[i]])>=0.45, arr.ind = TRUE)[, 1]]) )

Causal_lncRmR_LSCC <- lapply( seq_along(Causal_Score_LSCC), function(i) 
    cbind(colnames(Causal_Score_LSCC[[i]])[which(abs(Causal_Score_LSCC[[i]])>=0.45, arr.ind = TRUE)[, 2]], 
    rownames(Causal_Score_LSCC[[i]])[which(abs(Causal_Score_LSCC[[i]])>=0.45, arr.ind = TRUE)[, 1]]) )

Causal_lncRmR_OvCa <- lapply( seq_along(Causal_Score_OvCa), function(i) 
    cbind(colnames(Causal_Score_OvCa[[i]])[which(abs(Causal_Score_OvCa[[i]])>=0.45, arr.ind = TRUE)[, 2]], 
    rownames(Causal_Score_OvCa[[i]])[which(abs(Causal_Score_OvCa[[i]])>=0.45, arr.ind = TRUE)[, 1]]) )

Causal_lncRmR_PrCa <- lapply( seq_along(Causal_Score_PrCa), function(i) 
    cbind(colnames(Causal_Score_PrCa[[i]])[which(abs(Causal_Score_PrCa[[i]])>=0.45, arr.ind = TRUE)[, 2]], 
    rownames(Causal_Score_PrCa[[i]])[which(abs(Causal_Score_PrCa[[i]])>=0.45, arr.ind = TRUE)[, 1]]) )

## Calculate the number of positive and negative lncRNA-mRNA pairs
Causal_Score_GBM_Effect <- lapply( seq_along(Causal_Score_GBM), function(i) 
    Causal_Score_GBM[[i]][abs(Causal_Score_GBM[[i]])>=0.45] )
Global_Causal_Score_GBM_Effect <- unlist(Causal_Score_GBM_Effect)
GBM_Num_Negative <- length(which(Global_Causal_Score_GBM_Effect<0))
GBM_Num_Positive <- length(which(Global_Causal_Score_GBM_Effect>0))

Causal_Score_LSCC_Effect <- lapply( seq_along(Causal_Score_LSCC), function(i) 
    Causal_Score_LSCC[[i]][abs(Causal_Score_LSCC[[i]])>=0.45] )
Global_Causal_Score_LSCC_Effect <- unlist(Causal_Score_LSCC_Effect)
LSCC_Num_Negative <- length(which(Global_Causal_Score_LSCC_Effect<0))
LSCC_Num_Positive <- length(which(Global_Causal_Score_LSCC_Effect>0))

Causal_Score_OvCa_Effect <- lapply( seq_along(Causal_Score_OvCa), function(i) 
    Causal_Score_OvCa[[i]][abs(Causal_Score_OvCa[[i]])>=0.45] )
Global_Causal_Score_OvCa_Effect <- unlist(Causal_Score_OvCa_Effect)
OvCa_Num_Negative <- length(which(Global_Causal_Score_OvCa_Effect<0))
OvCa_Num_Positive <- length(which(Global_Causal_Score_OvCa_Effect>0))

Causal_Score_PrCa_Effect <- lapply( seq_along(Causal_Score_PrCa), function(i) 
    Causal_Score_PrCa[[i]][abs(Causal_Score_PrCa[[i]])>=0.45] )
Global_Causal_Score_PrCa_Effect <- unlist(Causal_Score_PrCa_Effect)
PrCa_Num_Negative <- length(which(Global_Causal_Score_PrCa_Effect<0))
PrCa_Num_Positive <- length(which(Global_Causal_Score_PrCa_Effect>0))

## Global lncRNA-mRNA causal relationships in GBM, LSCC, OvCa and PrCa
Global_Causal_lncRmR_GBM <- unique(do.call(rbind, Causal_lncRmR_GBM))
Global_Causal_lncRmR_LSCC <- unique(do.call(rbind, Causal_lncRmR_LSCC))
Global_Causal_lncRmR_OvCa <- unique(do.call(rbind, Causal_lncRmR_OvCa))
Global_Causal_lncRmR_PrCa <- unique(do.call(rbind, Causal_lncRmR_PrCa))

Global_Causal_lncRmR_GBM_graph <- make_graph(c(t(Global_Causal_lncRmR_GBM)), directed = TRUE)
Global_Causal_lncRmR_LSCC_graph <- make_graph(c(t(Global_Causal_lncRmR_LSCC)), directed = TRUE)
Global_Causal_lncRmR_OvCa_graph <- make_graph(c(t(Global_Causal_lncRmR_OvCa)), directed = TRUE)
Global_Causal_lncRmR_PrCa_graph <- make_graph(c(t(Global_Causal_lncRmR_PrCa)), directed = TRUE)

## Experimentally validated lncRNA-mRNA interactions from NPInter_3.0, LncRNADisease_2015 and LncRNA2Target_1.2 
Validated_lncRmR <- read.csv("NPInter_3.0+LncRNADisease_2015+LncRNA2Target_1.2.csv", header=TRUE, sep=",")
Validated_lncRmR_graph <- make_graph(c(t(Validated_lncRmR)), directed = TRUE)

Validated_GBM <- Global_Causal_lncRmR_GBM_graph %s% Validated_lncRmR_graph
Validated_LSCC <- Global_Causal_lncRmR_LSCC_graph %s% Validated_lncRmR_graph
Validated_OvCa <- Global_Causal_lncRmR_OvCa_graph %s% Validated_lncRmR_graph
Validated_PrCa <- Global_Causal_lncRmR_PrCa_graph %s% Validated_lncRmR_graph

## Gene list of module-specific lncRNA-mRNA causal relationships in GBM, LSCC, OvCa and PrCa
Causal_lncRmR_GBM_Modulelist <- lapply(seq_along(Causal_lncRmR_GBM), function(i) unique(c(Causal_lncRmR_GBM[[i]][, 1], Causal_lncRmR_GBM[[i]][, 2])))
Causal_lncRmR_LSCC_Modulelist <- lapply(seq_along(Causal_lncRmR_LSCC), function(i) unique(c(Causal_lncRmR_LSCC[[i]][, 1], Causal_lncRmR_LSCC[[i]][, 2])))
Causal_lncRmR_OvCa_Modulelist <- lapply(seq_along(Causal_lncRmR_OvCa), function(i) unique(c(Causal_lncRmR_OvCa[[i]][, 1], Causal_lncRmR_OvCa[[i]][, 2])))
Causal_lncRmR_PrCa_Modulelist <- lapply(seq_along(Causal_lncRmR_PrCa), function(i) unique(c(Causal_lncRmR_PrCa[[i]][, 1], Causal_lncRmR_PrCa[[i]][, 2])))

## Survival analysis
SurvData_GBM <- read.csv("GBM_survival.csv", header=TRUE, sep=",")
SurvData_LSCC <- read.csv("LSCC_survival.csv", header=TRUE, sep=",")
SurvData_OvCa <- read.csv("OvCa_survival.csv", header=TRUE, sep=",")
SurvData_PrCa <- read.csv("PrCa_survival.csv", header=TRUE, sep=",")
Surv_Analysis_GBM <- moduleSurvival(Causal_lncRmR_GBM_Modulelist, ExpData_GBM, SurvData_GBM)
Surv_Analysis_LSCC <- moduleSurvival(Causal_lncRmR_LSCC_Modulelist, ExpData_LSCC, SurvData_LSCC)
Surv_Analysis_OvCa <- moduleSurvival(Causal_lncRmR_OvCa_Modulelist, ExpData_OvCa, SurvData_OvCa)
Surv_Analysis_PrCa <- moduleSurvival(Causal_lncRmR_PrCa_Modulelist, ExpData_PrCa, SurvData_PrCa)

## GO and KEGG enrichment analysis
Enrich_Analysis_GBM <- moduleEA(Causal_lncRmR_GBM_Modulelist)
Enrich_Analysis_LSCC <- moduleEA(Causal_lncRmR_LSCC_Modulelist)
Enrich_Analysis_OvCa <- moduleEA(Causal_lncRmR_OvCa_Modulelist)
Enrich_Analysis_PrCa <- moduleEA(Causal_lncRmR_PrCa_Modulelist)

## Disease enrichment analysis
GBM_genes <- read.csv("GBM_genes.csv", header=FALSE, sep=",")
LSCC_genes <- read.csv("LSCC_genes.csv", header=FALSE, sep=",")
OvCa_genes <- read.csv("OvCa_genes.csv", header=FALSE, sep=",")
PrCa_genes <- read.csv("PrCa_genes.csv", header=FALSE, sep=",")
Disease_Analysis_GBM <- moduleDE(ExpData_GBM, GBM_genes, Causal_lncRmR_GBM_Modulelist)
Disease_Analysis_LSCC <- moduleDE(ExpData_LSCC, LSCC_genes, Causal_lncRmR_LSCC_Modulelist)
Disease_Analysis_OvCa <- moduleDE(ExpData_OvCa, OvCa_genes, Causal_lncRmR_OvCa_Modulelist)
Disease_Analysis_PrCa <- moduleDE(ExpData_PrCa, PrCa_genes, Causal_lncRmR_PrCa_Modulelist)

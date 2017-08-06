#################################################################################
## MSLCRN: a novel step-wise method for infer-ring module-specific lncRNA-mRNA  #
## causal regulatory network in human cancer                                    #
## Aug 6th, 2017, written by Junpeng Zhang                                      #
#################################################################################

## Load required R packages 
library(WGCNA)
library(flashClust)
library(ParallelPC)
library(bnlearn)
library(pcalg)
library(parallel)
library(survival)
library(igraph)
library(clusterProfiler)
library(DOSE)

## Transpose expression dataset (rows are genes, columns are samples) 
## into expression dataset (rows are samples, columns are genes).
expDataTranspose <- function(datacsv){
    
    ExpData <- read.csv(datacsv, header = FALSE, sep=",")
    ExpData <- t(ExpData)
    ExpData_Transpose <- matrix( as.numeric(ExpData[-1, -1]), nrow = nrow(ExpData)-1 )
    colnames(ExpData_Transpose) <- ExpData[, -1][1, ]
    rownames(ExpData_Transpose) <- ExpData[-1, ][, 1]
    return(ExpData_Transpose)
}

## Identification of co-expressed gene modules using WGCNA method
## ExpData: rows are samples, columns are genes.
moduleWGCNA <- function(ExpData, RsquaredCut = 0.9){
    
    Optimalpower <- pickSoftThreshold(ExpData, RsquaredCut = RsquaredCut)$powerEstimate
    adjacencymatrix <- adjacency(ExpData, power = Optimalpower)
    dissTOM <- TOMdist(adjacencymatrix)
    hierTOM <- flashClust(as.dist(dissTOM), method = "average")

    # The function cutreeDynamic colors each gene by the branches 
    # that result from choosing a particular height cut-off.
    # grey is reserved to color genes that are not part of any module.
    colorh <- cutreeDynamic(hierTOM, method="tree") + 1
    StandColor <- c("grey", standardColors(n = NULL))
    colorh <- unlist(lapply(1:length(colorh), function(i) StandColor[colorh[i]]))
    colorlevels <- unique(colorh)
    colorlevels <- colorlevels[-which(colorlevels=="grey")]

    Modulegenes <- lapply(1:length(colorlevels), function(i) 
        colnames(ExpData)[ which(colorh==colorlevels[i]) ])
    
    return(Modulegenes)
}

## By using ParallelPC R package, we generate lncRNA-mRNA causal relationships
## from matched lncRNA and mRNA expression data.
causalScore <- function(ExpData, Modulegenes, lncR, mR, 
    num.ModulelncRs = 2, num.ModulemRs = 2, pcmethod = "parallel", alpha = 0.01, num.cores = 2){
    
    ExpDataNames <- colnames(ExpData)
    CoExpData <- lapply( 1:length(Modulegenes), function(i) ExpData[, which(ExpDataNames %in% Modulegenes[[i]])] )
    cause <- lapply( 1:length(Modulegenes), function(i) which(colnames(CoExpData[[i]]) %in% ExpDataNames[lncR]) )
    effect <- lapply( 1:length(Modulegenes), function(i) which(colnames(CoExpData[[i]]) %in% ExpDataNames[mR]) )
    cause_length <- unlist( lapply( 1:length(Modulegenes), function(i) length(cause[[i]]) ) )
    effect_length <- unlist( lapply( 1:length(Modulegenes), function(i) length(effect[[i]]) ) )
    
    index <- which(cause_length>=num.ModulelncRs & effect_length>=num.ModulemRs)
    Update_CoExpData <- lapply(index, function(i) CoExpData[[i]])
    Update_cause <- lapply(index, function(i) cause[[i]])
    Update_effect <- lapply(index, function(i) effect[[i]])

    Causal_Score <- lapply( 1:length(Update_CoExpData), function(i) 
        IDA_parallel(Update_CoExpData[[i]], Update_cause[[i]], Update_effect[[i]], 
        pcmethod = pcmethod, alpha = alpha, num.cores = num.cores) )
    
    return(Causal_Score)
}

## Survival analysis of modules
moduleSurvival <- function(Modulelist, ExpData, SurvData, devidePercentage = 0.5, plot = FALSE) {

    names(Modulelist) <- seq_along(Modulelist)    
    Modulelist <- Modulelist[sapply(Modulelist,length)>0]
    ExpDataNames <- colnames(ExpData)    
    myfit <- list()
    LogRank <- list()    

    for (i in seq_along(Modulelist)) {
        Interin_Data <- cbind(SurvData[, seq(2, 3)], ExpData[, which(ExpDataNames %in% Modulelist[[i]])])
        Interin_Data <- na.omit(Interin_Data)

        try_mm <- try(coxph(Surv(time, status) ~ ., data = data.frame(Interin_Data)),
            silent = TRUE)
        if ("try-error" %in% class(try_mm))
            next

        mm <- coxph(Surv(time, status) ~ ., data = data.frame(Interin_Data))

        Risk_score <- predict(mm, newdata = data.frame(Interin_Data), type = "risk")

        group <- rep("NA", dim(Interin_Data)[1])
        group[Risk_score > quantile(Risk_score, probs = devidePercentage)] <- "High"
        group[Risk_score <= quantile(Risk_score, probs = devidePercentage)] <- "Low"

        Data <- cbind(Interin_Data[, seq_len(2)], group)
        myfit[[i]] <- survfit(Surv(time, status) ~ group, data = Data)

        sdf <- survdiff(Surv(time, status) ~ group, data = Data)
        sdf.p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
        HR <- (sdf$obs[1]/sdf$exp[1])/(sdf$obs[2]/sdf$exp[2])
        HRlow95 <- exp(log(HR) - qnorm(0.975) * sqrt(1/sdf$exp[1] + 1/sdf$exp[2]))
        HRup95 <- exp(log(HR) + qnorm(0.975) * sqrt(1/sdf$exp[1] + 1/sdf$exp[2]))

        LogRank[[i]] <- c(sdf$chisq, sdf.p.val, HR, HRlow95, HRup95)
    }

    if (plot) {
        for (i in seq_along(myfit)) {
	    if (!is.null(LogRank[[i]])) {
                dev.new()
                plot(myfit[[i]], lty = 1, col = c("red", "green"), main = paste("Module", i), xlab = "Time (Months)",
                    ylab = "Probability of survival")

                legend("topright", legend = c("High risk group", "Low risk group"), lty = seq_len(2),
                    col = c("red", "green"))
	    }
        }
    }

    LogRank_res <- do.call(rbind, LogRank)

    if (length(myfit) >= 1) {        
        colnames(LogRank_res) <- c("Chi-square", "p-value", "HR", "HRlow95", "HRup95")
	names(LogRank) <- names(Modulelist)
	LogRank[sapply(LogRank, is.null)] <- NULL
        rownames(LogRank_res) <- paste("Module", names(LogRank))
    }

    return(LogRank_res)
}

## Functional GO and KEGG enrichment analysis of modules
moduleEA <- function(Modulelist, ont = "BP", KEGGorganism = "hsa", OrgDb = "org.Hs.eg.db", 
    padjustvaluecutoff = 0.05, padjustedmethod = "BH") {
    
    names(Modulelist) <- seq_along(Modulelist)    
    Modulelist <- Modulelist[sapply(Modulelist,length)>0]

    entrezIDs <- lapply(seq_along(Modulelist), function(i) bitr(Modulelist[[i]], fromType = "SYMBOL",
        toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID)

    entrezIDs <- lapply(seq_along(Modulelist), function(i) as.character(entrezIDs[[i]]))
    names(entrezIDs) <- names(Modulelist)
    
    enrichGOs <- lapply(seq_along(Modulelist), function(i) enrichGO(entrezIDs[[i]], OrgDb = OrgDb,
        ont = ont, pvalueCutoff = padjustvaluecutoff, pAdjustMethod = padjustedmethod))
    names(enrichGOs) <- names(Modulelist)

    enrichKEGGs <- lapply(seq_along(Modulelist), function(i) enrichKEGG(entrezIDs[[i]], organism = KEGGorganism,
        pvalueCutoff = padjustvaluecutoff, pAdjustMethod = padjustedmethod))
    names(enrichKEGGs) <- names(Modulelist)

    return(list(enrichGOs, enrichKEGGs))
}

## Disease enrichment analysis using hypergeometric distribution test 
moduleDE <- function(ExpData, Diseasegenes, Modulelist) {

    names(Modulelist) <- seq_along(Modulelist)    
    Modulelist <- Modulelist[sapply(Modulelist,length)>0]

    B <- dim(ExpData)[2]
    N <- length(intersect(colnames(ExpData), as.matrix(Diseasegenes)))
    M <- unlist(lapply(seq_along(Modulelist), function(i) length(Modulelist[[i]])))
    x <- unlist(lapply(seq_along(Modulelist), function(i) length(intersect(Modulelist[[i]], as.matrix(Diseasegenes)))))    
    p.value <- 1 - phyper(x - 1, N, B - N, M)
    names(p.value) <- paste("Module", names(Modulelist))
    return(p.value)
}



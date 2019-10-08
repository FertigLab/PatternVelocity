#Compute PatternMarkerStatistics using CoGAPS functions

library(CoGAPS)

Gene_patterns <- read.csv('Geneweights10xHighVar.nP80.csv')

which(duplicated(Gene_patterns$X))
which(Gene_patterns$X == Gene_patterns[1519,1])
dif <- Gene_patterns[1518,-1] - Gene_patterns[1519,-1]
sum(abs(dif))
#The difference between the values is very low
Gene_patterns <- Gene_patterns[-1519,]
row.names(Gene_patterns) <- Gene_patterns$X
Gene_patterns$X   <- NULL


#From CoGAPS
patternMarkersCalculation <- function(mat, threshold, lp)
{
    nMeasurements <- nrow(mat)
    nPatterns <- ncol(mat)
    rowMax <- t(apply(mat, 1, function(x) x / max(x)))
    
    if (!is.na(lp))
    {
        if (length(lp) != nPatterns)
        {
            warning("lp length must equal the number of columns of the Amatrix")
        }
        sstat <- apply(rowMax, 1, function(x) sqrt(t(x-lp) %*% (x-lp)))
        ssranks <- rank(sstat)
        ssgenes.th <- names(sort(sstat,decreasing=FALSE,na.last=TRUE))
    }
    else
    {
        # determine which measurements are most associated with each pattern
        sstat <- matrix(NA, nrow=nMeasurements, ncol=nPatterns, dimnames=dimnames(mat))
        ssranks <- matrix(NA, nrow=nMeasurements, ncol=nPatterns, dimnames=dimnames(mat))
        ssgenes <- matrix(NA, nrow=nMeasurements, ncol=nPatterns, dimnames=NULL)
        for (i in 1:nPatterns)
        {
            lp <- rep(0,nPatterns)
            lp[i] <- 1
            sstat[,i] <- unlist(apply(rowMax, 1, function(x) sqrt(t(x-lp) %*% (x-lp))))
            ssranks[,i] <- rank(sstat[,i])
            ssgenes[,i] <- names(sort(sstat[,i], decreasing=FALSE, na.last=TRUE))
        }
        if (threshold=="cut")
        {
            geneThresh <- sapply(1:nPatterns, function(x) min(which(ssranks[ssgenes[,x],x] > apply(ssranks[ssgenes[,x],],1,min))))
            ssgenes.th <- sapply(1:nPatterns, function(x) ssgenes[1:geneThresh[x],x])
        }
        else if (threshold=="all")
        {
            pIndx <- unlist(apply(sstat, 1, which.min))
            gBYp <- list()
            for (i in sort(unique(pIndx)))
            {
                gBYp[[i]] <- sapply(strsplit(names(pIndx[pIndx==i]),"[.]"),function(x) x[[1]][1])
            }
            ssgenes.th <- lapply(1:max(sort(unique(pIndx))), function(x)
                ssgenes[which(ssgenes[,x] %in% gBYp[[x]]),x])
        }
        else
        {
            stop("invalid thresold argument")
        }
    }
    return(list("PatternMarkers"=ssgenes.th, "PatternRanks"=ssranks,
        "PatternMarkerScores"=sstat))       
}


patMarStat <- patternMarkersCalculation(Gene_patterns, threshold = "all",lp = NA)
patmarStatCut <- patternMarkersCalculation(Gene_patterns, threshold = "cut", lp = NA)


write.csv(file = 'patMarStatCutScores.csv',x = patmarStatCut$PatternMarkerScores)
write.csv(patmarStatCut$PatternRanks,'patMarStatCutRanks.csv')
write.csv(patMarStat$PatternMarkerScores,'patMarStatScores.csv')
write.csv(patMarStat$PatternRanks,'patMarStatRanks.csv')


patMarStatMarkers <- data.frame(pattern = 1:80,markerGenes = NA)
patMarStatCutMarkers <- data.frame(pattern = 1:80,markerGenes = NA)
for(i in 1:80){
    patMarStatMarkers$markerGenes[i] <- paste0(as.vector(patMarStat$PatternMarkers[[i]]),collapse = ',')
    patMarStatCutMarkers$markerGenes[i] <- paste0(as.vector(patmarStatCut$PatternMarkers[[i]]),collapse = ',')
    }


write.csv(file = 'patMarStatCutMarkers.csv', x = patMarStatCutMarkers)
write.csv(file = 'patMarStatMarkers.csv',x = patMarStatMarkers)




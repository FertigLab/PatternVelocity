# +
## Compute PatternMarkerStatistic Velocity 
## Args:
## velocity - Anndata object containing velocity information (gene x cell)
## pmarkScores - PatternMarkerScores (gene x pattern)
## pMat - CoGAPS Pattern matrix (pattern x cell)


#edited version of JHPCE commands ran by Gaurav
#Zeinab made many changes so that the code runs with data given
# -

require('data.table')
velocity <- fread('dataonlyVelocity_Horizontal.csv') 

row.names(velocity) = velocity$index
velocity$index = NULL

pmarkScores <- fread('patMarStatScores.csv')
rownames(pmarkScores) = pmarkScores$V1
pmarkScores$V1 = NULL

pMat <- fread('P_cogaps.csv') #May be a different P-matrix than the one Gaurav used (splitCnamesPmat?!)
rownames(pMat) = pMat$V1
pMat$V1 = NULL

# +
#Exctract data from genes common to the Velocity matrix and the PatternMarkerStatistic matrix
commonGenes = intersect( colnames(velocity), rownames(pmarkScores))

velocity_common_genes = subset(velocity, select=commonGenes)
rownames(velocity_common_genes) = rownames(velocity)
# -


pmarkScores_common = pmarkScores[rownames(pmarkScores) %in% commonGenes,]
rownames(pmarkScores_common) = rownames(pmarkScores)[rownames(pmarkScores) %in% commonGenes]


# +
#match pMat arcodes with the format of velocity barcodes
pMat_matched_colnames = gsub('\\.\\d+', '', colnames(pMat))
colnames(pMat) = pMat_matched_colnames


Index_fmt <- strsplit(rownames(velocity_common_genes),'x',fixed=T)
Index_fmt <- gsub(':','.',Index_fmt)
rownames(velocity_common_genes) = Index_fmt

# +
#Fetch data from common cells in Pattern matrix and Velocity matrix

commonCells = intersect(rownames(velocity_common_genes), colnames(pMat))

pMat_common = subset(pMat, select=commonCells)
velocity_common = velocity_common_genes[rownames(velocity_common_genes) %in% commonCells,]
rownames(velocity_common) = rownames(velocity_common_genes)[rownames(velocity_common_genes) %in% commonCells]
# -

dim(pMat_common)
dim(velocity_common)
dim(pmarkScores_common)

# +
#Zeinab's guess as to what Gaurav did:
veloPmsMat = as.matrix(velocity_common) %*% as.matrix(pmarkScores_common)
rownames(veloPmsMat) = rownames(velocity_common)

#####
veloPmsWgtMat <- veloPmsMat * as.matrix(t(pMat_common))
pmsSum = colSums(pmarkScores_common)
patVelo1 <- t(t(veloPmsWgtMat)/pmsSum)
#write.csv(veloPmsWgtMat,file='velocityPatmarkScoresPatWeightMatrix.csv')
#patVelo1 <- t(t(veloPmsWgtMat)/pmsSum)
#write.csv(patVelo1,file='patternVelocityIteration1.csv')

# +
write.csv(veloPmsMat,file='patternVelocity_Horizontal.csv')


# -








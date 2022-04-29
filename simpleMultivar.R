simpleMultivar  <- function(...) {


# load libraries

library(lattice)
library(APAFstats)
library(nlme)
# library(gplots)


PCA <- function (data, labelValue, scaleR = FALSE, scaleC = TRUE, k = min(dim(data)) - 
    1) 
{
    if (k > min(dim(data) - 1)) 
        Warning("The number of components was too large compared to the data and was adjusted accordingly")
    k <- min(k, min(dim(data)) - 1)
    if (scaleR) {
        row.nrm <- apply(data, 1, sd)
        row.nrm <- pmax(row.nrm, 1e-04)
        data <- sweep(data, 1, row.nrm, FUN = "/")
    }
    result <- try(prcomp(data, retx = TRUE, scale = scaleC), 
        silent = TRUE)
    if (inherits(result, "try-error")) 
        Error("Failed to Calculate Principal Components")
    componentVariances <- result$sdev^2
    componentLoadings <- result$rotation[, 1:k]
    componentScores <- result$x[, 1:k]
    totalVariance <- sum(componentVariances)
    componentVariances <- componentVariances[1:k]
    z <- componentScores
    plot(cloud(z[, 1] ~ z[, 3] + z[, 2], groups = as.factor(labelValue), 
        auto.key = list(points = TRUE, pch = 19, space = "right"), 
        xlab = "PC 3", ylab = "PC 2", zlab = "PC 1", distance = 0.1, 
        main = "Projection in the space of the first 3 princial components"))
    value <- list(componentVariances = componentVariances, componentScores = componentScores, 
        componentLoadings = componentLoadings, summary = summary(result))
    value
}


############
# parameters
############

# initialize to defaults
colscheme <- "CyanMagenta"
logData <- FALSE
ignoreCols <- 4
imputeMissing <- FALSE


args <- list(...)

for(i in 1:length(args)) {
flag <- substring(args[[i]], 0, 2)
value <- substring(args[[i]], 3, nchar(args[[i]]))


if(flag=='-f') fname <- value;
if(flag=='-i') ignoreCols <- as.numeric(value);
if(flag=='-l') logData <- value;
if(flag=='-c') colscheme <- value;
if(flag=='-z') imputeMissing <- value;


} 


logData <-  (logData == "yes") 
imputeMissing <- (imputeMissing == "yes")

	





data <- read.csv(fname)
dat <- data[, -c(1:ignoreCols)]
rownames(dat) <- paste(1:nrow(data), data[,ignoreCols])
samples <- colnames(dat)
Group <- as.factor( gsub("(.*)[-\\.](.*)", "\\1", samples))



if ( colscheme == "CyanMagenta" ) {
	colkey = cm.colors(256)
} else {
	colkey = redgreen(256)[256:1]
}


############
# levelplot
############



png("Levelplot.png", 2000, 2000, res=200)
print(levelplot(as.matrix(t(dat[nrow(dat):1,])), scales=list(x=list(rot=45)) ), title="Data Levelplot")
dev.off()




##############################
# 1. Heatmap - log and raw
#	
###############################

if (imputeMissing) {
	hdat <- dat 
	hdat[is.na(hdat)] <- 0
}


hdat <- as.matrix(na.omit(dat))
if (logData) hdat <- log( hdat + 0.01 )


png("HeatmapEuclidean.png", 2000, 2000, res=300)
h.res <- heatmap(hdat, margins=c(10,20) , 
    col = colkey, 
	ColSideColors = rainbow(nlevels(Group))[Group],
	cexRow=.5, cexCol=1.1)
legend("topright", fill=rainbow(nlevels(Group)), legend=levels(Group))

dev.off()

png("HClustEuclidean.png", 2000, 2000, res=300)
HClust(t(hdat), method="complete", clabel=Group,  metric="euclidean")
dev.off()


ordered.mat <- hdat[h.res$rowInd, h.res$colInd][nrow(hdat):1,]
write.csv(ordered.mat, file="HeatmapEuclideanData.csv")



cordist <- function(x) { as.dist((1-cor(t(x)))/2) }


png("HeatmapCordist.png", 2000, 2000, res=300)

hcor.res <- heatmap(hdat, margins=c(10,20) , 
    col = colkey,   dist=cordist,
	ColSideColors = rainbow(nlevels(Group))[Group],
	cexRow=.5, cexCol=1.1)
legend("topright", fill=rainbow(nlevels(Group)), legend=levels(Group))

dev.off()


png("HClustCordist.png", 2000, 2000, res=300)
HClust(t(hdat), method="complete", clabel=Group,  metric="pearsonCorrelation")
dev.off()


ordered.cor <- hdat[hcor.res$rowInd, hcor.res$colInd][nrow(na.omit(dat)):1,]
write.csv(ordered.cor, file="HeatmapCorrelationData.csv")
levelplot(t(hdat[h.res$rowInd, h.res$colInd]), scales=list(x=list(rot=45)))




#########
# 2. PCA
#########



png("PCA.png")
pca.res <- try( PCA(t(na.omit(log(dat+0.1))), Group, k=5, scaleC=FALSE) )
dev.off()

if (!inherits(pca.res, "try-error")) {

	ld <- pca.res$componentLoadings[,1:3]
	scores <- pca.res$componentScores[,1:3]
	write.csv(data.frame(ld, hdat), file="Loadings.csv")
	write.csv(scores, file="Scores.csv")

z <- pca.res$componentScores

png("PCA 2d - all.png", 2000, 2000, res=300)

layout(matrix(1:4, nrow=2))
plot(z[,1], z[,2], col=rainbow(nlevels(as.factor(Group)))[as.factor(Group)], pch=20, xlab="PC1", ylab="PC2")
text(z[,1], z[,2], samples, pos=3, cex=.5)
plot(z[,1], z[,3], col=rainbow(nlevels(as.factor(Group)))[as.factor(Group)], pch=20, xlab="PC1", ylab="PC3")
text(z[,1], z[,3], samples, pos=3, cex=.5)
plot(z[,2], z[,3], col=rainbow(nlevels(as.factor(Group)))[as.factor(Group)], pch=20, xlab="PC2", ylab="PC3")
text(z[,2], z[,3], samples, pos=3, cex=.5)
boxplot(log(dat+0.1),par(las=2), main="Boxplots of log data")

dev.off()


}



#######################################
# Grouped data - just see values better
#######################################



Values <- as.vector(as.matrix(dat))
Variables <- paste(1:nrow(data), data[,ignoreCols])

dlong <- data.frame(Values, Variables=rep(Variables, ncol(dat)), 
	 Group=rep(Group, each=nrow(dat)), Individual=rep(colnames(dat), each=nrow(dat)) )


gd <- groupedData(log( Values + 0.1) ~ Group | Variables, data=dlong)


png("Overall log Data.png", 2000, 2000, res=200) 
print(plot(gd, group=Group))
dev.off()


}




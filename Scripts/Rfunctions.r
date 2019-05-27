# filter genes with too many missing data and/or zero variance
geneFilter <- function(exprData)
{
  gsg <- goodSamplesGenes(exprData, verbose = 3)
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
    {printFlush(paste("Removing genes:", paste(names(expr6C)[!gsg$goodGenes], collapse = ", ")))}
    if (sum(!gsg$goodSamples)>0)
    {printFlush(paste("Removing samples:", paste(rownames(expr6C)[!gsg$goodSamples], collapse = ", ")))}
    # Remove the offending genes and samples from the data:
    return(exprData[gsg$goodSamples, gsg$goodGenes])
  }
  else
  {return(exprData)}
}

# TODO: filter genes with too many missing data and/or zero variance for multiple datasets
geneFilterMS <- function(exprData)
{
  gsg <- goodSamplesGenesMS(exprData, verbose = 3)
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
    {printFlush(paste("Removing genes:", paste(names(expr6C)[!gsg$goodGenes], collapse = ", ")))}
    if (sum(!gsg$goodSamples)>0)
    {printFlush(paste("Removing samples:", paste(rownames(expr6C)[!gsg$goodSamples], collapse = ", ")))}
    # Remove the offending genes and samples from the data:
    return(exprData[gsg$goodSamples, gsg$goodGenes])
  }
  else
  {return(exprData)}
}




sampleTree <- function(exprData, label = sample)
{
  rownames(exprData) <- samples[match(rownames(exprData), samples$sample), label]
  sampleTree = hclust(dist(exprData), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  # clusters according accession
}

plotSoftThresholdChoices <- function(exprData, maxSoftThrs = 20, title)
{
  # Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 12, to=maxSoftThrs, by=2))
  # Call the network topology analysis function
  sft <- pickSoftThreshold(exprData, powerVector = powers, verbose = 5)
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence", title));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity", title))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

# Module Connectivity calculation
summarizeConnectivity <- function(expr, modules)
{
  moduleConnectivity <- data.frame('Module name'= character(), 'KWithin'= numeric(), 'KOut'= numeric(), 'Quality module'= numeric())
  for (module in modules)
  {
    geneIdx <- colnames(expr$data)[expr$mergedColors == module]
    Kwithin <- mean(expr$geneConnectivity$kWithin[rownames(expr$geneConnectivity) %in% geneIdx])
    KOut <- mean(expr$geneConnectivity$kOut[rownames(expr$geneConnectivity) %in% geneIdx])
    Module.Quality <- Kwithin/KOut
    lineiwant <- data.frame('Module name'= module, 'KWithin'= Kwithin, 'KOut'= KOut, 'Quality module'= Module.Quality)
    moduleConnectivity <- rbind(moduleConnectivity, lineiwant)
  }
  return(moduleConnectivity)
}


#create color-coded table of the intersection counts
# Truncate p values smaller than 10^{-50} to 10^{-50}
Heatmap_overlap<- function(pTable, CountTbl, expr) {
  pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
  pTable[pTable>50 ] = 50 
  # Marginal counts (really module sizes)
  ModTotal<- list(apply(CountTbl, 1, sum), apply(CountTbl, 2, sum))
  # Actual plotting
  sizeGrWindow(10,7)
  #pdf(file = "/Volumes/nordborg/pub/forPieter/WGCNA/Results/Overlap modules 6 vs16.pdf", wi = 10, he = 7);
  par(mfrow=c(1,1));
  par(cex = 1.0);
  par(mar=c(8, 10.4, 2.7, 1)+0.3)
  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  labeledHeatmap(Matrix = pTable,
                 xLabels = paste(" ", expr[[2]]$Modules),
                 yLabels = paste(" ", expr[[1]]$Modules),
                 colorLabels = TRUE,
                 xSymbols = paste("16 ", expr[[2]]$Modules, ": ", ModTotal[[2]], sep=""),
                 ySymbols = paste("6 ", expr[[1]]$Modules, ": ", ModTotal[[1]], sep=""),
                 textMatrix = CountTbl,
                 colors = greenWhiteRed(100)[50:100],
                 main = "Correspondence of 6 set-specific and 16 set-specific modules",
                 cex.text = 0.5, cex.lab = 0.5, setStdMargins = FALSE);
}

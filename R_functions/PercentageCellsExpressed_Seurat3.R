#Pull genes of interest for cca
# This function calculates percentage of cells expressing a certain gene 
# and returns a dataframe containing gene names, fraction of cells and number of cells expressing each gene.
# SerautObj = a Seuraut 3 object
# genes.use = a vector of genes of interest
# ident.use = optional; a vector of clusters of interest;
#             default is NULL, which uses all cells for calculation
# group.by = optional; default is NULL;
#            Set group.by = "time" to calculate for each individual data set separately in a combined cca object 
#
#
#  Connor
# addded threshold arg (th)
# added bin
# made groupby work
#
# future ref: (bpcsubset@assays$RNA@counts["VSX1", ] > 0) & (bpcsubset@assays$RNA@counts["OTX2", ] > 0) for two genes





cellsExpressed <- function(SerautObj, genes.use, ident.use=NULL, group.by=NULL, th = 0, bin = T, v = F) {
  ### Initialize values 
  cellsExpressed.out <- c()
  percentage.out <- c()
  condition.out <- c()
  sd.out <- c()
  
  # ident.out <- rep(ident, length(genes.use))
  if (is.null(group.by)) {
    for (gene in genes.use) {
      ##  Index of row of gene
      #print(rownames(SerautObj@data))
      index <- grep(paste0("^", gene, "$"), rownames(SerautObj@assays$RNA@data))
      ## If 
      if (length(ident.use) == 0) {
        ## Pull out Vector of epression of said gene
        geneVector <- SerautObj@assays$RNA@data[index,]
      } else if (length(ident.use) > 0) {
        geneVector <- SerautObj@assays$RNA@data[index,grep(paste(ident.use, collapse ="|"), Idents(SerautObj))]
        if (v == T){
          print(length(geneVector))
          }
      }
      
      if (bin == T){
        ## Number of cells whos expression is above th
        cellsExpressed <- sum(geneVector > th)
        if (v == T){
          message(sprintf("%s Cells Expressed: %s", gene, cellsExpressed))
        }
      }
      else if(bin == F){
        #average expression
        cellsExpressed <- sum(geneVector)
        std <- sd(geneVector)
        sd.out <- c(sd.out, std)
        
        if (v == T){
          message(sprintf("%s Cells Expressed: %s", gene, cellsExpressed))
        }
      }
      else{
        message("Bin must be TRUE or FALSE")
        return(-1)
        }
      
      
      ## Fraction of cells expressed
      percentage <- cellsExpressed / length(geneVector)
      if (v == T){
        message(sprintf("Fraction %s Cells Expressed: %s", gene, percentage))
      }
      
      percentage.out <- c(percentage.out, percentage)
      cellsExpressed.out <- c(cellsExpressed.out, cellsExpressed)
    }
    if (bin == F){
      dfOUT <- data.frame(gene=genes.use, Fraction=percentage.out, Condition=condition.out, cellsExpressed=cellsExpressed.out, stdev = sd.out)
    }
    else{
      dfOUT<- data.frame(gene=genes.use, Fraction=percentage.out, cellsExpressed=cellsExpressed.out)
    }
  } 
  
  else if (!is.list(group.by)){
    groupingVars <- unique(SerautObj@meta.data[[group.by]])
    for (condition in groupingVars) {
      cellsOfCondition <- grep(condition, SerautObj@meta.data[[group.by]])
      for (gene in genes.use) {
        cellsOfIdent <- grep(paste(ident.use, collapse ="|"), Idents(SerautObj))
        cellsOfIntrest <- intersect(cellsOfIdent, cellsOfCondition)
        ##  Index of row of gene
        #print(rownames(SerautObj@data))
        index <- grep(paste0("^", gene, "$"), rownames(SerautObj@assays$RNA@data))
        ## If 
        if (length(ident.use) == 0) {
          ## Pull out Vector of epression of said gene
          geneVector <- SerautObj@assays$RNA@data[index,cellsOfIntrest]
        } else if (length(ident.use) > 0) {
          geneVector <- SerautObj@assays$RNA@data[index,cellsOfIntrest]
          if (v == T){
            print(length(geneVector))
          }
        }
        
        
        if (bin == T){
          ## Number of cells whos expression is above th
          cellsExpressed <- sum(geneVector > th)
          if (v == T){
            message(sprintf("%s Cells Expressed: %s", gene, cellsExpressed))
          }
        }
        else if(bin == F){
          #average expression
          cellsExpressed <- sum(geneVector)
          
          std <- sd(geneVector)
          sd.out <- c(sd.out, std)
          
          if (v == T){
            message(sprintf("%s Sum Expression: %s", gene, cellsExpressed))
            message(sprintf("%s STDEV Expression: %s", gene, std))
          }
        }
        else{
          message("Bin must be TRUE or FALSE")
          return(-1)
        }
        
        
        ## Fraction of cells expressed
        percentage <- cellsExpressed / length(geneVector)
        if (v == T){
          message(sprintf("Fraction %s Cells Expressed: %s", gene, percentage))
        }
        
        percentage.out <- c(percentage.out, percentage)
        cellsExpressed.out <- c(cellsExpressed.out, cellsExpressed)
        condition.out <-c(condition.out, condition)
      }
    }
    print(groupingVars)
    if (bin == F){
      dfOUT <- data.frame(gene=genes.use, Fraction=percentage.out, Condition=condition.out, cellsExpressed=cellsExpressed.out, stdev = sd.out)
    }
    else{
      dfOUT <- data.frame(gene=genes.use, Fraction=percentage.out, Condition=condition.out, cellsExpressed=cellsExpressed.out)
    }
  }
  else{
    message("this function does not work for lists at this time")
    return(-1)
  }
  return(dfOUT)
}




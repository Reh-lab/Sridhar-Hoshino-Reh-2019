library(Seurat)
 
dev_avg <- function(Se = org_Se, dev_genes = c("GAP43", "POU4F2", "POU4F1", "POU4F3", "NEFL", "NEFM", "ATOH7"),
                    Ctype = "RGC1", bin = FALSE, th = .0001, V = FALSE, Verb = FALSE){ 
  
  "This function is designed to be a manual version of pseudotime. It takes the following information:
  [1] a list of genes
  [2] a celltype
  [3] a 'treatment' column (ex time)
  [4] a Seurat object to operate on
  
  
  Limitations:
  
  R does interesting things with the $ and @ syntax which makes it hard to subset dynamically, this leads to the following     issues:
  [1] the cell type column must be called cell.type
  [2] the treatment column must be called time
  
  Additionally this code does not yet seperate by two treatment feilds (i.e. organoid and day) but that will happen soon
  "
  
  "Subset by cell type"
  
  # first we zoom in on the cell type
  sub = subset(Se, subset = cell.type ==  cType)
  
  if (V == TRUE){
    print("subset info (type)")
  }
  
  
  "Creat output df"
  df <- data.frame(matrix(ncol = length(unique(Se@meta.data$'time')),
                          nrow = length(dev_genes)))
  
  rownames(df) <- dev_genes
  names(df) <- unique(Se@meta.data$'time')
  
  if (V == TRUE){
    print("Dataframe stucture is:")
    print(df)
  }
  
  
  
  "Iterate over different timepoints"
  
  for(i in unique(Se@meta.data$time)){
    
    if (V == TRUE){
      print(i)
      #print(Se)
    }
    
    
    
    "To subset dynamically idents seems to need to be set to the column of intrest"
    Idents(sub) <- sub$time
    subtr <- subset(sub, idents = i)
    
    
    if (V == TRUE){
      print("subset info (type + time)")
      print(subtr)
    }
    
    
    
    "Now that we have subsetted by time and celltype we can get population size
    to normalize for population size later"
    num_cells = dim(subtr@assays$RNA@data)[2]
    
    
    if (V == TRUE){
      print("number of cells:")
      print(num_cells)
    }
    
    
    
    
    
    "Initialize avg_lv for first run, and clear avg_lv for later runs"
    avg_lv <- c()
    
    if (V ==TRUE){
      print("avg_lv cleared?")
      print(avg_lv)
    }
    
    
    "Iterate over the genes"
    for(gene in dev_genes){
      
      "see gene"
      if (V == TRUE){
        print(gene)
      }
      
      
      
      
      "Are we binning or are we averaging"
      if (bin == TRUE){
        #get the cells where there is expression above a threshold
        tot = sum(subtr@assays$RNA@data[gene,] > th)
        print(tot)
      }
      
      else if (bin == FALSE){
        tot = sum(subtr@assays$RNA@data[gene, ])
      }
      
      else{
        print("Error: please set bin to TRUE or FALSE")
        return()
      }
      
      
      
      
      "Normalize
      link normalized value and gene
      add to list
      "
      avg_ex = tot/num_cells
      names(avg_ex) <- gene
      avg_lv <- c(avg_lv, avg_ex)
      
    }
    
    "To view the list with values"
    if ((V == TRUE) | (Verb == TRUE)){
      print(i)
      print(avg_lv)
    }
    
    
    "load values to df"
    df[, i] <- avg_lv
    
    
  }
  
  return(df)
  
}
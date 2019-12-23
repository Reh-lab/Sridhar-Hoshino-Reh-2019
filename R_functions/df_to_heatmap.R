#functions for taking a dataframe of the folowing stucture
#
# gene  condition value
#
# and converting to the stucture
#
#       conditions... ... ...
# genes val val val ...


row_scale <- function(x){
  x <- x/sum(x)
}

gene_frame_to_heatmap <- function(df, sc = F, strip_cols = T){
  df <- reshape(df, timevar = "Condition", idvar = "gene", dir = "wide")
  rownames(df) <- df$gene
  df <- df[,-1]
  
  if (strip_cols == T){
    new_col= c()
    
    for (col in colnames(df)){
      #print(col)
      new_col <- c(new_col, substr(col, 10, 40))
    }
    
    #print(new_col)
    
    colnames(df) <- new_col
    }
  
  
  if(sc == T){
    df = t(apply(df, 1, row_scale))
  }
  
  return(df)
}